#!/usr/bin/env python3
"""Reconstruct every saved V10 diagnostic with a separate parser and Decimal arithmetic."""
import argparse
import copy
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys

from reference_measurement_union_sources import bind, verify


MAX_DOUBLE = Decimal.from_float(sys.float_info.max)
HEADER = dict(fields=['iter','prior','likelihood','posterior'],nested=True,format='MCON',
              version='0.2',projectScalarSchema='native-cjson-explicit-special-values-v6')


def load(text):
    def pairs(items):
        out = {}
        for key,value in items:
            assert key not in out
            out[key] = value
        return out
    def invalid(value):
        raise ValueError('Unencoded nonfinite: '+value)
    return json.loads(text,parse_float=Decimal,object_pairs_hook=pairs,parse_constant=invalid)


def number(value):
    assert type(value) in (int,Decimal)
    value = Decimal(value)
    assert value.is_finite()
    return value


def close(observed, expected, relative=Decimal('2e-13'), absolute=Decimal('2e-13')):
    observed = number(observed)
    assert abs(observed-expected) <= max(absolute,relative*abs(expected))


def check(row,mu,scale,reference=None):
    assert set(row)=={'iter','statistics//','parameters//','numericParameterQuality//'}
    assert type(row['iter']) is int and row['iter']>=0
    stats = row['statistics//'];quality=stats['__project_scalar_v6_quality__']
    values = {k:v for k,v in stats.items() if k!='__project_scalar_v6_quality__'}
    assert not values or set(values)=={'prior','likelihood','posterior'}
    for value in values.values():number(value)
    assert quality==dict(numericLeafCount=len(values),nonfinite=[],literalNullPaths=[])
    assert set(row['parameters//'])=={'S1/'}
    state = row['parameters//']['S1/']
    assert set(state)=={'latentLogAlpha','derivedAlpha','categoryRates','laplaceLocation','laplaceScale','latentLogDensity'}
    x=number(state['latentLogAlpha'])
    assert number(state['laplaceLocation'])==mu and number(state['laplaceScale'])==scale
    with localcontext() as ctx:
        ctx.prec=90
        expected=x.exp()
        overflow=expected>MAX_DOUBLE
        density=-(Decimal(2)*scale).ln()-abs(x-mu)/scale
        close(state['latentLogDensity'],density)
        if overflow:
            assert state['derivedAlpha']=='__project_scalar_v6__:positive_infinity'
            annotations=[dict(path=['S1/','derivedAlpha'],kind='positive_infinity')]
        else:
            close(state['derivedAlpha'],expected,absolute=Decimal(0))
            annotations=[]
        assert row['numericParameterQuality//']==dict(numericLeafCount=9,nonfinite=annotations,literalNullPaths=[])
        rates=state['categoryRates'];assert len(rates)==4
        rates=[number(rate) for rate in rates];assert min(rates)>=0
        close(sum(rates)/4,Decimal(1),relative=Decimal('1e-10'),absolute=Decimal('1e-10'))
        if x>=100:assert rates==[Decimal(1)]*4
    if reference is not None:
        assert row['iter']==reference['iter'] and row['statistics//']==reference['statistics//']
        assert state['derivedAlpha']==reference['parameters//']['S1/']['ASRV.Gamma:alpha']
    return overflow


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--receipt',type=Path,required=True)
    args=p.parse_args();assert not args.receipt.exists()
    qp=Path('metadata/baliphy_log_alpha_v10_software_validation_20261005_v1.json')
    tp=Path('metadata/baliphy_log_alpha_v10_software_transport_20261005_v1.json')
    qualified,transport=[json.loads(path.read_text()) for path in [qp,tp]]
    assert qualified['status']=='passed_reversible_latent_log_alpha_v10_full_source_and_paired_native_controls'
    assert transport['original_tool_terminal_exit_code']==0
    verify(transport['source_hashes']);pins=dict(transport['source_hashes'])
    all_rows=[];prior_rows=0;byte_checks=0
    for outcome in qualified['outcomes']:
        current=Path(outcome['native_receipt']).parent/'independent-chain-1'
        previous=Path('data/software_audits/baliphy-joint-fasta-v7-20261004-v1/native')/outcome['prior']/'attempt-0001/independent-chain-1'
        mu,scale={'broad':(Decimal(0),Decimal(2)),'centered':(Decimal(0),Decimal(1)),'package':(Decimal(6),Decimal(2))}[outcome['prior']]
        for name in outcome['byte_identical_files']:
            assert (current/name).read_bytes()==(previous/name).read_bytes();byte_checks+=1
        rows=[load(line) for line in (current/'C1.P1.log-alpha-samples.jsonl').read_text().splitlines()]
        original=[load(line) for line in (current/'C1.log.json').read_text().splitlines()]
        assert rows[0]==original[0]==HEADER
        rows,original=rows[1:],original[1:]
        assert len(rows)==len(original)==21 and [r['iter'] for r in rows]==list(range(21))
        for row,reference in zip(rows,original):assert not check(row,mu,scale,reference)
        all_rows.extend(rows);prior_rows+=len(rows)
    native=Path('data/software_audits/baliphy-log-alpha-v10-20261005-v1/deterministic/attempt-0001')
    deterministic=[load(line) for line in (native/'stdout.log').read_text().splitlines()]
    assert len(deterministic)==57 and [r['iter'] for r in deterministic]==list(range(57))
    grid=list(map(Decimal,['-20','-10','-6','-3','0','6','20','50','100','500','700','709','709.7','709.78','709.79','710','750','1000','5000']))
    overflows=0
    for start,mu,scale in [(0,Decimal(0),Decimal(2)),(19,Decimal(0),Decimal(1)),(38,Decimal(6),Decimal(2))]:
        subset=deterministic[start:start+19]
        assert [row['parameters//']['S1/']['latentLogAlpha'] for row in subset]==grid
        overflows+=sum(check(row,mu,scale) for row in subset)
    assert prior_rows==63 and byte_checks==18 and overflows==15
    # Real saved records are checked after single-field semantic changes;
    # expected failures are separate from immutable native output files.
    good=all_rows[0];special=deterministic[14];controls=[]
    mutations={
        'latent_value':lambda r:r['parameters//']['S1/'].update(latentLogAlpha=number(r['parameters//']['S1/']['latentLogAlpha'])+1),
        'derived_value':lambda r:r['parameters//']['S1/'].update(derivedAlpha=Decimal(2)),
        'prior_location':lambda r:r['parameters//']['S1/'].update(laplaceLocation=Decimal(1)),
        'prior_scale':lambda r:r['parameters//']['S1/'].update(laplaceScale=Decimal(3)),
        'latent_density':lambda r:r['parameters//']['S1/'].update(latentLogDensity=Decimal(0)),
        'category_rates':lambda r:r['parameters//']['S1/']['categoryRates'].__setitem__(0,Decimal(4)),
        'leaf_count':lambda r:r['numericParameterQuality//'].update(numericLeafCount=8),
        'finite_as_boolean':lambda r:r['parameters//']['S1/'].update(latentLogAlpha=True)}
    for name,mutate in mutations.items():
        bad=copy.deepcopy(good);mutate(bad)
        try:check(bad,Decimal(0),Decimal(2))
        except (AssertionError,ValueError):controls.append('reject_'+name)
        else:raise AssertionError('Mutation accepted: '+name)
    for name,mutate in [('missing_overflow_tag',lambda r:r['numericParameterQuality//'].update(nonfinite=[])),
                       ('wrong_overflow_state',lambda r:r['parameters//']['S1/'].update(derivedAlpha='__project_scalar_v6__:negative_infinity'))]:
        bad=copy.deepcopy(special);mutate(bad)
        try:check(bad,Decimal(0),Decimal(2))
        except (AssertionError,ValueError):controls.append('reject_'+name)
        else:raise AssertionError('Mutation accepted: '+name)
    for path in [qp,tp,Path(__file__)]:bind(pins,path)
    verify(pins)
    result=dict(status='passed_independent_V10_latent_alpha_and_native_output_reconstruction',
        checked_utc=datetime.now(timezone.utc).isoformat(),native_prior_rows=63,deterministic_rows=57,
        overflow_rows_retained=15,original_file_byte_checks=18,decimal_precision=90,
        semantic_rejection_controls=controls,source_hashes=pins,scientific_eligibility=False,
        posterior_qualified=False,original_chain_overflow_cause_proven=False,
        scope='Separate strict JSON/Decimal parser and90digit arithmetic check every saved latent draw, '
              'derived alpha, prior density/category state/quality tag and original paired bytes. Ten '
              'semantic mutations rejected. Shared original sources, SHA functions and native outputs '
              'remain explicit. This does not independently compile all405sources or qualify a full '
              'inference horizon, installed formatter repair, old review arrays or ancestral posterior.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
