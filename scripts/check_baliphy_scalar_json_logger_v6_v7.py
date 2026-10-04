#!/usr/bin/env python3
"""Qualify reversible V6 scalar rendering across the complete source grid.

Six paired native model runs compare the same seeds/settings before/after the
logger change. A separate native Value probe covers finite and nonfinite fields.
Every run starts in a new attempt root; no existing failed attempt is restarted.
"""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import run_attempt, sha
from baliphy_scalar_json_logger_v6c import HELPER, SCHEMA, QUALITY_KEY, transform, restore, validate_record
from baliphy_joint_node_logger_v5 import fresh_seeds, validate_frame
from independent_native_ancestral_alignment import fasta_records, tree_labels
from independent_joint_ancestral_frames import decode, write_arrays, verify_arrays
from read_baliphy_scalar_json_v6b import load, read_record, compare_tsv
from reference_measurement_union_sources import verify


FIXTURE_SEEDS = {20269001,20269002,20269003}
SEED_NAMESPACE = 'fungal-scalar-cjson-explicit-nonfinite-future-20261004-v6'


def seed_census(bindings):
    used = set()
    paths = ['metadata/baliphy_reference_sampler_qualification_plan_20261003.json',
             'metadata/baliphy_joint_node_logger_future_models_20261003_v4.json',
             'metadata/baliphy_joint_node_logger_future_models_20261003_v5.json',
             'metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json',
             'metadata/baliphy_stack_followup_plan_20261003_v1.json']
    for name in paths:
        plan = json.loads(Path(name).read_text());verify(plan['pins'])
        rows_path = Path(plan.get('jobs',plan.get('future_roles')))
        rows = json.loads(rows_path.read_text())
        for row in rows:
            used.add(row['chain']['seed'])
            for k in ['source_seed','earlier_source_seed','current_short_sampler_seed']:
                if k in row:used.add(row[k])
        bindings[name]=sha(name);bindings[str(rows_path)]=sha(rows_path)
    # Earlier software fixture seeds were deliberately repeated for paired
    # encoding comparisons. Collect their actual configurations, not guessed IDs.
    for base_path in sorted(Path('data/software_audits').glob('baliphy-*')):
        if base_path.name.startswith('baliphy-scalar-json-v6-'):continue
        base = base_path.name
        for path in sorted((Path('data/software_audits')/base).rglob('configuration.json')):
            command = json.loads(path.read_text()).get('command',[])
            if '--seed' in command:used.add(int(command[command.index('--seed')+1]))
            bindings[str(path)]=sha(path)
    assert FIXTURE_SEEDS.isdisjoint(used)
    used.update(FIXTURE_SEEDS)
    return used


def rejected(action):
    try:action()
    except (AssertionError,ValueError,KeyError,TypeError):return
    raise AssertionError('Malformed record accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    root=a.output.resolve();root.mkdir(exist_ok=False)
    pins={}
    gate=Path('metadata/baliphy_joint_node_logger_software_validation_20261003_v5.json')
    qualified=json.loads(gate.read_text());verify(qualified['source_hashes'])
    assert qualified['status']=='passed_joint_node_logger_v5_cjson_full_source_and_native_qualification'
    assert (qualified['full_programs_checked'],qualified['full_roles_checked'])==(405,1620)
    pins[str(gate)]=sha(gate)
    old_manifest=Path('metadata/baliphy_joint_node_logger_future_models_20261003_v5.json')
    old=json.loads(old_manifest.read_text());verify(old['pins'])
    roles_path=Path(old['future_roles']);roles=json.loads(roles_path.read_text())
    assert len(roles)==1620
    programs={r['chain']['program']:r['chain']['program_sha256'] for r in roles}
    assert len(programs)==405
    for name,d in programs.items():
        assert sha(name)==d
        original=Path(name).read_text();corrected=transform(original)
        assert restore(corrected)==original
        assert corrected.count('runContextAction action context')==1
        assert corrected.count(';addLogger')==original.count(';addLogger')
        pins[name]=d
    source=Path(next(iter(programs))).read_text()
    rejected(lambda:transform(transform(source)))
    rejected(lambda:transform(source.replace(';logParamsJSON <- if jsonEnabled then jsonLogger ','')))
    rejected(lambda:transform(source+'\nimport qualified Data.JSON.Encoding as E'))
    used=seed_census(pins)
    ids=[r['chain']['chain_id'] for r in roles]
    seeds=fresh_seeds(ids,used,SEED_NAMESPACE)
    assert seeds==fresh_seeds(list(reversed(ids)),used,SEED_NAMESPACE)
    assert len(seeds)==1620 and set(seeds.values()).isdisjoint(used)
    for q in [Path(__file__),HELPER,old_manifest,roles_path,
              *[Path('scripts')/n for n in ['baliphy_scalar_json_logger_v6c.py','read_baliphy_scalar_json_v6b.py',
                'baliphy_joint_node_logger_v5.py','baliphy_joint_node_logger_v3.py',
                'independent_joint_ancestral_frames.py','independent_native_ancestral_alignment.py',
                'ancestral_chain_attempt.py','reference_measurement_union_sources.py']]]:pins[str(q)]=sha(q)
    binary=Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    api=binary.parent.parent/'lib/bali-phy/haskell'
    # Bind every installed Haskell source used by the native compilation closure.
    libraries=list(sorted(api.rglob('*.hs')))
    baseline=Path('data/software_audits/baliphy-joint-node-logger-20261003-v5').resolve()
    fixture=Path('data/software_audits/baliphy-joint-node-logger-20261003-v4').resolve()
    observed={k:v.replace('-','') for k,v in fasta_records((fixture/'alignment.faa').read_text().splitlines()).items()}
    def native(name,text,seed=None):
        program=root/(name+'.hs');program.write_text(text)
        paths=[Path('/usr/bin/prlimit'),binary,program,*libraries]
        if seed is None:paths.append(root/'nonfinite-values.txt')
        command=[str(paths[0]),'--as='+str(12*2**30),'--cpu=300','--fsize='+str(64*2**20),'--',str(binary)]
        if seed is not None:
            command+=['--seed',str(seed)]
            paths += [fixture/'alignment.faa',fixture/'tree.nwk']
        command+=['run',str(program)]
        if seed is not None:command+=['--iterations','20','--log-format','json,tsv','--name','independent-chain']
        config=dict(command=command,timeout_seconds=300,pins={str(q):sha(q) for q in paths})
        target=root/'native'/name;assert not target.exists()
        receipt=run_attempt(target,config);result=json.loads(receipt.read_text())
        assert result['exit_code']==0, str(receipt)
        pins.update(config['pins']);pins[str(receipt)]=sha(receipt)
        pins[str(target/'configuration.json')]=sha(target/'configuration.json')
        for q,d in result['artifacts'].items():pins[str(receipt.parent/q)]=d
        print('scalar_v6_native',name,'exited_zero',flush=True)
        if seed is None:return receipt.parent
        directories=list(receipt.parent.glob('independent-chain-*'));assert len(directories)==1
        return directories[0]
    nonfinite_input=root/'nonfinite-values.txt'
    nonfinite_input.write_text('inf\n-inf\nnan\n')
    probe_imports='''{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ExtendedDefaultRules #-}
module Main where
import Probability
import Probability.Logger
import MCMC.Types (runContextAction)
import Compiler.RealFloat (isNaN, isInfinite, isNegativeZero, encodeFloat, isDenormalized, decodeFloat)
import qualified Data.JSON as J
import qualified Data.JSON.Encoding as E
import qualified Data.JSON.Types.Internal as JSONInternal
import Data.JSON ((.=))
import qualified Data.Text as Text
import qualified Data.Text.IO as T
import System.IO
'''
    probe=probe_imports+HELPER.read_text()+'''
probeValuesV6 :: [Double]
probeValuesV6 = [2.34e-10,2.34e-11,2.34e-20,2.34e-50,2.34e-100,2.34e-101,2.34e10,2.34e20,0,1,0.25,4,-2.34e-10,2.34e-30,-2.34e20,(encodeFloat 1 (-1074)),(1e308 * 1.5),(negate (encodeFloat 0 0))]
fields plus minus nan = [(J.toJSONKey "finite",J.Array (map J.FNumber probeValuesV6)),
          (J.toJSONKey "special/key/",J.Object [(J.toJSONKey "weird/~",J.Array
            [J.FNumber plus,J.FNumber minus,J.FNumber nan,J.Null,J.Bool True,J.String (Text.pack "literal")])])]
checks = [if isDenormalized x then decodeFloat x == (1, -1074)
          else let y = (read (Text.unpack (J.cjsonToText (J.toCJSON x))) :: Double)
               in x == y && (x /= 0 || isNegativeZero x == isNegativeZero y)
          | x <- probeValuesV6]

main = do
  contents <- T.readFile PROBE_INPUT_PLACEHOLDER
  let [plus, minus, nan] = map read (lines (Text.unpack contents)) :: [Double]
      loggedFields = fields plus minus nan
  T.putStrLn (J.fromEncoding (projectV6EncodeRecord 0 (projectV6ContextValue loggedFields) loggedFields))
  T.putStrLn (J.cjsonToText (J.toCJSON checks))
  T.putStrLn (J.cjsonToText (J.toCJSON (projectV6Sanitize (J.Object loggedFields))))
  T.putStrLn (J.cjsonToText (J.toCJSON [isNegativeZero (probeValuesV6 !! 8), isNegativeZero (last probeValuesV6)]))
'''
    probe=probe.replace('PROBE_INPUT_PLACEHOLDER',json.dumps(str(nonfinite_input)))
    directory=native('value-probe',probe)
    lines=(directory/'stdout.log').read_text().splitlines();assert len(lines)==4
    record=load(lines[0]);assert load(lines[1])==[True]*18
    assert load(lines[3])==[False,True]
    assert record['parameters//']==load(lines[2])
    assert record['parameters//']['finite'][15]==float.fromhex('0x0.0000000000001p-1022')
    first=validate_record(record);second=read_record(record)
    assert first['context']['numeric_leaves']==first['parameters']['numeric_leaves']==21
    assert len(first['context']['nonfinite'])==len(first['parameters']['nonfinite'])==3
    assert len(second['nonfinite_reviews'])==6 and not second['finite_record']
    assert set(first['parameters']['nonfinite'].values())=={'nan','positive_infinity','negative_infinity'}
    mutations=[]
    for section,auditkey in [('statistics//',QUALITY_KEY),('parameters//','numericParameterQuality//')]:
        for case in ['missing_tag','wrong_kind','duplicate_tag','negative_index','bool_index','wrong_count','null_overlap','missing_literal','tag_finite','untagged_null','nonstandard_number','wrong_placeholder']:
            bad=copy.deepcopy(record)
            audit=bad['statistics//'][auditkey] if section=='statistics//' else bad[auditkey]
            if case=='missing_tag':audit['nonfinite'].pop()
            elif case=='wrong_kind':audit['nonfinite'][0]['kind']='unknown'
            elif case=='duplicate_tag':audit['nonfinite'].append(copy.deepcopy(audit['nonfinite'][0]))
            elif case=='negative_index':audit['nonfinite'][0]['path'][-1]=-1
            elif case=='bool_index':audit['nonfinite'][0]['path'][-1]=True
            elif case=='wrong_count':audit['numericLeafCount']+=1
            elif case=='null_overlap':audit['literalNullPaths'].append(audit['nonfinite'][0]['path'])
            elif case=='missing_literal':audit['literalNullPaths']=[]
            elif case=='tag_finite':audit['nonfinite'][0]['path']=['finite',0]
            elif case=='untagged_null':bad[section]['finite'][0]=None
            elif case=='nonstandard_number':bad[section]['finite'][0]=float('inf')
            else:bad[section]['special/key/']['weird/~'][0]='__project_scalar_v6__:negative_infinity'
            rejected(lambda:validate_record(bad));rejected(lambda:read_record(bad))
            mutations.append(section+':'+case)
    pairs=[];frames=[];mapped=[]
    for index,prior in enumerate(['broad','centered','package']):
        old_path=baseline/(prior+'-numeric-only.hs');pins[str(old_path)]=sha(old_path)
        before=native(prior+'-v5',old_path.read_text(),20269001+index)
        after=native(prior+'-v6',transform(old_path.read_text()),20269001+index)
        unchanged=['C1.log','C1.log.column-map.json','runtime-tree.nwk','C1.P1.fastas','C1.P1.site-property-samples.jsonl']
        for name in unchanged:assert (before/name).read_bytes()==(after/name).read_bytes(),(prior,name)
        old_records=[load(line) for line in (before/'C1.log.json').read_text().splitlines()[1:]]
        records=[load(line) for line in (after/'C1.log.json').read_text().splitlines()[1:]]
        assert len(records)==len(old_records)==21
        for row in records:
            primary=validate_record(row);independent=read_record(row)
            assert primary['finite_record']==independent['finite_record'] is True
        mapped.append(dict(prior=prior,**compare_tsv(after)))
        labels,tips=tree_labels((after/'runtime-tree.nwk').read_text())
        joint=[load(line) for line in (after/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
        assert [r['iter'] for r in joint]==[0,10,20]
        for row in joint:
            sequences=fasta_records(row['alignmentLines']);assert set(sequences)==labels
            checked=validate_frame(row,row['iter'],sequences,tips)
            independently,arrays=decode(row,row['iter'],observed,sorted(labels),{n:n for n in sorted(labels-tips)})
            assert checked['ancestral_pairs']==independently['ancestral_pairs']
            output=root/(prior+'-'+str(row['iter'])+'.npz')
            write_arrays(output,arrays);verify_arrays(output,arrays);pins[str(output)]=sha(output)
            frames.append(dict(prior=prior,iteration=row['iter'],joint=checked,independent=independently))
        pairs.append(dict(prior=prior,seed=20269001+index,old_directory=str(before),new_directory=str(after),
                          byte_identical_files=unchanged,probability_model_changed=False))
    assert len(frames)==9 and sum(r['joint']['native_ancestors'] for r in frames)==36
    verify(pins)
    result=dict(status='passed_scalar_json_v6_full_source_and_paired_native_qualification',
        checked_utc=datetime.now(timezone.utc).isoformat(),full_programs_checked=405,full_roles_checked=1620,
        exactly_reversible_from_v5=True,context_action_evaluated_once=True,native_runs=7,
        paired_prior_checks=pairs,scalar_rows_checked=63,mapped_scalar_comparisons=mapped,
        joint_frames_checked=9,ancestral_records_checked=36,serialized_array_checks=9,
        native_finite_value_roundtrips=17,native_smallest_subnormal_construction_and_python_readback=1,native_positive_and_negative_zero_signs_verified=True,native_nonfinite_tags_checked=6,native_literal_null_tags_checked=2,
        malformed_record_cases_rejected_by_both_readers=mutations,future_seed_namespace=SEED_NAMESPACE,
        deterministic_disjoint_future_seeds=1620,forbidden_seed_count=len(used),fixture_seed_values=sorted(FIXTURE_SEEDS),
        project_scalar_schema=SCHEMA,source_hashes=pins,full_grid_native_execution_launched=False,
        probability_model_changed=False,installed_binary_or_library_changed=False,original_outputs_changed=False,
        gpu=False,new_cost_usd=0,scientific_eligibility=False,posterior_qualified=False,
        scope='Full405program/1620role static reversibility; all three priors paired at fresh fixture seeds '
        'for20iterations preserve original TSV/mapping/tree/ancestral alignment/joint frames byte-exact. '
        '63scalar rows read by two validators and independently mapped to TSV at rel2e-13 abs0; '
        '17native finite probes roundtrip exactly (including zero sign), plus one exact native smallest-subnormal construction/Python readback, six explicit nonfinite tagged string placeholders and '
        'two literal-null tags,24malformed records rejected by both readers. No full-grid startup/horizon, '
        'crash repair, recovered historical numbers, posterior adequacy or biological aim completion.')
    with a.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','paired_prior_checks','mapped_scalar_comparisons']},indent=2),flush=True)


if __name__=='__main__':main()
