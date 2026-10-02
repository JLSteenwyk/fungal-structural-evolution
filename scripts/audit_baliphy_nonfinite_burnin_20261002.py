#!/usr/bin/env python3
"""Locate every flagged nonfinite cell in the full original scalar census."""
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def main():
    source=Path('metadata/baliphy_full_native_scalar_range_census_20261002.json')
    r=json.loads(source.read_text()); assert r['status']=='complete_full_original_baliphy_native_scalar_range_census'
    table=Path('docs/tables/baliphy_full_native_scalar_ranges_20261002.tsv'); assert sha(table)==r['artifacts'][str(table)]
    with table.open() as f:
        paths={row['chain_id']:row for row in csv.DictReader(f,delimiter='\t') if row['variable']=='ASRV.Gamma:alpha'}
    flagged=[c for c in r['chain_summaries'] if c['nonfinite_by_variable']]
    assert len(flagged)==r['chains_with_nonfinite_scalar_values']==134
    events=[];bindings={str(source):sha(source),str(table):sha(table),str(Path(__file__)):sha(__file__)}
    for chain in flagged:
        row=paths[chain['chain_id']];path=Path(row['source_log']);digest=row['source_log_sha256']
        assert sha(path)==r['source_hashes'][str(path)]==digest;bindings[str(path)]=digest
        found=0
        with path.open() as f:
            for values in csv.DictReader(f,delimiter='\t'):
                for variable,raw in values.items():
                    if math.isfinite(float(raw)):continue
                    found+=1;events.append(dict(chain_id=chain['chain_id'],status=chain['status'],iteration=int(values['iter']),variable=variable,raw_value=raw))
        assert found==sum(chain['nonfinite_by_variable'].values())
    assert len(events)==sum(r['nonfinite_by_variable'].values())==335
    counts={str(c):sum(e['iteration']>c for e in events) for c in [250,500]}
    assert counts=={'250':0,'500':0}
    for p,d in bindings.items():assert sha(p)==d
    result=dict(status='verified_full_original_nonfinite_scalar_burnin_location',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_original_chains=1620,flagged_chains=134,nonfinite_cells=335,earliest_iteration=min(e['iteration'] for e in events),latest_iteration=max(e['iteration'] for e in events),
        nonfinite_cells_after_original_cutoffs=counts,events=events,source_hashes=bindings,scientific_eligibility=False,
        scope='The full1620chain original census proves allnonfinite cells occur in these134logs. Every335raw nonfinite cell located before bothoriginalburn-in cutoffs; none enters the reportedpostburn scalar diagnostics. This doesnot certify initialization, numericalmodel adequacy, posterior convergence or length/category mixing; original diagnostics already qualifyzeroquartets. No thresholds changed or logsmodified. Live memory-recovery attempts not combined with originals.')
    with Path('metadata/baliphy_full_nonfinite_burnin_audit_20261002.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['events','source_hashes','scope']},indent=2))


if __name__=='__main__':main()
