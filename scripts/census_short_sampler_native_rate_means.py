#!/usr/bin/env python3
"""Census every saved native rate frame without altering or accepting damaged values."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
import math
from pathlib import Path

from audit_selected_taxon_identity_snapshot_v2 import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    plan=json.loads(a.plan.read_text());completion=Path(plan['completion']);closed=json.loads(completion.read_text())
    assert closed['status']=='complete_verified_full_reference_short_sampler_qualification'
    assert closed['checked_sampler_attempts']==1620 and closed['unsuccessful_sampler_attempts']==0
    assert sha(closed['full_hash_archive'])==closed['full_hash_archive_sha256']
    path=Path(plan['output'])/'dispositions.json';rows=json.loads(path.read_text());assert len(rows)==1620
    bindings={str(a.plan):sha(a.plan),str(completion):sha(completion),str(path):sha(path),str(Path(__file__)):sha(__file__)}
    bad=[];frames=0;maximum=0
    for row in rows:
        receipt=Path(row['native_receipt']);assert sha(receipt)==row['native_receipt_sha256']
        r=json.loads(receipt.read_text());assert r['exit_code']==0
        directories=list(receipt.parent.glob('independent-chain-*'));assert len(directories)==1
        data=directories[0]/'C1.P1.site-property-samples.jsonl'
        digest=r['artifacts'][str(data.relative_to(receipt.parent))];assert sha(data)==digest
        bindings[str(receipt)]=sha(receipt);bindings[str(data)]=digest
        records=[json.loads(line) for line in data.read_text().splitlines()]
        assert [f['iter'] for f in records]==[0,10,20]
        for f in records:
            rates=f['properties']['rate'];assert len(rates)==4 and all(len(values)==20 for values in rates)
            assert all(all(type(v) in [int,float] and math.isfinite(v) and v>=0 for v in values) and all(v==values[0] for v in values) for values in rates)
            values=[values[0] for values in rates];mean=math.fsum(values)/4;error=abs(mean-1)
            frames+=1;maximum=max(maximum,error)
            if not math.isclose(mean,1,rel_tol=1e-10,abs_tol=1e-10):
                bad.append(dict(chain_id=row['chain_id'],family=row['family'],iteration=f['iter'],category_rates=values,
                    mean=mean,absolute_error=error,source=str(data),scalar_log=str(directories[0]/'C1.log')))
    assert frames==4860
    for path,digest in bindings.items():assert sha(path)==digest,path
    result=dict(status='complete_read_only_saved_native_rate_mean_census',checked_utc=datetime.now(timezone.utc).isoformat(),
        native_roles=1620,full_native_frames=frames,strict_mean_one_failures=len(bad),maximum_absolute_deviation=maximum,
        affected_family_counts=dict(Counter(r['family'] for r in bad)),failed_records=bad,source_hashes=bindings,
        tolerance=dict(rel_tol=1e-10,abs_tol=1e-10),historical_float_corruption_assessment_complete=False,scientific_eligibility=False,
        scope='Every original saved frame/receipt hash checked and normalization check retained exactly. Mean-check failures are not a complete exponent-corruption census; rates not reconstructed, renormalized or adopted. No native restart, source/output mutation or biological acceptance.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['failed_records','source_hashes','scope']},indent=2))


if __name__=='__main__':main()
