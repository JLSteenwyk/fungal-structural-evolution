#!/usr/bin/env python3
"""Bind the fully audited case search to complete experimental metadata retrieval."""
import json
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    source=Path('results/experimental_structures/whole-domain-case-sequence-search-20260927-v1')
    audit=Path('results/experimental_structures/whole-domain-case-sequence-readback-20260927-v1/receipt.json');a=json.loads(audit.read_text())
    assert a['status']=='passed_full_case_experimental_sequence_search_readback'
    rp=source/'receipt.json';assert a['source_hashes'][str(rp)]==sha(rp)
    r=json.loads(rp.read_text());p=source/'polymer_entity_ids.txt';assert sha(p)==r['artifacts'][p.name]
    ids=p.read_text().split();assert len(ids)==len(set(ids))==a['unique_entities']==1205
    entries={x.rsplit('_',1)[0] for x in ids}
    out=Path('results/experimental_structures/whole-domain-case-metadata-inputs-20260927-v1');out.mkdir(exist_ok=False)
    target=out/p.name;target.write_bytes(p.read_bytes());assert sha(target)==sha(p)
    result=dict(status='complete_audited_case_experimental_entity_inventory',unique_polymer_entities=len(ids),entries=len(entries),source_hashes={str(rp):sha(rp),str(p):sha(p),str(audit):sha(audit)},script_sha256=sha(__file__),artifacts={target.name:sha(target)},resources=dict(cpus=1,memory_gib=2,swap_gib=0,http_workers=2,requests=len(ids)+len(entries),output_gib=.5,planning_minutes=[10,90],paid_resources=False),scope='Every experimental entity returned for all 39 case queries retained; no quality or favorable-geometry selection before metadata retrieval.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
