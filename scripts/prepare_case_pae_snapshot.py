#!/usr/bin/env python3
"""Bind all opposing-scale case models to exact sequence/version PAE requests."""
import csv,hashlib,json
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    base=Path('results/structural_comparisons');sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name;assert sha(p)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    cases=base/'whole-domain-case-dossiers-20260927-v1'
    with checked(cases,'whole_reference_links.tsv').open() as f:links=list(csv.DictReader(f,delimiter='\t'))
    wanted={r['triad_id'] for r in links};expected={}
    with checked(base/'whole-protein-common-residues-20260927-v1','model_triads.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['triad_id'] in wanted:
                for role in ['a','b','reference']:
                    key=r[role+'_model'],int(r[role+'_version']);digest=r[role+'_sequence_sha256']
                    assert key not in expected or expected[key]==digest;expected[key]=digest
    found={}
    for name in ['duplication-alignment-inputs-20260926-v1','duplication-reference-alignment-inputs-20260926-v1']:
        for line in checked(base/name,'inputs.jsonl').open():
            r=json.loads(line);key=r['model_id'],r['version']
            if key in expected and r['mask']=='full':
                assert hashlib.sha256(r['sequence'].encode()).hexdigest()==expected[key]
                assert r['status']=='ready' and len(r['sequence'])==r['original_length'] and r['original_positions']==list(range(1,r['original_length']+1))
                assert sha(r['path'])==r['sha256'];sources[r['path']]=r['sha256']
                m=dict(model_id=key[0],version=key[1],length=r['original_length'],sequence_sha256=expected[key],pae_url=f'https://alphafold.ebi.ac.uk/files/{key[0]}-predicted_aligned_error_v{key[1]}.json')
                assert key not in found or found[key]==m;found[key]=m
    assert set(found)==set(expected) and len(found)==39
    output=base/'whole-domain-case-pae-inputs-20260927-v1';output.mkdir(exist_ok=False)
    models=[found[key] for key in sorted(found)];p=output/'model_provenance.json';p.write_text(json.dumps(models,indent=2)+'\n');assert json.loads(p.read_text())==models
    cache=Path('data/structures/pae');cached=sum((cache/f"{r['model_id']}-v{r['version']}.receipt.json").exists() for r in models)
    receipt=dict(status='complete_case_pae_request_snapshot',models=len(models),existing_cache_receipts=cached,source_hashes=sources,script_sha256=sha(__file__),
        resources=dict(cpus=1,memory_gib=2,swap_gib=0,http_workers=2,output_gib=0.5,planning_minutes=[1,30],maximum_matrix_length=max(r['length'] for r in models),total_matrix_entries=sum(r['length']**2 for r in models),paid_resources=False),
        artifacts={p.name:sha(p)},scope='All 39 exact accession/version/sequence models from 13 cases; local full-coordinate sequences and hashes verified. Cache presence is not validation. No new folding, alternate-version substitution or account required; retrieval failures must remain explicit.')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
