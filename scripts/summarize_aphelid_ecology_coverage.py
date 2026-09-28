"""Verify preserved45-taxon results and summarize expanded source coverage."""
import csv
import json
from pathlib import Path
from ancestral_chain_attempt import sha,write_json


def load(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter='\t'))

sources={};taxa={};pairs={};summaries={}
for method in ['afdb','esmfold']:
    root=Path(f'results/ecology/qualified-{method}-aphelids47-overlap-20260928-v1')
    previous=Path(f'results/ecology/qualified-{method}-reviewed45-overlap-20260927-v1')
    receipt=json.loads((root/'receipt.json').read_text());old=json.loads((previous/'receipt.json').read_text())
    audit=Path(f'metadata/qualified_{method}_aphelids47_ecology_readback_20260928.json')
    ar=json.loads(audit.read_text());assert ar['status']=='passed_full_qualified_ecology_overlap_matrix_readback'
    assert ar['producer_receipt_sha256']==sha(root/'receipt.json')
    for base,r in [(root,receipt),(previous,old)]:
        sources[str(base/'receipt.json')]=sha(base/'receipt.json')
        for name,h in r['artifacts'].items():assert sha(base/name)==h
    sources[str(audit)]=sha(audit)
    for name,keys in [('taxa.tsv',['taxon_id']),('pairs.tsv',['taxon_a','taxon_b']),('pair_marker_columns.tsv',['taxon_a','taxon_b','marker']),('groups.tsv',['provisional_group'])]:
        old_rows=load(previous/name);new_rows=load(root/name)
        indexed={tuple(r[k] for k in keys):r for r in new_rows};assert len(indexed)==len(new_rows)
        for row in old_rows:assert indexed[tuple(row[k] for k in keys)]==row
    taxa[method]={r['taxon_id'] for r in load(root/'taxa.tsv') if int(r['eligible_markers'])>0}
    pairs[method]={tuple(sorted([r['taxon_a'],r['taxon_b']])) for r in load(root/'pairs.tsv') if int(r['markers_with_at_least_50_shared_columns'])>0}
    summaries[method]=dict(taxa_with_markers=len(taxa[method]),qualified_pairs=len(pairs[method]),new_taxa=[r for r in load(root/'taxa.tsv') if r['taxon_id'] in ['F1498855','F1243177']])
result=dict(status='complete47_taxon_ecology_coverage_and_preservation_readback',taxa=47,pairs=1081,sources=sources,
    methods=summaries,taxa_with_either_source=len(taxa['afdb']|taxa['esmfold']),pairs_with_either_source=len(pairs['afdb']|pairs['esmfold']),
    pairs_with_both_sources=len(pairs['afdb']&pairs['esmfold']),pairs_with_neither_source=1081-len(pairs['afdb']|pairs['esmfold']),
    script_sha256=sha(__file__),scope='All prior45-taxon outputs preserved. Source-specific qualified coverage only, not independent contrasts, power or ecological effects. Source union is bookkeeping; no cross-source observations pooled.')
write_json(Path('metadata/aphelids47_ecology_coverage_completed_20260928.json'),result)
print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))
