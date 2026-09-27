#!/usr/bin/env python3
"""Compare verified source-specific coverage without combining structural observations."""
import csv,json,hashlib
from collections import Counter
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
sources={};bindings={};preserved={}
for source,oldpath in [('afdb','results/ecology/qualified-afdb-recovered-overlap-20260927-v1'),('esmfold','results/ecology/qualified-esmfold-overlap-20260922-v1')]:
    root=Path(f'results/ecology/qualified-{source}-reviewed45-overlap-20260927-v1');rp=root/'receipt.json';proof=Path(f'metadata/qualified_{source}_reviewed45_ecology_readback_20260927.json');r=json.loads(rp.read_text());pr=json.loads(proof.read_text())
    assert pr['status']=='passed_full_qualified_ecology_overlap_matrix_readback' and pr['producer_receipt_sha256']==sha(rp)
    for name,h in r['artifacts'].items():assert sha(root/name)==h;bindings[str(root/name)]=h
    for path in [rp,proof]:bindings[str(path)]=sha(path)
    sources[source]={name:rows(root/(name+'.tsv')) for name in ['taxa','pairs','groups']}
    old=Path(oldpath);orr=json.loads((old/'receipt.json').read_text());bindings[str(old/'receipt.json')]=sha(old/'receipt.json')
    for name,key in [('taxa',lambda r:r['taxon_id']),('pairs',lambda r:(r['taxon_a'],r['taxon_b'])),('groups',lambda r:r['provisional_group'])]:
        path=old/(name+'.tsv');assert sha(path)==orr['artifacts'][path.name];bindings[str(path)]=sha(path)
        current={key(r):r for r in sources[source][name]};previous=rows(path)
        assert all(current[key(r)]==r for r in previous)
        preserved[source+'_'+name]=len(previous)
by_source={s:{r['taxon_id']:r for r in data['taxa']} for s,data in sources.items()}
assert by_source['afdb'].keys()==by_source['esmfold'].keys()
output=[];patterns=Counter()
for taxon,a in by_source['afdb'].items():
    e=by_source['esmfold'][taxon];assert (a['species_name'],a['state'])==(e['species_name'],e['state'])
    af=int(a['eligible_markers']);es=int(e['eligible_markers']);category='both_sources' if af and es else 'afdb_only' if af else 'esmfold_only' if es else 'neither_source';patterns[category]+=1
    output.append(dict(taxon_id=taxon,species_name=a['species_name'],state=a['state'],afdb_eligible_markers=af,esmfold_eligible_markers=es,coverage_category=category))
pairpatterns=Counter()
a={(r['taxon_a'],r['taxon_b']):r for r in sources['afdb']['pairs']};e={(r['taxon_a'],r['taxon_b']):r for r in sources['esmfold']['pairs']};assert a.keys()==e.keys()
for k in a:
    af=int(a[k]['markers_with_at_least_50_shared_columns'])>0;es=int(e[k]['markers_with_at_least_50_shared_columns'])>0
    pairpatterns['both_sources' if af and es else 'afdb_only' if af else 'esmfold_only' if es else 'neither_source']+=1
path=Path('metadata/reviewed45_ecology_source_coverage_20260927.tsv');assert not path.exists()
with path.open('w') as f:
    w=csv.DictWriter(f,list(output[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(output)
for p,h in bindings.items():assert sha(p)==h
result=dict(status='complete_verified_ecology_source_coverage_comparison',source_bindings=bindings,script_sha256=sha(__file__),preserved_original_rows=preserved,taxon_coverage=dict(patterns),pair_coverage_at_least_one_50_column_marker=dict(pairpatterns),table_sha256=sha(path),scope='All original 32-taxon, 496-pair and nine-group outputs preserved exactly within each source. Coverage categories only: both_sources does not mean the same marker or homologous retained columns across sources. No merged prediction-source structural analysis, independent transitions or ecological effects.')
Path('metadata/reviewed45_ecology_source_coverage_receipt_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_bindings'},indent=2))
