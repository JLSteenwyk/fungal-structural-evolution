#!/usr/bin/env python3
"""Describe taxon/family concentration across the full verified matching grid."""
import argparse,csv,gzip,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    p=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in p['pins'].items():assert sha(path)==h,path
    verify()
    for label in ['graph','selection','balance']:
        r=json.loads(Path(p[label]+'/receipt.json').read_text());proof=json.loads(Path(p[label+'_proof']).read_text())
        assert proof['producer_receipt_sha256']==sha(p[label]+'/receipt.json')
        assert proof['status']==p['proof_status'][label]
        for name in p['required_artifacts'][label]:assert sha(p[label]+'/'+name)==r['artifacts'][name]
    nodes={};baseline=defaultdict(Counter)
    for line in open(p['graph']+'/target_nodes.jsonl'):
        n=json.loads(line);assert n['node_id'] not in nodes
        nodes[n['node_id']]=(n['guide'],n['taxon_id'],n['family'])
        for dimension,key in [('taxon',n['taxon_id']),('family',n['family'])]:baseline[n['guide'],dimension][key]+=1
    coverage={}
    for row in csv.DictReader(open(p['balance']+'/selection_coverage.tsv'),delimiter='\t'):
        k=row['guide'],row['policy'],row['scenario_id'];assert k not in coverage;coverage[k]=row
    counts=defaultdict(Counter);total=0
    with gzip.open(p['selection']+'/selections.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            guide,taxon,family=nodes[row['target_id']];k=guide,row['policy'],row['scenario_id'];assert k in coverage
            counts[(*k,'taxon')][taxon]+=1;counts[(*k,'family')][family]+=1;total+=1
    assert total==json.loads(Path(p['selection']+'/receipt.json').read_text())['selected_records']
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False)
    nf=['guide','policy','scenario_id','dimension','group_id','original_targets','selected_targets','share_of_selected','fraction_of_original_group_selected']
    sf=['guide','policy','scenario_id','dimension','original_targets','selected_targets','original_groups','selected_groups','groups_without_matches','maximum_group_count','top_five_share','inverse_concentration']
    rows_written=0
    with (out/'group_counts.tsv').open('w') as f,(out/'concentration.tsv').open('w') as g:
        w=csv.DictWriter(f,nf,delimiter='\t',lineterminator='\n');w.writeheader();v=csv.DictWriter(g,sf,delimiter='\t',lineterminator='\n');v.writeheader()
        for k,cov in sorted(coverage.items()):
            for dimension in ['taxon','family']:
                c=counts[(*k,dimension)];b=baseline[k[0],dimension];n=sum(c.values());original=sum(b.values())
                assert n==int(cov['matched_target_records']) and original==int(cov['all_target_records'])
                assert len(c)==int(cov['matched_taxa' if dimension=='taxon' else 'matched_families'])
                common=dict(zip(['guide','policy','scenario_id'],k),dimension=dimension)
                for group,count in sorted(c.items()):
                    assert count<=b[group]
                    w.writerow(dict(common,group_id=group,original_targets=b[group],selected_targets=count,share_of_selected=count/n,fraction_of_original_group_selected=count/b[group]));rows_written+=1
                v.writerow(dict(common,original_targets=original,selected_targets=n,original_groups=len(b),selected_groups=len(c),groups_without_matches=len(b)-len(c),maximum_group_count=max(c.values(),default=0),top_five_share=sum(sorted(c.values(),reverse=True)[:5])/n if n else '',inverse_concentration=n*n/sum(x*x for x in c.values()) if n else ''))
    verify()
    result=dict(status='complete_selected_control_concentration_pending_readback',plan_sha256=ph,selected_records=total,strata=len(coverage),summary_rows=2*len(coverage),group_rows=rows_written,artifacts={x.name:sha(x) for x in out.iterdir()},scope='Selected target records aggregated separately by taxon and family in every guide/policy/scenario; group table contains represented groups only, summaries explicitly count original groups with no matches. Original means all modeled terminal targets in the guide. Inverse concentration is N squared / sum(group count squared), a descriptive equivalent number of equally represented groups, NOT an independent sample size or phylogenetic correction. No structural outcomes or effects used.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
