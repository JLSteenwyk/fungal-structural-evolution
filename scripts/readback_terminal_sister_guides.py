#!/usr/bin/env python3
"""Check the complete guide union against both independently audited inventories."""
import argparse,csv,json,sqlite3,time
from collections import Counter
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha


def expected_row(a,b):
    rows={'profile':a,'mafft':b};present={g for g,r in rows.items() if r is not None}
    first=next(r for r in rows.values() if r is not None)
    result={k:first[k] for k in ['gene_a','gene_b','taxon_a','taxon_b']}
    result['presence']='both' if len(present)==2 else next(iter(present))+'_only'
    copied=('family','gene_node','candidate_status','model_coverage','parent_reported_duplication','sequence_tip_a','sequence_tip_b','sequence_pair_distance','model_a','version_a','model_b','version_b')
    candidates=set()
    for guide,row in rows.items():
        for field in copied:result[guide+'_'+field]=row[field] if row is not None else ''
        if row is not None:
            assert all(row[k]==first[k] for k in ['gene_a','gene_b','taxon_a','taxon_b'])
            if row['candidate_status']=='cross_taxon_unreported_candidate' and row['model_coverage'] in {'identical_model','two_distinct_models'}:candidates.add(guide)
    for name,field in [('candidate_class_agrees','candidate_status'),('family_label_agrees','family')]:
        result[name]=str(int(a[field]==b[field])) if len(present)==2 else ''
    result['modeled_candidate_in_either_guide']=str(int(bool(candidates)))
    result['modeled_candidate_in_both_guides']=str(int(len(candidates)==2))
    result['both_parents_unreported']=str(int(len(present)==2 and all(r['parent_reported_duplication']=='0' for r in rows.values())))
    return result


def merge_sources(left,right):
    left,right=iter(left),iter(right);a=next(left,None);b=next(right,None)
    while a is not None or b is not None:
        ka=a[:2] if a else None;kb=b[:2] if b else None
        if a is not None and b is not None and ka==kb:
            yield json.loads(a[2]),json.loads(b[2]);a=next(left,None);b=next(right,None)
        elif a is not None and (b is None or ka<kb):
            yield json.loads(a[2]),None;a=next(left,None)
        else:
            yield None,json.loads(b[2]);b=next(right,None)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.config.read_text());ch=sha(args.config)
    def check_config():
        assert sha(args.config)==ch
        for p,h in config['pins'].items():assert sha(p)==h,p
    check_config();dep=config['dependency']
    while True:
        try:
            proc=psutil.Process(dep['pid'])
            if proc.create_time()!=dep['created'] or proc.status()==psutil.STATUS_ZOMBIE:break
            assert proc.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    check_config();plan=json.loads(Path(config['plan']).read_text());root=Path(plan['output']);source=Path(plan['source'])
    receipt=json.loads((root/'receipt.json').read_text());sr=json.loads((source/'receipt.json').read_text());audit=json.loads(Path(plan['readback']).read_text())
    assert receipt['status']=='complete_terminal_sister_guide_comparison_pending_readback'
    assert receipt['plan_sha256']==sha(config['plan'])
    assert audit['status']=='passed_full_terminal_sister_inventory_readback'
    assert receipt['source_receipt_sha256']==audit['producer_receipt_sha256']==sha(source/'receipt.json')
    assert receipt['source_readback_sha256']==sha(plan['readback'])
    bindings={str(root/'receipt.json'):sha(root/'receipt.json'),str(source/'receipt.json'):sha(source/'receipt.json'),plan['readback']:sha(plan['readback'])}
    for directory,r in [(root,receipt),(source,sr)]:
        bindings.update({str(directory/name):h for name,h in r['artifacts'].items()})
    def verify():
        check_config()
        for p,h in bindings.items():assert sha(p)==h,p
    verify()
    con=sqlite3.connect('file:'+str((root/'pair_sources.sqlite').resolve())+'?mode=ro',uri=True);con.execute('PRAGMA cache_size=-65536')
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    checked={}
    for guide in ['profile','mafft']:
        n=0
        with (source/(guide+'_terminal_sisters.tsv')).open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                found=con.execute(f'SELECT record FROM {guide} WHERE a=? AND b=?',(row['gene_a'],row['gene_b'])).fetchall()
                assert len(found)==1 and json.loads(found[0][0])==row
                n+=1
        assert con.execute(f'SELECT count(*) FROM {guide}').fetchone()[0]==n==next(g['pairs'] for g in sr['guides'] if g['guide']==guide)
        checked[guide]=n;print(guide,n,'source rows verified',flush=True)
    total=selected=0;presence=Counter();classes=Counter();counts=Counter()
    left=con.execute('SELECT a,b,record FROM profile ORDER BY a,b');right=con.execute('SELECT a,b,record FROM mafft ORDER BY a,b')
    with (root/'guide_comparison.tsv').open() as full,(root/'modeled_candidate_union.tsv').open() as subset:
        fullrows=iter(csv.DictReader(full,delimiter='\t'));selectedrows=iter(csv.DictReader(subset,delimiter='\t'))
        for a,b in merge_sources(left,right):
            expected=expected_row(a,b);row=next(fullrows,None);assert row==expected
            total+=1;presence[row['presence']]+=1;classes[row['profile_candidate_status']+'|'+row['mafft_candidate_status']]+=1
            if row['modeled_candidate_in_either_guide']=='1':
                assert next(selectedrows,None)==row;selected+=1;counts[row['presence']]+=1
                counts['both_candidate']+=int(row['modeled_candidate_in_both_guides'])
                counts['both_candidate_parents_unreported']+=int(row['modeled_candidate_in_both_guides'])*int(row['both_parents_unreported'])
            if total%100000==0:print(total,'union pairs verified',flush=True)
        assert next(fullrows,None) is None and next(selectedrows,None) is None
    con.close()
    assert total==receipt['union_pairs'] and selected==receipt['modeled_candidate_union']
    assert dict(presence)==receipt['presence'] and dict(classes)==receipt['classification_cross_tab'] and dict(counts)==receipt['modeled_candidate_counts']
    verify()
    result=dict(status='passed_full_terminal_sister_guide_comparison_readback',config_sha256=ch,producer_receipt_sha256=sha(root/'receipt.json'),source_rows_checked=checked,union_pairs=total,modeled_candidate_union=selected,presence=dict(presence),modeled_candidate_counts=dict(counts),checker_sha256=sha(__file__),scope='Every SQLite source record checked against original audited TSVs; independent sorted-cursor union reconstructs every output field and exact modeled subset, including absent guides and disagreements. No biological orthology or matched-control eligibility established.')
    target=Path(config['output']);assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
