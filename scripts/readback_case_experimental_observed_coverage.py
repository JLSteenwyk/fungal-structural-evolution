#!/usr/bin/env python3
"""Recompute all resolved-query coverage counts from frozen residue correspondences."""
import csv,gzip,json,subprocess,time
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import psutil
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results/experimental_structures')
ROOT=BASE/'whole-domain-case-observed-coverage-20260927-v1'
OUTPUT=BASE/'whole-domain-case-observed-coverage-readback-20260927-v1'
CASE=['family','gene_a','gene_b','pfam_accession']
FIELDS=CASE+['role','interval_id','sequence_id','entity_id','context_index','screen']


def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    sources={};lp=Path('metadata/case_experimental_observed_coverage_launch_20260927.json');launch=json.loads(lp.read_text());sources[str(lp)]=sha(lp)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact observed coverage producer',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0') and sha(launch['cmdline'][1])==launch['script_sha256']
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name;assert sha(p)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    pairs={}
    with gzip.open(checked(BASE/'whole-domain-case-residue-pairs-20260927-v1','residue_correspondences.jsonl.gz'),'rt') as f:
        for r in map(json.loads,f):pairs[r['sequence_id'],r['entity_id'],str(r['context_index'])]=r['query_subject_position_pairs']
    ca_path=BASE/'whole-domain-case-ca-readback-20260927-v1/receipt.json';ca=json.loads(ca_path.read_text());assert ca['status']=='complete_full_case_CA_raw_atom_and_grid_readback';sources[str(ca_path)]=sha(ca_path)
    mapping=BASE/'whole-domain-case-ca-mapping-20260927-v1';config=mapping/'config.json';assert sha(config)==ca['source_hashes'][str(config)];excluded=json.loads(config.read_text())['deferred_entries']
    observed=defaultdict(dict)
    for path,digest in ca['source_hashes'].items():
        if not path.endswith('.residues.tsv.gz'):continue
        p=Path(path);assert sha(p)==digest;sources[path]=digest;entry=p.name.split('.')[0]
        with gzip.open(p,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                key=(r['model_number'],r['label_asym_id']);entity=entry+'_'+r['entity_id'];observed[entity].setdefault(key,set())
                if r['CA_status']=='unambiguous_full_occupancy_CA':observed[entity][key].add(int(r['label_seq_id']))
    old=BASE/'whole-domain-case-domain-coverage-20260927-v1';original=rows(checked(old,'alignment_interval_screens.tsv'));prior_summary=rows(checked(old,'case_summary.tsv'))
    needed={r['interval_id'] for r in original};bounds={}
    for line in checked(Path('results/structural_comparisons/duplication-domain-inputs-20260926-v1'),'inputs.jsonl').open():
        r=json.loads(line)
        if r['mask']=='full' and r['interval_id'] in needed:bounds[r['interval_id']]=r
    required={}
    for r in rows(checked(old,'shared_entity_screens.tsv')):
        ck=tuple(r[k] for k in CASE)
        if ck in required:assert required[ck]==int(r['required_role_intervals'])
        required[ck]=int(r['required_role_intervals'])
    expected=set()
    for r in original:
        if r['entity_id'].rsplit('_',1)[0] in excluded:continue
        assert r['entity_id'] in observed
        for model,chain in observed[r['entity_id']]:expected.add(tuple(r[k] for k in FIELDS)+(model,chain))
    groups=defaultdict(lambda:{'domain':set(),'both':set()});count=0;cache={}
    with checked(ROOT,'observed_alignment_interval_screens.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=tuple(r[k] for k in FIELDS)+(r['model_number'],r['label_asym_id']);assert key in expected;expected.remove(key);count+=1
            ck=tuple(r[k] for k in CASE);bounds_row=bounds[r['interval_id']]
            memo=tuple(r[k] for k in ['sequence_id','entity_id','context_index','model_number','label_asym_id','interval_id'])
            if memo not in cache:
                present=observed[r['entity_id']][r['model_number'],r['label_asym_id']]
                matched=[q for q,s in pairs[r['sequence_id'],r['entity_id'],r['context_index']] if s in present]
                inside=sum(bounds_row['start']<=q<=bounds_row['end'] for q in matched);cache[memo]=(inside,len(matched)-inside)
            inside,outside=cache[memo];length=bounds_row['end']-bounds_row['start']+1;remaining=bounds_row['original_length']-length
            assert [int(r[k]) for k in ['observed_domain_residues','observed_outside_residues','domain_length','outside_length']]==[inside,outside,length,remaining]
            n,c=(int(x[1:]) for x in r['screen'].split('_'));anchor=inside>=n and Fraction(inside,length)>=Fraction(c,100);both=anchor and outside>=n and remaining>0 and Fraction(outside,remaining)>=Fraction(c,100)
            assert int(r['domain_pass'])==anchor and int(r['domain_and_outside_pass'])==both
            g=groups[ck,r['screen'],r['entity_id'],r['model_number'],r['label_asym_id']]
            if anchor:g['domain'].add((r['role'],r['interval_id']))
            if both:g['both'].add((r['role'],r['interval_id']))
    assert not expected
    qualified=defaultdict(lambda:{'domain':set(),'both':set()});seen=set()
    for r in rows(checked(ROOT,'shared_chain_model_screens.tsv')):
        ck=tuple(r[k] for k in CASE);key=(ck,r['screen'],r['entity_id'],r['model_number'],r['label_asym_id']);assert key in groups and key not in seen;seen.add(key)
        g=groups[key];assert int(r['required_role_intervals'])==required[ck]
        assert int(r['passing_domain_role_intervals'])==len(g['domain']) and int(r['passing_both_role_intervals'])==len(g['both'])
        for name,field in [('domain','shared_domain_pass'),('both','shared_domain_and_outside_pass')]:
            assert len(g[name])<=required[ck];passed=len(g[name])==required[ck];assert int(r[field])==passed
            if passed:qualified[ck,r['screen']][name].add((r['entity_id'],r['model_number'],r['label_asym_id']))
    assert seen==set(groups)
    expected_summaries={tuple(r[k] for k in CASE)+(r['screen'],) for r in prior_summary}
    for r in rows(checked(ROOT,'case_summary.tsv')):
        ck=tuple(r[k] for k in CASE);key=(*ck,r['screen']);assert key in expected_summaries;expected_summaries.remove(key)
        for name in ['domain','both']:
            values=qualified[ck,r['screen']][name]
            assert int(r['shared_'+name+'_chain_models'])==len(values)
            assert int(r['shared_'+name+'_entities'])==len({v[0] for v in values})
            assert int(r['shared_'+name+'_entries'])==len({v[0].rsplit('_',1)[0] for v in values})
    assert not expected_summaries
    for p,digest in sources.items():assert sha(p)==digest
    OUTPUT.mkdir(exist_ok=False);result=dict(status='complete_full_observed_case_coverage_readback',terminal_state=state,source_hashes=sources,script_sha256=sha(__file__),interval_screen_rows=count,chain_model_screen_rows=len(groups),case_screen_rows=len(prior_summary),scope='Every expected alignment/interval/screen expanded across deposited chain/models; all counts independently reconstructed from frozen residue pairs and audited CA positions. Exact rational thresholds, complete shared-role qualification, and all zero-inclusive case summaries verified. No geometry or experimental independence claim.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
