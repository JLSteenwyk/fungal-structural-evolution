#!/usr/bin/env python3
"""Compare all four audited largest-family reconstruction alternatives."""
import csv,itertools,json,subprocess,time
from pathlib import Path
import psutil
from Bio import SeqIO
from rapidfuzz.distance import Levenshtein
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/historian_largest_family_comparison_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();pins={str(pp):sha(pp)};records=[];reference_tips=None;reference_options=None
    mappings=list(csv.DictReader(Path(plan['nodes']).open(),delimiter='\t'))
    mappings=[r for r in mappings if r['guide']=='profile' and r['family']=='OG0000972' and r['dataset']=='whole'];assert len(mappings)==4
    for config in plan['collections']:
        ident=json.loads(Path(config['launch']).read_text())
        while True:
            state=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',ident['unit'],'-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
            if state['ActiveState'] in ['inactive','failed']:
                assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
                break
            assert state['ActiveState']=='active' and int(state['MainPID'])==ident['pid']
            try:
                p=psutil.Process(ident['pid']);assert p.create_time()==ident['created'] and p.cmdline()==ident['cmdline']
            except psutil.NoSuchProcess:
                time.sleep(1);continue
            time.sleep(30)
        ar=Path(config['audit']);audit=json.loads(ar.read_text());assert audit['expected_jobs']==2 and audit['dispositions']=={'independent_tip_tree_and_candidate_readback_passed':2}
        pins[str(ar)]=sha(ar)
        for f,h in audit['pins'].items():assert sha(f)==h
        for f,h in audit['artifacts'].items():assert sha(ar.parent/f)==h
        producer=json.loads(Path(config['producer_plan']).read_text())
        if reference_options is None:reference_options=producer['options']
        assert producer['options']==reference_options
        inputs=json.loads(Path(producer['input_receipt']).read_text())
        for job in inputs['jobs']:
            assert job['family']=='OG0000972' and job['proteins']==622
            folder=Path(producer['output'])/job['job_id'];rp=folder/'receipt.json';r=json.loads(rp.read_text());assert r['job']==job and r['plan_sha256']==sha(config['producer_plan'])
            assert audit['pins'][str(rp)]==sha(rp);pins[str(rp)]=sha(rp)
            for f,h in r['artifacts'].items():assert sha(folder/f)==h
            tips={r.id:str(r.seq).replace('-','') for r in SeqIO.parse(job['alignment'],'fasta')}
            if reference_tips is None:reference_tips=tips
            assert tips==reference_tips and len(tips)==622
            data=json.loads((folder/'reconstruction.json').read_text());children={}
            for parent,child,length in data['branches']:children.setdefault(parent,[]).append(child)
            ancestors={}
            for m in mappings:
                names=json.loads(m['output_nodes_json']);assert len(names)==1;node=names[0];todo=[node];desc=set()
                while todo:
                    n=todo.pop()
                    if n in children:todo.extend(children[n])
                    else:desc.add(n)
                assert desc==set(json.loads(m['retained_set_json']))
                ancestors[int(m['level'])]=data['rowData'][node].replace('-','').upper()
            floor=float(job['job_id'].split('-floor')[1].split('-resolution')[0])
            records.append(dict(job_id=job['job_id'],alignment=config['alignment'],floor=floor,tree_sha256=job['tree_sha256'],ancestors=ancestors))
    assert len(records)==4 and {(r['alignment'],r['floor']) for r in records}==set(itertools.product(['mafft','famsa'],[1e-9,1e-7]))
    rows=[]
    for a,b in itertools.combinations(records,2):
        same_floor=a['floor']==b['floor'];same_alignment=a['alignment']==b['alignment']
        if not (same_floor or same_alignment):continue
        if same_floor:assert a['tree_sha256']==b['tree_sha256']
        for level in range(4):
            x,y=a['ancestors'][level],b['ancestors'][level]
            rows.append(dict(comparison='alignment' if same_floor else 'branch_floor',job_a=a['job_id'],job_b=b['job_id'],level=level,length_a=len(x),length_b=len(y),edit_distance=Levenshtein.distance(x,y),same_sequence=x==y,assumed_root=level==3))
    assert len(rows)==16
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    with (out/'comparisons.tsv').open('w') as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    (out/'diagnostic_sequences.json').write_text(json.dumps(records,indent=2)+'\n')
    r=dict(status='complete_four_alternative_diagnostic_comparison',comparisons=16,changed=sum(not r['same_sequence'] for r in rows),maximum_edit_distance=max(r['edit_distance'] for r in rows),pins=pins,artifacts={p.name:sha(p) for p in out.iterdir()},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')

if __name__=='__main__':main()
