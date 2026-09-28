#!/usr/bin/env python3
"""Independent readback and paired sensitivity summaries of completed Historian jobs."""
import argparse, collections, csv, hashlib, itertools, json, math, re, time
from pathlib import Path
from Bio import Phylo, SeqIO
from rapidfuzz.distance import Levenshtein
from prepare_case_ancestral_neighborhoods import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    pp=Path('metadata/historian_famsa_memory_retry_plan_20260927.json');plan=json.loads(pp.read_text())
    inputs=json.loads(Path(plan['input_receipt']).read_text());jobs=inputs['jobs']
    for p,h in plan['pins'].items(): assert sha(p)==h,p
    root=Path(plan['output']);out=args.output;out.mkdir(parents=True,exist_ok=False)
    mapping_path=Path('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv')
    mappings=list(csv.DictReader(mapping_path.open(),delimiter='\t'))
    rows=[];ancestors={};pins={str(pp):sha(pp),str(mapping_path):sha(mapping_path),str(Path(__file__)):sha(__file__)}
    for job in jobs:
        folder=root/job['job_id'];rp=folder/'receipt.json'
        row=dict(job_id=job['job_id'],input_id=job['input_id'],status='pending_receipt')
        if not rp.exists(): rows.append(row);continue
        r=json.loads(rp.read_text());assert r['plan_sha256']==sha(pp) and r['job']==job
        pins[str(rp)]=sha(rp)
        for name,h in r['artifacts'].items():assert sha(folder/name)==h,(folder,name)
        row.update(producer_status=r['status'],exit_code=r['exit_code'],elapsed_seconds=r['elapsed_seconds'],peak_sampled_rss_bytes=r['peak_sampled_rss_bytes'])
        if r['exit_code']!=0:
            row['status']='recorded_unsuccessful_execution';rows.append(row);continue
        tree=Phylo.read(job['tree'],'newick')
        assert sha(job['tree'])==job['tree_sha256'] and sha(job['alignment'])==job['alignment_sha256']
        data=json.loads((folder/'reconstruction.json').read_text());seqs=data['rowData']
        edges=data['branches'];children=collections.defaultdict(list);parent={};length={}
        for a,b,v in edges:
            assert b not in parent and math.isfinite(v) and v>=0
            parent[b]=a;children[a].append(b);length[b]=v
        names=set(parent)|set(parent.values());roots=names-set(parent)
        assert roots=={data['root']}=={tree.root.name}
        assert names==set(seqs)=={n.name for n in tree.find_clades()}
        seen=set();stack=[data['root']]
        while stack:
            n=stack.pop();assert n not in seen;seen.add(n);stack.extend(children[n])
        assert seen==names
        assert all(len(children[n]) in (0,2) for n in names)
        branch_error=0.
        for n in tree.find_clades():
            assert set(children[n.name])=={c.name for c in n.clades}
            for c in n.clades:
                error=abs(length[c.name]-c.branch_length);branch_error=max(branch_error,error)
                assert error<=1e-5*max(1.,abs(c.branch_length))
        assert all(isinstance(s,str) and set(s.upper())<=set('ARNDCQEGHILKMFPSTWYVX-') for s in seqs.values())
        widths={len(s) for s in seqs.values()};assert len(widths)==1
        observed={r.id:str(r.seq).replace('-','').upper() for r in SeqIO.parse(job['alignment'],'fasta')}
        leaves={n for n in names if not children[n]};assert leaves==set(observed) and len(leaves)==job['proteins']
        unknown=0
        for n,s in observed.items():
            result=seqs[n].replace('-','').upper();assert len(s)==len(result)
            assert all(a==b or a=='X' for a,b in zip(s,result));unknown+=s.count('X')
        log=(folder/'run.log').read_text()
        likelihood=re.findall(r'Final Forward log-likelihood is ([-+0-9.eE]+), final alignment log-likelihood is ([-+0-9.eE]+)',log)
        assert len(likelihood)==1 and all(math.isfinite(float(v)) for v in likelihood[0])
        dataset='whole' if '-whole-' in job['input_id'] else 'domain'
        selected=[m for m in mappings if m['guide']=='profile' and m['family']==job['family'] and m['dataset']==dataset]
        assert len(selected)==4
        levels={}
        for m in selected:
            node=json.loads(m['output_nodes_json']);assert len(node)==1 and m['unique_output_node']=='1'
            node=node[0];assert node not in job['artificial_nodes']
            actual=set();todo=[node]
            while todo:
                n=todo.pop()
                if children[n]:todo.extend(children[n])
                else:actual.add(n)
            assert actual==set(json.loads(m['retained_set_json']))
            levels[int(m['level'])]=seqs[node].replace('-','').upper()
        ancestors[job['job_id']]=levels
        row.update(status='independent_tip_tree_and_candidate_readback_passed',tips=len(leaves),columns=next(iter(widths)),maximum_branch_serialization_error=branch_error,unknown_extant_residues=unknown,guide_relaxations=log.count('Zero forward likelihood'),reported_forward_log_likelihood=float(likelihood[0][0]),reported_alignment_log_likelihood=float(likelihood[0][1]))
        rows.append(row)
    paired=[]
    for input_id in sorted({j['input_id'] for j in jobs}):
        available=[j for j in jobs if j['input_id']==input_id and j['job_id'] in ancestors]
        for a,b in itertools.combinations(available,2):
            # Only compare one factor at a time: floor or star resolution.
            fa=a['job_id'].split('-floor')[1].split('-resolution')[0];fb=b['job_id'].split('-floor')[1].split('-resolution')[0]
            ra=a['job_id'].split('-resolution')[1];rb=b['job_id'].split('-resolution')[1]
            if fa!=fb and ra!=rb:continue
            comparison='branch_floor' if fa!=fb else 'star_resolution'
            for level in range(4):
                x=ancestors[a['job_id']][level];y=ancestors[b['job_id']][level]
                paired.append(dict(input_id=input_id,job_a=a['job_id'],job_b=b['job_id'],comparison=comparison,level=level,length_a=len(x),length_b=len(y),edit_distance=Levenshtein.distance(x,y),same_sequence=x==y,root_assumption=level==3))
    (out/'dispositions.json').write_text(json.dumps(rows,indent=2)+'\n')
    if paired:
        with (out/'paired_sequences.tsv').open('w') as h:
            w=csv.DictWriter(h,fieldnames=list(paired[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(paired)
    counts=dict(collections.Counter(r['status'] for r in rows))
    result=dict(status='partial_snapshot' if counts.get('pending_receipt',0) else 'all_job_dispositions_audited',expected_jobs=len(jobs),dispositions=counts,paired_node_comparisons=len(paired),changed_node_comparisons=sum(not r['same_sequence'] for r in paired),maximum_edit_distance=max((r['edit_distance'] for r in paired),default=None),observed_unix=time.time(),pins=pins,artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scope='Independent output graph, all known extant residues, tip lengths and candidate descendant sets checked; no numerical likelihood replay. Paired ungapped sequence edit distances are computational sensitivity, not posterior uncertainty or biological event counts. Assumed root retained explicitly. All frozen retry dispositions retained including pending/failures; original grid dispositions remain separate.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts','scope']}))

if __name__=='__main__':main()
