#!/usr/bin/env python3
"""Verify runtime tree identities against every saved diagnostic alignment sample."""
import argparse,csv,json,math,re
from collections import Counter
from io import StringIO
from pathlib import Path
from Bio import Phylo,SeqIO
from prepare_case_ancestral_neighborhoods import sha


def descendant_index(tree):
    index={}
    for n in tree.find_clades():
        key=tuple(sorted(t.name for t in n.get_terminals()));assert key not in index
        index[key]=n
    return index


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    pp=Path('metadata/baliphy_sample_mapping_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h,p
    inputs=json.loads(Path(plan['input_receipt']).read_text());root=Path(plan['output']);out=args.output;out.mkdir(parents=True,exist_ok=False)
    mp=Path('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv');mappings=list(csv.DictReader(mp.open(),delimiter='\t'))
    pins={str(pp):sha(pp),str(mp):sha(mp),str(Path(__file__)):sha(__file__)};rows=[];sample_rows=[]
    for job in inputs['jobs']:
        folder=root/job['job_id'];rp=folder/'receipt.json';row=dict(job_id=job['job_id'],status='pending_receipt')
        if not rp.exists():rows.append(row);continue
        r=json.loads(rp.read_text());assert r['plan_sha256']==sha(pp) and r['job']==job;pins[str(rp)]=sha(rp)
        for p,h in r['artifacts'].items():assert sha(folder/p)==h
        if r.get('exit_code')!=0:row['status']='unsuccessful_diagnostic_retained';rows.append(row);continue
        source=json.loads((folder/'source.json').read_text())
        assert sha(source['receipt'])==source['receipt_sha256']
        original=(Path(source['receipt']).parent/'BAliPhy.Main.hs').read_text();mapped=(folder/'MappedModel.hs').read_text()
        expected=';T.writeFile '+json.dumps(str((folder/'runtime-tree.nwk').resolve()))+' (writeNewick_rooted (addInternalLabels tree))\n'
        assert mapped==original.replace(';mcmcState <- makeMCMCState',expected+';mcmcState <- makeMCMCState')
        source_tree=Phylo.read(job['tree'],'newick');runtime=Phylo.read(folder/'runtime-tree.nwk','newick')
        si=descendant_index(source_tree);ri=descendant_index(runtime);assert set(si)==set(ri)
        for key,n in si.items():
            if n is not source_tree.root:assert abs(n.branch_length-ri[key].branch_length)<1e-10
        labels={n.name for n in runtime.find_clades()};assert len(labels)==len(ri)
        observed={r.id:str(r.seq).replace('-','').upper() for r in SeqIO.parse(job['alignment'],'fasta')}
        assert {n.name for n in runtime.get_terminals()}==set(observed) and len(observed)==job['proteins']
        paths=list(folder.glob('mapping-check-*/C1.P1.fastas'));assert len(paths)==1
        log=list(csv.DictReader((paths[0].parent/'C1.log').open(),delimiter='\t'));assert [int(r['iter']) for r in log]==list(range(21))
        for entry in log:
            for k in ['prior','likelihood','posterior']:assert math.isfinite(float(entry[k]))
            assert abs(float(entry['prior'])+float(entry['likelihood'])-float(entry['posterior']))<1e-7
        text=paths[0].read_text();blocks=re.split(r'^iterations = (\d+)\s*\n',text,flags=re.M);assert not blocks[0].strip()
        assert [int(blocks[i]) for i in range(1,len(blocks),2)]==[0,10,20]
        dataset='whole' if '-whole-' in job['input_id'] else 'domain'
        selected=[m for m in mappings if m['guide']=='profile' and m['family']==job['family'] and m['dataset']==dataset];assert len(selected)==4
        for i in range(1,len(blocks),2):
            iteration=int(blocks[i]);records=list(SeqIO.parse(StringIO(blocks[i+1]),'fasta'));seq={r.id:str(r.seq).upper() for r in records}
            assert len(records)==len(seq)==len(labels) and set(seq)==labels
            assert len({len(s) for s in seq.values()})==1
            assert all(set(s)<=set('ARNDCQEGHILKMFPSTWYVX-') for s in seq.values())
            for name,s in observed.items():
                actual=seq[name].replace('-','');assert len(actual)==len(s) and all(a==b or a=='X' for a,b in zip(s,actual))
            for m in selected:
                key=tuple(sorted(json.loads(m['retained_set_json'])));node=ri[key]
                sample_rows.append(dict(job_id=job['job_id'],iteration=iteration,level=m['level'],source_node=m['source_node'],runtime_node=node.name,ungapped_length=len(seq[node.name].replace('-','')),assumed_root=int(m['level'])==3))
        row.update(status='all_saved_samples_and_candidate_nodes_verified',samples=3,tips=len(observed),runtime_root=runtime.root.name)
        rows.append(row)
    (out/'dispositions.json').write_text(json.dumps(rows,indent=2)+'\n')
    if sample_rows:
        with (out/'candidate_samples.tsv').open('w') as h:
            w=csv.DictWriter(h,fieldnames=list(sample_rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sample_rows)
    counts=dict(Counter(r['status'] for r in rows))
    result=dict(status='partial_mapping_snapshot' if counts.get('pending_receipt') else 'all_mapping_dispositions_checked',dispositions=counts,expected_jobs=324,candidate_samples=len(sample_rows),pins=pins,artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scope='Each saved sample paired with same-process fixed runtime tree; all tips/known residues, rooted clades, branch lengths and candidate identities checked. Initial and early samples retained for integrity testing only, not converged posterior ensembles. Root remains an assumption; no independent likelihood replay.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(counts));print('candidate samples',len(sample_rows))
if __name__=='__main__':main()
