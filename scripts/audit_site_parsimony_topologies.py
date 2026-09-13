#!/usr/bin/env python3
"""Read back all topology-score artifacts and independently check one tree per marker."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
from Bio import Phylo,SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def set_scores(tree,seq):
    nodes=list(tree.find_clades(order='postorder'));alphabet=set('ACDEFGHIKLMNPQRSTVWY');answer=[]
    for j in range(len(next(iter(seq.values())))):
        sets={};scores={}
        for node in nodes:
            if node.is_terminal():
                c=seq[node.name][j];sets[node]=alphabet if c=='?' else {c};scores[node]=0
            else:
                counts=Counter(c for child in node.clades for c in sets[child]);maximum=max(counts.values())
                sets[node]={c for c,n in counts.items() if n==maximum};scores[node]=sum(scores[ch] for ch in node.clades)+len(node.clades)-maximum
        answer.append(scores[tree.root])
    return np.asarray(answer)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['assessment','diagnostic','inputs','fits','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable audit output')
    r=checked_receipt(a.assessment);d=checked_receipt(a.diagnostic);checked_receipt(a.inputs)
    c=json.loads((a.assessment/'config.json').read_text())
    if r['config_sha256']!=sha(a.assessment/'config.json') or any(c['source_receipts'][k]!=sha(getattr(a,k)/'receipt.json') for k in ['diagnostic','inputs','fits']):raise ValueError('Source lineage differs')
    reference={(x['marker'],x['paired_column_1based']):x for x in read_table(a.diagnostic/'site_parsimony_exposure.tsv')}
    combined=read_table(a.assessment/'site_topology_sensitivity.tsv');seen=set();audits=[];countvalues=0;checks=0;varcounts=Counter();expectedtrees=c['bootstrap_trees_per_marker'];combined_by_marker={}
    for row in combined:combined_by_marker.setdefault(row['marker'],[]).append(row)
    if set(combined_by_marker)!=set(r['marker_receipts']) or len(combined)!=d['sites']:raise ValueError('Full site universe differs')
    for marker,pin in sorted(r['marker_receipts'].items()):
        folder=a.assessment/marker
        if sha(folder/'receipt.json')!=pin:raise ValueError('Marker receipt changed')
        mr=checked_receipt(folder);rows=read_table(folder/'site_topology_sensitivity.tsv')
        if mr['config_sha256']!=r['config_sha256'] or mr['bootstrap_trees']!=expectedtrees or rows!=combined_by_marker[marker]:raise ValueError('Marker metadata differs')
        for row in rows:
            key=marker,row['paired_column_1based']
            if key in seen or any(row[k]!=v for k,v in reference[key].items()):raise ValueError('Original site fields changed')
            seen.add(key)
        archive=a.fits/marker/'aa.ufboot'
        if sha(archive)!=mr['bootstrap_tree_sha256']:raise ValueError('Tree archive changed')
        chosen=int(hashlib.sha256(('topology_score_audit_v1|'+marker).encode()).hexdigest(),16)%expectedtrees
        seqs={label:{x.id:str(x.seq) for x in SeqIO.parse(a.inputs/marker/file,'fasta')} for label,file in [('aa','aa.faa'),('3di','3di.faa')]}
        tree_count=0;selected=None
        for index,tree in enumerate(Phylo.parse(archive,'newick')):
            tips=[x.name for x in tree.get_terminals()]
            if len(tips)!=len(set(tips)) or set(tips)!=set(seqs['aa']) or set(tips)!=set(seqs['3di']):raise ValueError('Bootstrap tip universe differs')
            if index==chosen:selected=tree
            tree_count+=1
        if tree_count!=expectedtrees:raise ValueError('Tree count differs')
        with np.load(folder/'bootstrap_site_scores.npz') as arrays:
            if set(arrays.files)!={'aa','3di'}:raise ValueError('Unexpected stored arrays')
            for label in ['aa','3di']:
                v=arrays[label]
                if v.shape!=(expectedtrees,len(rows)) or not np.issubdtype(v.dtype,np.integer):raise ValueError('Score grid differs')
                lower=np.array([int(x[label+'_distinct_states'])-1 for x in rows]);upper=np.array([int(x['observed_taxa'])-1 for x in rows])
                if np.any(v<lower) or np.any(v>upper):raise ValueError('Score bounds failed')
                ordered=np.sort(v,axis=0);original=np.array([int(x[label+'_minimum_changes']) for x in rows])
                for j,row in enumerate(rows):
                    expected={'min':int(ordered[0,j]),'max':int(ordered[-1,j]),'fraction_equal_reference':float(np.count_nonzero(v[:,j]==original[j])/expectedtrees),'distinct_scores':len(set(v[:,j].tolist()))}
                    for name,q in [('p025',.025),('median',.5),('p975',.975)]:
                        position=q*(expectedtrees-1);lo=int(np.floor(position));hi=int(np.ceil(position));fraction=position-lo
                        expected[name]=float(ordered[lo,j])*(1-fraction)+float(ordered[hi,j])*fraction
                    if any(not np.isclose(float(row[label+'_topology_'+k]),value,rtol=1e-12,atol=1e-12) for k,value in expected.items()):raise ValueError('Reported topology summary differs')
                variable=int(np.count_nonzero(ordered[0]!=ordered[-1]));varcounts[label]+=variable
                if variable!=mr[label+'_sites_with_variable_scores']:raise ValueError('Marker variation count differs')
                independent=set_scores(selected,seqs[label])
                if not np.array_equal(independent,v[chosen]):raise ValueError('Independent topology scores differ')
                checks+=len(independent);countvalues+=v.size
        audits.append({'marker':marker,'bootstrap_trees':tree_count,'sites':len(rows),'independent_tree_index_1based':chosen+1,'independent_scores_checked':2*len(rows)})
        print('checked',marker,flush=True)
    if seen!=set(reference) or len(audits)!=r['markers'] or any(varcounts[x]!=r[x+'_sites_with_variable_scores'] for x in ['aa','3di']) or sum(x['bootstrap_trees'] for x in audits)!=r['bootstrap_topologies_evaluated']:raise ValueError('Aggregate counts differ')
    a.output.mkdir(parents=True);write_table(a.output/'marker_audit.tsv',audits)
    out={'status':'passed_full_topology_output_readback','assessment_receipt_sha256':sha(a.assessment/'receipt.json'),'script_sha256':sha(Path(__file__)),'markers':len(audits),'stored_scores_checked':countvalues,'independently_recomputed_scores':checks,'independent_trees_checked':len(audits),'scope':'All saved score dimensions/bounds, artifact pins, topology tip sets, original fields and summaries checked. One hash-selected bootstrap tree per marker independently rescored by a set-membership recurrence; not independent rescoring of all bootstrap trees. No claim of statistical calibration.','artifacts':{'marker_audit.tsv':sha(a.output/'marker_audit.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
