#!/usr/bin/env python3
"""Join per-site minimum character changes on audited topologies to exposure annotations."""
import argparse,csv,gzip,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import Phylo,SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def minimum_changes(tree,sequences,alphabet='ACDEFGHIKLMNPQRSTVWY'):
    tips=[x.name for x in tree.get_terminals()]
    if len(tips)!=len(set(tips)) or set(tips)!=set(sequences):raise ValueError('Tree/sequence taxa differ')
    widths={len(s) for s in sequences.values()}
    if len(widths)!=1:raise ValueError('Unequal sequence lengths')
    width=widths.pop();index={x:i for i,x in enumerate(alphabet)};cost={};infinity=len(tips)+1
    if len(index)!=len(alphabet):raise ValueError('Duplicate alphabet states')
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            v=np.full((width,len(alphabet)),infinity,dtype=np.int32)
            for j,state in enumerate(sequences[node.name]):
                if state=='?':v[j,:]=0
                elif state in index:v[j,index[state]]=0
                else:raise ValueError('Unknown observed character: '+state)
        else:
            v=np.zeros((width,len(alphabet)),dtype=np.int32)
            for child in node.clades:
                q=cost.pop(child)
                v+=np.minimum(q,q.min(axis=1)[:,None]+1)
        cost[node]=v
    return cost[tree.root].min(axis=1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','fits','audit','exposure','projection','review','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable output')
    inputs=checked_receipt(a.inputs);audit=checked_receipt(a.audit);exposure=checked_receipt(a.exposure);projection=checked_receipt(a.projection)
    fits=json.loads((a.fits/'receipt.json').read_text());config=json.loads((a.fits/'config.json').read_text())
    if audit['fit_receipt_sha256']!=sha(a.fits/'receipt.json') or fits['config_sha256']!=sha(a.fits/'config.json') or config['input_receipt_sha256']!=sha(a.inputs/'receipt.json') or exposure['source_receipt_sha256']!=sha(a.projection/'receipt.json') or projection['source_receipts']['inputs']!=sha(a.inputs/'receipt.json'):
        raise ValueError('Input lineage differs')
    caveats={r['marker']:r['status'] for r in json.loads(a.review.read_text())['records']}
    sites=defaultdict(dict);n=0
    with gzip.open(a.exposure/'normalized_paired_sites.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=r['marker'],int(r['paired_column_1based']);taxon=r['taxon_id']
            if taxon in sites[key]:raise ValueError('Duplicate site observation')
            sites[key][taxon]=(r['amino_acid'],r['3di_state'],r['model_id'],r['rsa_tien2013_theoretical'],r['rsa_miller1987']);n+=1
    if n!=exposure['rows']:raise ValueError('Exposure row count differs')
    fitted={r['marker'] for r in fits['results']};ready=[r for r in read_table(a.inputs/'marker_summary.tsv') if r['status']=='ready_for_inference']
    if fitted!={r['marker'] for r in ready} or len(ready)!=inputs['ready_markers']:raise ValueError('Fitted marker set differs')
    results=[];tree_pins={};totalobserved=0
    for m in ready:
        marker=m['marker'];folder=a.inputs/marker;path=a.fits/marker/'aa.treefile';fr=json.loads((a.fits/marker/'aa.receipt.json').read_text())
        if sha(path)!=fr['artifacts'][path.name]:raise ValueError('Changed fitted topology')
        tree_pins[marker]=sha(path);tree=Phylo.read(path,'newick')
        aa={r.id:str(r.seq) for r in SeqIO.parse(folder/'aa.faa','fasta')};di={r.id:str(r.seq) for r in SeqIO.parse(folder/'3di.faa','fasta')}
        scores=[minimum_changes(tree,x) for x in [aa,di]];cols=read_table(folder/'columns.tsv')
        if len(cols)!=len(scores[0]) or len(cols)!=len(scores[1]):raise ValueError('Column grid differs')
        for j,col in enumerate(cols):
            obs=sites.pop((marker,j+1));expected={t for t,s in aa.items() if s[j]!='?'}
            if set(obs)!=expected or expected!={t for t,s in di.items() if s[j]!='?'}:raise ValueError('Observed site cohort differs')
            if any(v[0]!=aa[t][j] or v[1]!=di[t][j] for t,v in obs.items()):raise ValueError('Character states differ')
            distinct=[len({v[k] for v in obs.values()}) for k in [0,1]]
            if any(not distinct[k]-1<=scores[k][j]<=len(obs)-1 for k in [0,1]):raise ValueError('Parsimony bounds failed')
            row={'marker':marker,'paired_column_1based':j+1,'matrix_column_1based':col['matrix_column_1based'],'observed_taxa':len(obs),'unique_models':len({v[2] for v in obs.values()}),'aa_distinct_states':distinct[0],'3di_distinct_states':distinct[1],'aa_minimum_changes':int(scores[0][j]),'3di_minimum_changes':int(scores[1][j]),'marker_review_status':caveats.get(marker,'no_current_copy_review_flag')}
            for label,k in [('tien2013',3),('miller1987',4)]:
                values=[float(v[k]) for v in obs.values() if v[k]!='']
                row[label+'_normalized_taxa']=len(values)
                for stat,quantile in [('q25',.25),('median',.5),('q75',.75)]:row[label+'_rsa_'+stat]=float(np.quantile(values,quantile)) if values else ''
            results.append(row);totalobserved+=len(obs)
        print(marker,'sites',len(cols),flush=True)
    if sites or totalobserved!=n:raise ValueError('Incomplete exposure projection')
    a.output.mkdir(parents=True);write_table(a.output/'site_parsimony_exposure.tsv',results)
    r={'status':'complete_site_parsimony_exposure_diagnostic','markers':len(ready),'sites':len(results),'taxon_site_observations':n,'source_receipts':{k:sha(getattr(a,k)/'receipt.json') for k in ['inputs','fits','audit','exposure','projection']},'marker_review_sha256':sha(a.review),'tree_sha256':tree_pins,'script_sha256':sha(Path(__file__)),'interpretation':'Unit-cost minimum character changes per site on the fitted AA topology, with unknown states unconstrained. Uses topology but ignores branch lengths, time, unequal substitution costs, homoplasy beyond its minimum and topology uncertainty. AA/3Di counts are alphabet-dependent lower bounds, not rates or selection tests. Exposure quantiles describe extant observed taxa and are not ancestral values or phylogenetically adjusted estimates. No particular branch/state assignment or coupling significance is inferred.','artifacts':{'site_parsimony_exposure.tsv':sha(a.output/'site_parsimony_exposure.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='tree_sha256'},indent=2))

if __name__=='__main__':main()
