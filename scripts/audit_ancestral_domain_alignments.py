#!/usr/bin/env python3
"""Verify 52 domain alignments and map every residue to full-protein coordinates."""
import argparse,csv,gzip,itertools,json,subprocess,time
from pathlib import Path
import numpy as np
import psutil
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha,read


def fasta(path):
    records=list(SeqIO.parse(path,'fasta'));data={r.id:str(r.seq) for r in records}
    assert len(records)==len(data)
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());digest=sha(args.plan)
    def verify():
        assert sha(args.plan)==digest
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    identity=plan['producer']
    while True:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,SubState,MainPID,Result,ExecMainStatus'],text=True).splitlines())
        if state['ActiveState']=='inactive':
            assert state['Result']=='success' and state['ExecMainStatus']=='0';break
        assert state['ActiveState'] in ['active','activating','deactivating'],state
        if state['SubState']=='running':
            p=psutil.Process(identity['pid']);assert int(state['MainPID'])==p.pid and p.create_time()==identity['created'] and p.cmdline()==identity['cmdline']
        print('waiting_for_all_52_domain_alignments',flush=True);time.sleep(30)
    verify();root=Path(plan['runs']);input_root=Path(plan['inputs'])
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_52_ancestral_domain_alignments_pending_independent_audit'
    assert receipt['plan_sha256']==sha(plan['producer_plan'])
    input_receipt=json.loads((input_root/'receipt.json').read_text())
    cases={r['family']+'-'+boundary:dict(r,boundary=boundary) for r in read(plan['cases']) for boundary in ['alignment','envelope']}
    assert len(cases)==26 and len(receipt['results'])==52
    bykey={(r['input_set'],r['method']):r for r in receipt['results']}
    assert set(bykey)=={(f,m) for f in cases for m in ['mafft','famsa']}
    summaries=[];columns=[];focal_pairs=[];audits=[];position_records=[]
    dispositions=read(input_root/'protein_dispositions.tsv')
    assert sha(input_root/'protein_dispositions.tsv')==input_receipt['artifacts']['protein_dispositions.tsv']
    bounds={(r['family']+'-'+r['boundary'],r['gene']):r for r in dispositions if r['status']=='extracted_single_target_hit'}
    full_root=Path(plan['full_proteins']);full_receipt=json.loads((full_root/'receipt.json').read_text())
    for family,case in cases.items():
        original_path=input_root/(family+'.faa');assert sha(original_path)==input_receipt['artifacts'][original_path.name]
        original=fasta(original_path);genes=sorted(original);maps={};pairs={}
        full_path=full_root/(case['family']+'.faa');assert sha(full_path)==full_receipt['artifacts'][full_path.name];full=fasta(full_path)
        assert set(genes)=={g for key,g in bounds if key==family}
        for method in ['mafft','famsa']:
            folder=root/(family+'-'+method);rp=folder/'receipt.json';r=json.loads(rp.read_text());assert r==bykey[family,method]
            assert r['status']=='complete_exact_sequence_preserving_alignment' and r['input_sha256']==sha(original_path)
            for name,h in r['artifacts'].items():assert sha(folder/name)==h
            aligned=fasta(folder/'alignment.faa');assert set(aligned)==set(original)
            assert all(aligned[g].replace('-','')==original[g] for g in genes)
            assert len({len(s) for s in aligned.values()})==1
            array=np.array([list(aligned[g]) for g in genes]);occupied=array!='-'
            positions=np.where(occupied,np.cumsum(occupied,axis=1),0)
            # Independently enumerate each sequence's nongap columns to check
            # every inferred sequence position, not merely column totals.
            for row,g in enumerate(genes):
                nongaps=[i for i,x in enumerate(aligned[g]) if x!='-']
                assert len(nongaps)==len(original[g])
                assert [int(positions[row,i]) for i in nongaps]==list(range(1,len(original[g])+1))
                assert ''.join(aligned[g][i] for i in nongaps)==original[g]
            counts=occupied.sum(axis=0);canonical=np.isin(array,list('ACDEFGHIKLMNPQRSTVWY')).sum(axis=0)
            source=read(folder/'column_coverage.tsv');assert len(source)==array.shape[1]==r['columns']
            assert array.shape[0]==r['proteins']
            for i,row in enumerate(source):
                assert int(row['column'])==i+1 and int(row['observed_residues'])==counts[i] and int(row['canonical_residues'])==canonical[i] and int(row['total_sequences'])==len(genes)
            for threshold in [50,70,90]:assert r[f'columns_at_least_{threshold}pct']==int((100*counts>=threshold*len(genes)).sum())
            for row,g in enumerate(genes):
                start,end=int(bounds[family,g]['start']),int(bounds[family,g]['end'])
                assert original[g]==full[g][start-1:end]
                positions[row]=np.where(occupied[row],positions[row]+start-1,0)
                for col in np.flatnonzero(occupied[row]):
                    pos=int(positions[row,col]);assert aligned[g][col]==full[g][pos-1]
                    position_records.append(dict(family=case['family'],boundary=case['boundary'],method=method,gene=g,column=int(col)+1,protein_position=pos,residue=aligned[g][col]))
            signatures=[tuple(int(x) for x in positions[:,i]) for i in range(array.shape[1])]
            assert len(set(signatures))==len(signatures) and all(any(s) for s in signatures)
            maps[method]={s:i+1 for i,s in enumerate(signatures)}
            ai,bi=genes.index(case['gene_a']),genes.index(case['gene_b'])
            pairs[method]={(int(positions[ai,i]),int(positions[bi,i])) for i in range(array.shape[1]) if positions[ai,i] and positions[bi,i]}
            audits.append(dict(family=family,method=method,proteins=len(genes),columns=len(signatures),residues=int(counts.sum()),unknown_residues=int((occupied.sum()-canonical.sum())),source_receipt_sha256=sha(rp)))
        a,b=maps['mafft'],maps['famsa'];common=set(a)&set(b)
        for method,other in [('mafft','famsa'),('famsa','mafft')]:
            for signature,col in sorted(maps[method].items(),key=lambda x:x[1]):
                columns.append(dict(family=family,method=method,column=col,observed_sequences=sum(x>0 for x in signature),total_sequences=len(genes),matching_other_column=maps[other].get(signature,''),exact_full_column_match=int(signature in maps[other])))
        ap,bp=pairs['mafft'],pairs['famsa']
        for pair in sorted(ap|bp):focal_pairs.append(dict(family=family,gene_a=case['gene_a'],gene_b=case['gene_b'],position_a=pair[0],position_b=pair[1],mafft=int(pair in ap),famsa=int(pair in bp)))
        summaries.append(dict(family=family,proteins=len(genes),mafft_columns=len(a),famsa_columns=len(b),exact_shared_columns=len(common),mafft_shared_fraction=len(common)/len(a),famsa_shared_fraction=len(common)/len(b),focal_mafft_pairs=len(ap),focal_famsa_pairs=len(bp),focal_shared_pairs=len(ap&bp),focal_union_pairs=len(ap|bp)))
    for filename,rows in [('alignment_audit.tsv',audits),('method_comparison.tsv',summaries),('column_correspondence.tsv',columns),('focal_pair_correspondence.tsv',focal_pairs)]:
        with (out/filename).open('w') as handle:
            w=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
        assert read(out/filename)==[{k:str(v) for k,v in row.items()} for row in rows]
    with gzip.open(out/'protein_position_maps.tsv.gz','wt') as handle:
        w=csv.DictWriter(handle,list(position_records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(position_records)
    with gzip.open(out/'protein_position_maps.tsv.gz','rt') as handle:
        count=0
        for observed,expected in itertools.zip_longest(csv.DictReader(handle,delimiter='\t'),position_records):
            assert observed is not None and expected is not None
            assert observed=={k:str(v) for k,v in expected.items()};count+=1
        assert count==len(position_records)
    verify()
    proof=dict(status='complete_52_ancestral_domain_alignment_readbacks_and_protein_correspondences',plan_sha256=digest,source_receipt_sha256=sha(root/'receipt.json'),terminal_state=state,alignments=len(audits),family_boundary_sets=len(summaries),protein_position_records=len(position_records),column_records=len(columns),focal_pair_records=len(focal_pairs),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Every domain input residue/copy, occupancy column and original protein coordinate checked; all position maps exported. Full-column and focal-pair comparisons use absolute protein positions within each boundary definition. Exact column matching is stringent, not a homology probability. Export family IDs include the boundary suffix. No cross-boundary inference, ancestral sequence, trimming or selection claim.')
    (out/'receipt.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof),flush=True)


if __name__=='__main__':main()
