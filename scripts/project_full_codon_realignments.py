#!/usr/bin/env python3
"""Project all verified local alignments and distinguish correspondence from mask loss."""
import argparse,csv,json,hashlib,itertools
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from codon_realignment_projection import project,correspondence,compare

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fasta(p):return {r.id:str(r.seq).upper() for r in SeqIO.parse(p,'fasta')}
def write(path,rows,fields=None):
    with path.open('w') as f:
        w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in p['pins'].items():assert sha(path)==h,path
    verify();inputs=Path(p['inputs']);aligned=Path(p['alignments']);ar=json.loads((aligned/'receipt.json').read_text());proof=json.loads(Path(p['alignment_readback']).read_text());ir=json.loads((inputs/'receipt.json').read_text())
    assert proof['status']=='passed_full_codon_group_protein_alignment_readback' and proof['producer_receipt_sha256']==sha(aligned/'receipt.json')
    assert ar['source_receipt_sha256']==sha(inputs/'receipt.json')
    for name,h in ir['artifacts'].items():assert sha(inputs/name)==h
    with (inputs/'cases.tsv').open() as f:cases=list(csv.DictReader(f,delimiter='\t'))
    receipts={r['case_id']:r['receipt_sha256'] for r in ar['case_receipts']};assert {c['case_id'] for c in cases}==set(receipts) and len(cases)==len(receipts)==1712
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False);summary=[];pairs=[];artifacts={}
    for case in cases:
        cid=case['case_id'];source=inputs/cid;local=aligned/cid;assert sha(local/'receipt.json')==receipts[cid];cr=json.loads((local/'receipt.json').read_text())
        for name,h in cr['artifacts'].items():assert sha(local/name)==h
        aa=fasta(local/'aligned_amino_acids.faa');dna=fasta(source/'full_codons.fna');original=json.loads((source/'original_residue_positions.json').read_text());assert set(aa)==set(original)
        positions,codons,keep=project(aa,dna,int(case['translation_table']));folder=out/cid;folder.mkdir()
        for name,seqs in [('codons.fna',{t:''.join(codons[t][i] for i in keep) for t in aa}),('amino_acids.faa',{t:''.join(aa[t][i] if positions[t][i] else '?' for i in keep) for t in aa})]:
            with (folder/name).open('w') as f:
                for t,s in seqs.items():f.write('>'+t+'\n'+s+'\n')
        write(folder/'columns.tsv',[dict(retained_codon_column_1based=j+1,local_alignment_column_1based=i+1,called_taxa=sum(positions[t][i]>0 for t in aa)) for j,i in enumerate(keep)],['retained_codon_column_1based','local_alignment_column_1based','called_taxa'])
        (folder/'raw_called_residue_positions.json').write_text(json.dumps(positions,separators=(',',':'))+'\n')
        for t,u in itertools.combinations(sorted(aa),2):
            old=correspondence(original[t],original[u]);raw=correspondence(positions[t],positions[u]);selected=correspondence([positions[t][i] for i in keep],[positions[u][i] for i in keep]);stats=compare(old,raw,selected)
            assert stats['original_pairs']==stats['original_preserved_retained']+stats['original_absent_from_raw']+stats['original_preserved_but_filtered']
            pairs.append(dict(case_id=cid,taxon_a=t,taxon_b=u,**stats))
        observed=[sum(positions[t][i]>0 for i in keep) for t in aa];status='passes_existing_coverage_gate' if len(keep)>=100 and min(observed)>0 and len(aa)>=4 else 'below_existing_coverage_gate'
        summary.append(dict(case_id=cid,translation_table=case['translation_table'],taxa=len(aa),original_columns=case['original_codon_columns'],local_alignment_columns=len(next(iter(aa.values()))),retained_codon_columns=len(keep),observed_taxon_codons=sum(observed),minimum_observed_codons=min(observed),marker_copy_caveat=case['marker_copy_caveat'],coverage_disposition=status,selection_eligibility='not_established'))
        for path in folder.iterdir():artifacts[str(path.relative_to(out))]=sha(path)
    write(out/'cases.tsv',summary);write(out/'pair_correspondence.tsv',pairs)
    for name in ['cases.tsv','pair_correspondence.tsv']:artifacts[name]=sha(out/name)
    verify();result=dict(status='complete_full_codon_realignment_projection_pending_readback',plan_sha256=ph,cases=len(summary),pairs=len(pairs),observed_taxon_codons=sum(r['observed_taxon_codons'] for r in summary),dispositions=dict(Counter(r['coverage_disposition'] for r in summary)),source_alignment_receipt_sha256=sha(aligned/'receipt.json'),source_alignment_readback_sha256=sha(p['alignment_readback']),artifacts=artifacts,scope='Exact source codons projected into every local alignment using recorded genetic code; only canonical translated codons count toward integer 80-percent occupancy. Unknowns are whole ??? triplets. All cases retained including coverage losses. Original residue-pair losses partitioned into absent from raw local alignment versus preserved there but filtered; new retained pairs may reflect original masking as well as homology shifts. Independent projection readback remains pending. No selection eligibility or alignment correctness claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))
if __name__=='__main__':main()
