#!/usr/bin/env python3
"""Recover every level-three case sequence with independent proteome matching."""
import csv,hashlib,json
from collections import defaultdict,Counter
from pathlib import Path
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha,read


def main():
    root=Path('results/ancestral/case-neighborhoods-20260927-v1')
    proof=Path('metadata/case_ancestral_neighborhoods_completed_20260927.json')
    audit=json.loads(proof.read_text());assert audit['status']=='passed_full_case_ancestral_neighborhood_readback'
    rp=root/'receipt.json';assert audit['source_receipt_sha256']==sha(rp)
    receipt=json.loads(rp.read_text())
    for name in ['neighborhoods.tsv','descendants.tsv']:
        assert sha(root/name)==receipt['artifacts'][name]
    selected=[r for r in read(root/'descendants.tsv') if r['level']=='3']
    wanted={r['gene'] for r in selected};families=defaultdict(set)
    for row in selected:families[row['family']].add(row['gene'])
    assert len(families)==13
    stage=Path('results/orthology/full-reconciliation-inputs-v1')
    staging=json.loads((stage/'receipt.json').read_text())
    manifest=stage/'profile/copied_files.tsv'
    assert sha(manifest)==next(x['manifest_sha256'] for x in staging['guides'] if x['guide']=='profile')
    expected={x['relative_path']:x['sha256'] for x in read(manifest)}
    wd=stage/'profile/Source/WorkingDirectory'
    pins={str(p):sha(p) for p in [proof,rp,root/'descendants.tsv',root/'neighborhoods.tsv',stage/'receipt.json',manifest]}
    def check_staged(name):
        path=wd/name;h=sha(path);assert h==expected['Source/WorkingDirectory/'+name];pins[str(path)]=h;return path
    taxa={}
    for line in check_staged('SpeciesIDs.txt').read_text().splitlines():
        number,name=line.split(': ',1);taxa[number]=Path(name).stem
    mapping={};by_species=defaultdict(dict)
    with check_staged('SequenceIDs.txt').open() as handle:
        for line in handle:
            native,accession=line.rstrip().split(': ',1)
            sid=native.split('_',1)[0];gene=taxa[sid]+'_'+accession
            if gene in wanted:
                assert gene not in mapping;mapping[gene]=native;by_species[sid][native]=gene
    assert set(mapping)==wanted
    sequences={};checks=[]
    for sid,targets in sorted(by_species.items(),key=lambda x:int(x[0])):
        path=check_staged('Species'+sid+'.fa')
        for record in SeqIO.parse(path,'fasta'):
            if record.id in targets:
                gene=targets[record.id];assert gene not in sequences;sequences[gene]=str(record.seq)
        qc=Path('data/qc_proteomes')/(taxa[sid]+'.faa');pins[str(qc)]=sha(qc)
        # Check identifiers and complete sequence strings against the QC proteome,
        # without using the native numerical mapping to parse its headers.
        observed=set()
        for record in SeqIO.parse(qc,'fasta'):
            gene=taxa[sid]+'_'+record.id
            if gene in targets.values():
                assert gene not in observed and str(record.seq)==sequences[gene];observed.add(gene)
        assert observed==set(targets.values())
        print(taxa[sid],len(observed),'exact_sequences_checked',flush=True)
    assert set(sequences)==wanted
    output=Path('results/ancestral/case-sequences-20260927-v1');output.mkdir(exist_ok=False)
    canonical=set('ACDEFGHIKLMNPQRSTVWY');rows=[];summaries=[]
    for family,genes in sorted(families.items()):
        path=output/(family+'.faa');path.write_text(''.join('>'+g+'\n'+sequences[g]+'\n' for g in sorted(genes)))
        assert {r.id:str(r.seq) for r in SeqIO.parse(path,'fasta')}=={g:sequences[g] for g in genes}
        for gene in sorted(genes):
            seq=sequences[gene];extra=Counter(x for x in seq if x not in canonical)
            rows.append(dict(family=family,gene=gene,native_id=mapping[gene],taxon_id=gene.split('_',1)[0],length=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),noncanonical_symbols=json.dumps(dict(sorted(extra.items())),sort_keys=True),noncanonical_residues=sum(extra.values()),internal_stops=seq[:-1].count('*'),terminal_stop=int(seq.endswith('*'))))
        lengths=[len(sequences[g]) for g in genes]
        summaries.append(dict(family=family,proteins=len(genes),taxa=len({g.split('_',1)[0] for g in genes}),residues=sum(lengths),minimum_length=min(lengths),maximum_length=max(lengths),unique_sequences=len({sequences[g] for g in genes}),noncanonical_proteins=sum(bool(set(sequences[g])-canonical) for g in genes)))
    for filename,records in [('sequence_inventory.tsv',rows),('family_summary.tsv',summaries)]:
        with (output/filename).open('w') as f:
            w=csv.DictWriter(f,list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
        assert read(output/filename)==[{k:str(v) for k,v in r.items()} for r in records]
    for path,h in pins.items():assert sha(path)==h
    result=dict(status='complete_exact_ancestral_case_extant_sequences',families=len(families),unique_genes=len(wanted),family_gene_records=len(rows),taxa=len(by_species),residues=sum(map(len,sequences.values())),noncanonical_proteins=sum(bool(set(s)-canonical) for s in sequences.values()),source_hashes=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in output.iterdir()},scope='All level-three descendant genes retained across both guides, without sequence deduplication, trimming or residue replacement. Every selected sequence exactly matches both checksum-bound native reconciliation input and taxon-specific QC proteome. Candidate clade inputs only; outside-clade context, alignment assessment and ancestral inference remain pending.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
