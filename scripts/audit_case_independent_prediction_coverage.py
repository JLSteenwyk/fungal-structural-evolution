#!/usr/bin/env python3
"""Audit exact sequence coverage in frozen predictors and prepare accession-specific controls."""
import csv,hashlib,json,re
from collections import defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    base=Path('results');sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name;assert sha(p)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    snapshot=base/'structural_comparisons/whole-domain-case-pae-inputs-20260927-v1'
    models=json.loads(checked(snapshot,'model_provenance.json').read_text());wanted={m['sequence_sha256'] for m in models};assert len(models)==len(wanted)==39
    inventories=[base/'structures/esmfold-all-completed-20260922-v1']+sorted((base/'structures').glob('esmfold-experimental-controls*'))
    found=defaultdict(list);census=[]
    for root in inventories:
        if not (root/'inventory.jsonl').is_file():continue
        count=0
        for line in checked(root,'inventory.jsonl').open():
            for r in json.loads(line).get('models',[]):
                count+=1
                if r['sequence_sha256'] in wanted:
                    assert sha(r['path'])==r['sha256'];found[r['sequence_sha256']].append(r)
        census.append(dict(inventory=str(root),models=count))
    experimental=defaultdict(list)
    with checked(base/'experimental_structures/sequence-screen-full-v1','sequence_correspondence.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            for digest in {r['model_sequence_sha256'],r['entity_sequence_sha256']} & wanted:experimental[digest].append(r)
    sequences={}
    for name in ['duplication-alignment-inputs-20260926-v1','duplication-reference-alignment-inputs-20260926-v1']:
        for line in checked(base/'structural_comparisons'/name,'inputs.jsonl').open():
            r=json.loads(line)
            if r['mask']=='full':
                digest=hashlib.sha256(r['sequence'].encode()).hexdigest()
                if digest in wanted:sequences[digest]=r['sequence']
    assert set(sequences)==wanted
    requests=[];coverage=[]
    for m in models:
        accession=re.fullmatch(r'AF-(.+)-F1',m['model_id']);assert accession
        requests.append(dict(m,uniprot_accession=accession.group(1)))
        coverage.append(dict(model_id=m['model_id'],version=m['version'],sequence_sha256=m['sequence_sha256'],length=m['length'],existing_esm_models=len(found[m['sequence_sha256']]),existing_experimental_screen_rows=len(experimental[m['sequence_sha256']])))
    out=base/'structural_comparisons/case-independent-control-inputs-20260927-v1';out.mkdir(exist_ok=False)
    (out/'model_provenance.json').write_text(json.dumps(requests,indent=2)+'\n')
    with (out/'coverage.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(coverage[0]),delimiter='\t');w.writeheader();w.writerows(coverage)
    missing={m['sequence_sha256']:sequences[m['sequence_sha256']] for m in models if not found[m['sequence_sha256']]}
    (out/'missing_independent_sequences.faa').write_text(''.join('>S'+k+'\n'+v+'\n' for k,v in sorted(missing.items())))
    from Bio import SeqIO
    assert {r.id[1:]:str(r.seq) for r in SeqIO.parse(out/'missing_independent_sequences.faa','fasta')}==missing
    with (out/'coverage.tsv').open() as f:assert list(csv.DictReader(f,delimiter='\t'))==[{k:str(v) for k,v in r.items()} for r in coverage]
    for path,digest in sources.items():assert sha(path)==digest
    receipt=dict(status='complete_existing_independent_case_coverage_audit',source_hashes=sources,script_sha256=sha(__file__),models=39,inventory_census=census,case_sequences_with_esm=sum(bool(found[k]) for k in wanted),case_sequences_in_existing_experimental_screen=sum(bool(experimental[k]) for k in wanted),missing_independent_sequences=len(missing),missing_sequence_residues=sum(map(len,missing.values())),artifacts={p.name:sha(p) for p in out.iterdir()},resources=dict(rcsb_http_workers=1,rcsb_requests=1,planning_minutes=[1,5],memory_gib=2,output_gib=.1,paid_resources=False),scope='Exact sequence hash audit of declared frozen local inventories only; absence here is not absence of experimental homologs. FASTA prepares future independent predictions but does not authorize or start paused GPUs. RCSB accession candidate nomination remains separate from sequence and coordinate validation.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
