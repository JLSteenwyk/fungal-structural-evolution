#!/usr/bin/env python3
"""Bind every local fitted alignment sequence to the audited FCS marker map."""
import csv
import hashlib
import json
from pathlib import Path
import re
from Bio import SeqIO


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    fits=Path('results/cds/local-mg94-diagnostics-20260927-v1')
    mapping=Path('results/qc/fcs-cds-overlap-v1')
    proof=Path('metadata/fcs_cds_overlap_audit_receipt.json')
    out=Path('results/qc/local-codon-fcs-exposure-20260927-v1')
    assert not out.exists()
    mr=json.loads((mapping/'receipt.json').read_text())
    pr=json.loads(proof.read_text())
    assert pr['status']=='passed_full_region_cds_intersection_and_protein_marker_join_audit'
    assert pr['mapping_receipt_sha256']==sha(mapping/'receipt.json')
    assert mr['artifacts']['marker_overlap_review.tsv']==sha(mapping/'marker_overlap_review.tsv')
    with (mapping/'marker_overlap_review.tsv').open() as f: raw=list(csv.DictReader(f,delimiter='\t'))
    index={(r['marker'],r['taxon_id']):r for r in raw}
    assert len(index)==len(raw)==59840
    receipt=json.loads((fits/'receipt.json').read_text())
    allrows=[];positive=[];characters=0;seen=set()
    for item in receipt['case_receipts']:
        case=item['case_id'];folder=fits/case
        assert case not in seen;seen.add(case)
        assert sha(folder/'receipt.json')==item['receipt_sha256']
        r=json.loads((folder/'receipt.json').read_text())
        assert sha(folder/'config.json')==r['config_sha256']
        config=json.loads((folder/'config.json').read_text())
        command=config['command'];alignment=Path(command[command.index('--alignment')+1])
        assert sha(alignment)==config['alignment_sha256']
        assert sha(folder/'fit.bf')==r['artifacts']['fit.bf']
        records=list(SeqIO.parse(alignment,'fasta'));seqs={s.id:str(s.seq) for s in records}
        assert len(seqs)==len(records)
        # Independently parse the characters actually embedded in the saved fit.
        text=(folder/'fit.bf').read_text()
        blocks=re.findall(r'(?im)^MATRIX\s*\n(.*?);',text,re.S|re.M)
        assert len(blocks)==1
        embedded={}
        for line in blocks[0].splitlines():
            if not line.strip():continue
            match=re.fullmatch(r"\s*'([^']+)'\s+([A-Za-z?\-]+)\s*",line)
            assert match,line
            taxon,seq=match.groups();assert taxon not in embedded;embedded[taxon]=seq
        assert embedded==seqs
        marker=case.split('__')[1]
        for taxon,seq in sorted(seqs.items()):
            source=index[marker,taxon]
            row=dict(case_id=case,marker=marker,taxon_id=taxon,protein_id=source['protein_id'],fcs_actions=source['fcs_actions'],fcs_overlap_status=source['fcs_overlap_status'],alignment_sha256=config['alignment_sha256'],fit_receipt_sha256=item['receipt_sha256'])
            allrows.append(row);characters+=len(seq)
            if {'EXCLUDE','FIX','TRIM'} & set(source['fcs_actions'].split(';')):positive.append(row)
    assert len(seen)==1632
    out.mkdir(parents=True)
    for name,rows in [('all_case_taxon_exposure.tsv',allrows),('codon_case_exposure.tsv',positive)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,list(allrows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result=dict(status='passed_full_local_fit_alignment_fcs_join',cases=len(seen),case_taxon_rows=len(allrows),embedded_alignment_characters_checked=characters,fcs_exposed_rows=len(positive),fcs_exposed_cases=len({r['case_id'] for r in positive}),source_fit_receipt_sha256=sha(fits/'receipt.json'),mapping_receipt_sha256=sha(mapping/'receipt.json'),mapping_audit_sha256=sha(proof),script_sha256=sha(Path(__file__)),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Every fitted alignment and embedded sequence checked against exact source hashes and the full audited marker map. Exposure means EXCLUDE/FIX/TRIM overlap; zero exposure does not establish absence of contamination or biological eligibility.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
