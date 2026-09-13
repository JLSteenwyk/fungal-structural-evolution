#!/usr/bin/env python3
"""Account for every completed codon case under the existing four-taxon FCS sensitivity gate."""
import argparse,csv,json,hashlib
from collections import defaultdict,Counter
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['fits','exposure','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable output')
    r=json.loads((a.fits/'receipt.json').read_text());e=json.loads((a.exposure/'receipt.json').read_text())
    if e['source_fit_receipt_sha256']!=sha(a.fits/'receipt.json'):raise ValueError('Exposure refers to different codon fits')
    ep=a.exposure/'codon_case_exposure.tsv'
    if sha(ep)!=e['artifacts'][ep.name]:raise ValueError('Changed exposure table')
    flags=defaultdict(set)
    with ep.open() as f:
        for x in csv.DictReader(f,delimiter='\t'):
            if {'EXCLUDE','FIX','TRIM'} & set(x['fcs_actions'].split(';')):flags[x['case_id']].add(x['taxon_id'])
    expected={x['case_id']:x['receipt_sha256'] for x in r['case_receipts']}
    if set(flags)-set(expected):raise ValueError('Unknown exposed case')
    rows=[];proofs=[]
    for case,digest in sorted(expected.items()):
        folder=a.fits/case;cp=folder/'config.json';rp=folder/'receipt.json';receipt=json.loads(rp.read_text());c=json.loads(cp.read_text())
        if sha(rp)!=digest or sha(cp)!=receipt['config_sha256']:raise ValueError('Changed case source')
        command=c['command'];alignment=Path(command[command.index('--alignment')+1])
        if sha(alignment)!=c['alignment_sha256']:raise ValueError('Changed alignment')
        records=list(SeqIO.parse(alignment,'fasta'));taxa={x.id for x in records}
        if len(taxa)!=len(records) or flags[case]-taxa:raise ValueError('Case taxon grid mismatch')
        remain=taxa-flags[case]
        status='unchanged_by_declared_fcs_omission' if not flags[case] else 'below_existing_four_taxon_minimum' if len(remain)<4 else 'requires_coverage_information_screen_and_new_tree'
        rows.append({'case_id':case,'original_taxa':len(taxa),'omitted_taxa':';'.join(sorted(flags[case])),'remaining_taxa':len(remain),'remaining_taxon_ids':';'.join(sorted(remain)),'status':status})
        if flags[case]:proofs.append({'case_id':case,'fit_receipt_sha256':sha(rp),'alignment_path':str(alignment),'alignment_sha256':sha(alignment)})
    a.output.mkdir(parents=True)
    with (a.output/'case_disposition.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result={'status':'complete_full_codon_case_fcs_omission_eligibility_accounting','cases':len(rows),'status_counts':dict(Counter(x['status'] for x in rows)),'source_fit_receipt_sha256':sha(a.fits/'receipt.json'),'source_exposure_receipt_sha256':sha(a.exposure/'receipt.json'),'script_sha256':sha(Path(__file__)),'changed_case_proofs':proofs,'artifacts':{'case_disposition.tsv':sha(a.output/'case_disposition.tsv')},'interpretation':'Apply the existing at-least-four-taxon gate before any codon refit. Below-threshold cases remain unestimable under this project policy, not null effects. No three-taxon policy exception, new alignment, inferred tree or selection test. Unchanged by these flags does not establish eligibility under all other QC criteria.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['changed_case_proofs','artifacts']},indent=2))


if __name__=='__main__':main()
