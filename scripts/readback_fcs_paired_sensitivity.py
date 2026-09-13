#!/usr/bin/env python3
"""Verify the exact omitted observations and every retained paired character."""
import argparse,csv,hashlib,json
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['baseline','sensitivity','mapping','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new readback receipt')
    r=json.loads((a.sensitivity/'receipt.json').read_text());br=json.loads((a.baseline/'receipt.json').read_text());mr=json.loads((a.mapping/'receipt.json').read_text())
    if r['source_input_receipt_sha256']!=sha(a.baseline/'receipt.json') or r['source_mapping_receipt_sha256']!=sha(a.mapping/'receipt.json'):raise ValueError('Source binding differs')
    for root,receipt in [(a.baseline,br),(a.sensitivity,r),(a.mapping,mr)]:
        for name,digest in receipt['artifacts'].items():
            if sha(root/name)!=digest:raise ValueError('Changed artifact')
    flagged={(x['marker'],x['taxon_id']) for x in rows(a.mapping/'marker_overlap_review.tsv') if set(x['fcs_actions'].split(';')) & {'EXCLUDE','FIX','TRIM'}}
    old={x['marker']:x for x in rows(a.baseline/'marker_summary.tsv')};new={x['marker']:x for x in rows(a.sensitivity/'marker_summary.tsv')}
    if set(old)!=set(new):raise ValueError('Marker universe differs')
    omitted=set();changed=0;characters=0;removed_columns=0
    for marker,summary in old.items():
        if summary['status']!='ready_for_inference':
            if new[marker]['status']!='baseline_not_ready_no_refit':raise ValueError('Changed unready marker')
            continue
        base={label:{x.id:str(x.seq) for x in SeqIO.parse(a.baseline/marker/(label+'.faa'),'fasta')} for label in ['aa','3di']}
        omit={t for m,t in flagged if m==marker and t in base['aa']};remain=set(base['aa'])-omit
        omitted.update((marker,t) for t in omit)
        if not omit:
            if new[marker]['status']!='unchanged_reuse_baseline_fit' or (a.sensitivity/marker).exists():raise ValueError('Unchanged marker unexpectedly refitted')
            continue
        if len(remain)<4:
            if new[marker]['status']!='insufficient_taxa_or_columns_after_omission':raise ValueError('Underpopulated marker not flagged')
            continue
        changed+=1
        col=rows(a.sensitivity/marker/'columns.tsv');indices=[int(x['baseline_paired_column_1based'])-1 for x in col]
        expected=[i for i in range(len(next(iter(base['aa'].values())))) if any(base['aa'][t][i]!='?' for t in remain)]
        if indices!=expected:raise ValueError('Retained column grid differs')
        removed_columns+=len(next(iter(base['aa'].values())))-len(indices)
        for label in ['aa','3di']:
            fresh={x.id:str(x.seq) for x in SeqIO.parse(a.sensitivity/marker/(label+'.faa'),'fasta')}
            if set(fresh)!=remain:raise ValueError('Retained taxon grid differs')
            for taxon,sequence in fresh.items():
                if sequence!=''.join(base[label][taxon][i] for i in indices):raise ValueError('Retained character changed')
                characters+=len(sequence)
    actual={(x['marker'],x['taxon_id']) for x in rows(a.sensitivity/'omitted_marker_observations.tsv')}
    if actual!=omitted or len(actual)!=r['omitted_marker_observations'] or changed!=r['ready_markers']:raise ValueError('Omission/refit universe differs')
    result={'status':'passed_exact_fcs_omission_and_retained_character_readback','source_receipt_sha256':sha(a.sensitivity/'receipt.json'),'script_sha256':sha(Path(__file__)),'baseline_markers':len(old),'changed_ready_markers':changed,'omitted_observations':len(omitted),'retained_characters_checked_both_alphabets':characters,'now_all_missing_columns_removed':removed_columns,'interpretation':'Every emitted taxon/column/character matches the exact baseline projection after the declared FCS action filter. Unchanged and previously unready marker dispositions verified. No inference or contamination confirmation.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
