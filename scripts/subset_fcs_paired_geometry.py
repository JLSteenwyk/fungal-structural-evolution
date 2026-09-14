#!/usr/bin/env python3
"""Reuse exact remaining-pair geometry when an FCS omission leaves columns unchanged."""
import argparse,csv,json,hashlib
from collections import Counter
from itertools import combinations
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline-inputs','geometry','geometry-audit','sensitivity','readback','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable output')
    gr=json.loads((a.geometry/'receipt.json').read_text());ga=json.loads((a.geometry_audit/'receipt.json').read_text());sr=json.loads((a.sensitivity/'receipt.json').read_text());rr=json.loads(a.readback.read_text())
    if ga['source_receipts']['comparisons']!=sha(a.geometry/'receipt.json') or ga['status']!='passed_complete_paired_grid_character_and_sampled_geometry_readback':raise ValueError('Missing baseline geometry audit')
    if gr['source_receipts']['inputs']['sha256']!=sha(a.baseline_inputs/'receipt.json') or sr['source_input_receipt_sha256']!=sha(a.baseline_inputs/'receipt.json'):raise ValueError('Baseline lineage mismatch')
    if rr['source_receipt_sha256']!=sha(a.sensitivity/'receipt.json') or rr['status']!='passed_exact_fcs_omission_and_retained_character_readback' or rr['now_all_missing_columns_removed']!=0:raise ValueError('Recompute geometry if columns changed or readback missing')
    for root,receipt in [(a.geometry,gr),(a.sensitivity,sr)]:
        for name,digest in receipt['artifacts'].items():
            if sha(root/name)!=digest:raise ValueError('Changed artifact')
    ready={r['marker']:r for r in rows(a.sensitivity/'marker_summary.tsv') if r['status']=='ready_for_inference'};taxa={};expected=set()
    for marker in ready:
        taxa[marker]={r.id for r in SeqIO.parse(a.sensitivity/marker/'aa.faa','fasta')}
        expected.update((marker,x,y) for x,y in combinations(sorted(taxa[marker]),2))
    a.output.mkdir(parents=True);seen=set();counts=Counter();source_counts=Counter();field_checks=0
    for name,status in [('paired_site_geometry.tsv','accepted'),('coverage_exclusions.tsv','excluded')]:
        with (a.geometry/name).open() as src,(a.output/name).open('w',newline='') as dst:
            reader=csv.DictReader(src,delimiter='\t');writer=csv.DictWriter(dst,reader.fieldnames,delimiter='\t',lineterminator='\n');writer.writeheader()
            for row in reader:
                source_counts[status]+=1;key=row['marker'],row['taxon_a'],row['taxon_b']
                if key not in expected:continue
                if key in seen:raise ValueError('Duplicate retained pair')
                seen.add(key);new=dict(row);new['tree_taxa']=str(len(taxa[row['marker']]))
                if int(new['tree_alignment_columns'])!=int(ready[row['marker']]['retained_columns']):raise ValueError('Alignment columns changed')
                writer.writerow(new);counts[status]+=1
        # Full output/source readback: only the cohort taxon-count annotation may change.
        output=iter(rows(a.output/name));emitted=0
        with (a.geometry/name).open() as src:
            for row in csv.DictReader(src,delimiter='\t'):
                if (row['marker'],row['taxon_a'],row['taxon_b']) not in expected:continue
                actual=next(output);row['tree_taxa']=str(len(taxa[row['marker']]))
                if actual!=row:raise ValueError('Subset changes source measurement')
                field_checks+=len(row);emitted+=1
        if next(output,None) is not None or emitted!=counts[status]:raise ValueError('Wrong output row count')
    if seen!=expected:raise ValueError('Incomplete remaining-taxon pair grid')
    if source_counts['accepted']!=gr['accepted_taxon_pairs'] or source_counts['excluded']!=gr['excluded_taxon_pairs']:raise ValueError('Baseline count mismatch')
    source_receipts=dict(gr['source_receipts']);source_receipts['inputs']={'path':str(a.sensitivity),'sha256':sha(a.sensitivity/'receipt.json')}
    result={'status':'complete_exact_retained_pair_geometry_subset','markers':len(ready),'accepted_taxon_pairs':counts['accepted'],'excluded_taxon_pairs':counts['excluded'],'source_receipts':source_receipts,'baseline_geometry_receipt_sha256':sha(a.geometry/'receipt.json'),'baseline_geometry_audit_sha256':sha(a.geometry_audit/'receipt.json'),'sensitivity_readback_sha256':sha(a.readback),'script_sha256':sha(Path(__file__)),'all_output_source_fields_checked':field_checks,'interpretation':'Exact baseline measurements retained for all remaining pairs of changed sensitivity markers. No coordinates, masks or alignment columns changed; only tree_taxa cohort annotations updated. Complete pair grid and every output field read back against baseline. Fresh sampled geometry audit and new fitted-tree path benchmark remain required. Unchanged-marker geometry remains in the baseline. No new geometry inference or branch result.','artifacts':{n:sha(a.output/n) for n in ['paired_site_geometry.tsv','coverage_exclusions.tsv']}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['artifacts','source_receipts']},indent=2))


if __name__=='__main__':main()
