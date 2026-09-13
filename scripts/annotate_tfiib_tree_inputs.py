#!/usr/bin/env python3
"""Join exact family-tree input identities to domain evidence and marker selection."""
import argparse,json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use an immutable annotation output')
    base=ROOT/'results/phylogeny/tfiib-domain-pairs-v1';br=checked_receipt(base)
    source=ROOT/'results/domains/tfiib-candidates-v1';checked_receipt(source)
    if br['candidate_receipt_sha256']!=sha(source/'receipt.json'):raise ValueError('Candidate lineage differs')
    manifest=ROOT/'metadata/analysis_manifest.tsv';taxa={r['taxon_id']:r for r in read_table(manifest)}
    records=list(SeqIO.parse(base/'paired_domains.faa','fasta'));ids={r.id for r in records}
    candidates=read_table(base/'candidate_pair_audit.tsv');eligible={r['entry_id']:r for r in candidates if r['pair_status']=='eligible_ordered_pair'}
    if len(ids)!=len(records) or ids!=set(eligible) or len(ids)!=br['aligned_proteins']:raise ValueError('Input tip identity differs')
    hybrid_config=ROOT/'config/curated_hybrid_evidence.json';hybrids={r['taxon_id'] for r in json.loads(hybrid_config.read_text())['records']}
    rows=[];group=defaultdict(list)
    for tip in sorted(ids):
        c=eligible[tip];t=taxa[c['taxon_id']];selected=json.loads(c['selected_busco_markers_json']);label='BRF1_detected' if int(c['BRF1_hits']) else 'BRF1_not_detected'
        row={**c,'species_name':t['species_name'],'study_role':t['study_role'],'lineage_group':t['lineage'].split(';')[0],'BRF1_evidence_class':label,'selected_focal_marker': '4986044at2759' in selected,'curated_hybrid_taxon':c['taxon_id'] in hybrids,'interpretation':'GA-profile detection class, not validated function or resolved orthology. Non-detection is not proven domain loss; distinct entries need gene/subgenome reconciliation.'}
        rows.append(row);group[c['taxon_id']].append(row)
    taxon_rows=[]
    for taxon,entries in sorted(group.items()):
        counts=Counter(r['BRF1_evidence_class'] for r in entries);selected=[r for r in entries if r['selected_focal_marker']]
        taxon_rows.append({'taxon_id':taxon,'species_name':taxa[taxon]['species_name'],'aligned_candidate_entries':len(entries),'BRF1_detected_entries':counts['BRF1_detected'],'BRF1_not_detected_entries':counts['BRF1_not_detected'],'both_evidence_classes_present':len(counts)==2,'selected_focal_entries':len(selected),'selected_focal_classes':';'.join(sorted({r['BRF1_evidence_class'] for r in selected})),'curated_hybrid_taxon':taxon in hybrids})
    a.output.mkdir(parents=True);write_table(a.output/'tree_tip_annotations.tsv',rows);write_table(a.output/'taxon_class_coverage.tsv',taxon_rows)
    result={'status':'complete_exact_tfiib_input_tip_annotation','input_receipt_sha256':sha(base/'receipt.json'),'candidate_receipt_sha256':sha(source/'receipt.json'),'manifest_sha256':sha(manifest),'hybrid_config_sha256':sha(hybrid_config),'script_sha256':sha(Path(__file__)),'tree_input_tips':len(rows),'taxa':len(taxon_rows),'class_counts':dict(Counter(r['BRF1_evidence_class'] for r in rows)),'taxa_with_both_classes':sum(r['both_evidence_classes_present'] for r in taxon_rows),'selected_focal_class_counts':dict(Counter(r['BRF1_evidence_class'] for r in rows if r['selected_focal_marker'])),'interpretation':'Annotations of the exact 1040 input tips, usable for combined and individual-repeat trees. No fitted topology, clade association, functional identity or ancestral duplication is inferred.','artifacts':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
