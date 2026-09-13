#!/usr/bin/env python3
"""Reapply taxon identity policies to source-specific paired marker availability."""
import argparse,json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable output')
    identity=ROOT/'results/phylogeny/taxon-identity-inputs-v1';ir=checked_receipt(identity)
    mp=ROOT/'metadata/analysis_manifest.tsv'
    if sha(mp)!=ir['manifest_sha256']:raise ValueError('Manifest changed')
    taxa={r['taxon_id']:r for r in read_table(mp)};membership=read_table(identity/'taxon_membership.tsv')
    policies={'full_frozen_panel':set(taxa)}
    for r in membership:
        policies.setdefault(r['policy'],set())
        if r['retained']=='True':policies[r['policy']].add(r['taxon_id'])
    marker_rows=[];lineage_rows=[];summaries=[];source_hashes={}
    for source,folder in [('AlphaFold','paired-inputs-gdm-expanded-v1'),('ESMFold','paired-inputs-esmfold-partial-v1')]:
        base=ROOT/'results/phylogeny'/folder;r=checked_receipt(base);source_hashes[source]=sha(base/'receipt.json')
        ready={m['marker'] for m in read_table(base/'marker_summary.tsv') if m['status']=='ready_for_inference'}
        eligible=defaultdict(set)
        for row in read_table(base/'taxon_coverage.tsv'):
            if row['taxon_eligible']=='True':eligible[row['marker']].add(row['taxon_id'])
        for marker in ready:
            actual={x.id for x in SeqIO.parse(base/marker/'aa.faa','fasta')}
            if actual!=eligible[marker]:raise ValueError('Coverage table and alignment identities differ')
        for policy,keep in policies.items():
            usable={m:ts&keep for m,ts in eligible.items() if m in ready and len(ts&keep)>=4}
            observed=set().union(*usable.values()) if usable else set()
            summaries.append({'source':source,'policy':policy,'retained_panel_taxa':len(keep),'ready_markers':len(usable),'taxa_with_usable_marker':len(observed),'usable_taxon_marker_cells':sum(map(len,usable.values())),'lost_ready_markers':';'.join(sorted(ready-set(usable)))})
            for marker in sorted(ready):marker_rows.append({'source':source,'policy':policy,'marker':marker,'original_eligible_taxa':len(eligible[marker]),'retained_eligible_taxa':len(eligible[marker]&keep),'passes_four_taxon_gate':marker in usable})
            for role,group in sorted({(r['study_role'],r['lineage'].split(';')[0]) for r in taxa.values()}):
                members={t for t in keep if taxa[t]['study_role']==role and taxa[t]['lineage'].split(';')[0]==group}
                lineage_rows.append({'source':source,'policy':policy,'role':role,'lineage_group':group,'retained_panel_taxa':len(members),'taxa_with_usable_marker':len(members&observed),'usable_taxon_marker_cells':sum(len(ts&members) for ts in usable.values())})
    a.output.mkdir(parents=True)
    for name,rows in [('summary.tsv',summaries),('marker_coverage.tsv',marker_rows),('lineage_coverage.tsv',lineage_rows)]:write_table(a.output/name,rows)
    result={'status':'complete_identity_policy_paired_coverage','identity_input_receipt_sha256':sha(identity/'receipt.json'),'paired_input_receipt_sha256':source_hashes,'script_sha256':sha(Path(__file__)),'summaries':summaries,'interpretation':'Source-specific availability after taxon exclusions and reapplication of the minimum-four-taxon gate. Original per-taxon confidence/coverage eligibility retained; ready-family FASTA identities verified. No new paired alignments, tree fits, power assessment or species validation implied.','artifacts':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
