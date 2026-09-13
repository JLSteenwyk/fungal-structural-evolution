#!/usr/bin/env python3
"""Map exact-sequence reference availability to project taxa without assuming source identity."""
import argparse,json
from collections import defaultdict,Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['references','snapshot','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable coverage output')
    rr=checked_receipt(a.references);checked_receipt(a.snapshot)
    manifest=ROOT/'metadata/analysis_manifest.tsv';taxa={r['taxon_id']:r for r in read_table(manifest)}
    entries={r['entry_id']:r for r in read_table(a.references/'entries.tsv')};refs=defaultdict(list)
    for r in read_table(a.references/'entities.tsv'):
        if entries[r['entry_id']]['methodology']=='experimental':refs[r['model_id'],r['model_sequence_sha256']].append(r)
    covered=defaultdict(list);links=[];matchedmodels=set()
    for r in read_table(a.snapshot/'marker_structure_links.tsv'):
        hits=refs.get((r['model_id'],r['sequence_sha256']),[])
        if not hits:continue
        matchedmodels.add(r['model_id']);covered[r['taxon_id']].extend(hits)
        links.append({'taxon_id':r['taxon_id'],'marker':r['marker'],'protein_id':r['protein_id'],'model_id':r['model_id'],'sequence_sha256':r['sequence_sha256'],'experimental_entries':len({h['entry_id'] for h in hits}),'experimental_entities':len({h['entity_id'] for h in hits}),'interpretation':'Exact sequence reference candidate; deposited organism may differ from project taxon; no species-specific experimental validation implied.'})
    output=[]
    for taxon,t in taxa.items():
        hits=covered.get(taxon,[]);matching=[r for r in links if r['taxon_id']==taxon]
        output.append({'taxon_id':taxon,'species_name':t['species_name'],'study_role':t['study_role'],'lineage':t['lineage'],'reference_candidate_models':len({r['model_id'] for r in matching}),'reference_candidate_markers':len({r['marker'] for r in matching}),'reference_candidate_entries':len({r['entry_id'] for r in hits}),'reference_candidate_entities':len({r['entity_id'] for r in hits})})
    a.output.mkdir(parents=True);write_table(a.output/'taxon_coverage.tsv',output);write_table(a.output/'marker_reference_links.tsv',links)
    result={'status':'complete_exact_sequence_reference_taxon_coverage','reference_receipt_sha256':sha(a.references/'receipt.json'),'snapshot_receipt_sha256':sha(a.snapshot/'receipt.json'),'manifest_sha256':sha(manifest),'script_sha256':sha(Path(__file__)),'project_taxa':len(output),'covered_taxa':len(covered),'covered_models':len(matchedmodels),'covered_marker_taxon_links':len(links),'covered_roles':dict(Counter(taxa[t]['study_role'] for t in covered)),'covered_major_lineages':dict(Counter(taxa[t]['lineage'].split(';')[0] for t in covered)),'interpretation':'Exact sequence candidate availability projected through the frozen prediction snapshot. Identical sequences can link several project taxa to one deposited reference, whose source organism can differ. This is neither taxon-specific experimental evidence nor independent benchmark sample size. Missing references in this accession-based inventory do not imply no experimental homolog exists.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
