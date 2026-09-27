#!/usr/bin/env python3
"""Preserve species-level source classifications without inferring binary absence."""
import csv,json,hashlib
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs):
    with open(p,'w') as f:
        w=csv.DictWriter(f,list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def main():
    config=Path('config/ecology_species_source_review_20260927.json');c=json.loads(config.read_text())
    for p,h in c['pins'].items():assert sha(p)==h,p
    old=rows(c['previous']);sources={r['species_name']:r for r in rows(c['source'])};identity={r['species_name']:r for r in rows(c['identity'])}
    proof=json.loads(Path(c['source_readback']).read_text());assert proof['table_sha256']==sha(c['source']) and proof['status']=='passed_raw_xlsx_xml_readback'
    assert json.loads(Path(c['previous_receipt']).read_text())['table_sha256']==sha(c['previous'])
    assert set(c['decisions'])==set(sources)-{r['species_name'] for r in old}
    out=list(old);review=[]
    for name,d in c['decisions'].items():
        s=sources[name];i=identity[name];assert s['source_ecology']==d['source_label']
        assert i['selected_assembly_accession']==s['selected_assembly_accession']
        review.append(dict(taxon_id=s['taxon_id'],species_name=name,source_label=d['source_label'],disposition=d['disposition'],retained_state=d.get('state',''),note=d['note'],sample_link_status=i['status']))
        if d['disposition']!='retain_published_species_classification':continue
        row=dict.fromkeys(old[0],'')
        row.update(taxon_id=s['taxon_id'],assembly_accession=s['selected_assembly_accession'],trait='ecological_guild_as_classified_by_source',species_name=name,state=d['state'],source_species_name=s['source_species_name'],source_doi=s['source_doi'],source_url=s['source_url'],source_locator=f"Supplementary Dataset 2, {s['source_sheet']}, Excel row {s['source_excel_row']}",evidence='Species-level classification explicitly reported in the supplementary table.',provisional_transition_group='unassigned_'+s['taxon_id'],scope_note=d['note']+' Published species classification; selected isolate equivalence unresolved ('+i['status']+'). Categories are not exclusive. No binary ECM absence inferred; no independent transition assigned.',evidence_level='published_species_classification',genus_candidate_primary_lifestyle=s['genus_candidate_primary_lifestyle'],conflicts_with_genus_ECM_candidate=s['genus_ECM_disagreement'],selected_isolate_experimentally_verified='False',confirmatory_test_status='pending_taxonomy_phylogeny_and_replicated_transition_review')
        out.append(row)
    assert len(out)==45 and len({r['taxon_id'] for r in out})==45 and out[:32]==old
    output=Path(c['output']);dispositions=Path(c['dispositions']);assert not output.exists() and not dispositions.exists()
    write(output,out);write(dispositions,review)
    r=dict(status='complete_reviewed_species_classification_extension',statements=len(out),additions=len(out)-len(old),reviewed_pending_records=len(review),config_sha256=sha(config),script_sha256=sha(__file__),table_sha256=sha(output),artifacts={str(output):sha(output),str(dispositions):sha(dispositions)},scope='Original 32 statements preserved; 13 published species-level ecological labels added with isolate and transition uncertainty explicit. Yeast morphology kept separate; Punctularia pathogen coding withheld. No ECM-negative coding, selected-isolate phenotype validation or independent ecological effect claim.')
    Path(c['receipt']).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
