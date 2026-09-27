#!/usr/bin/env python3
"""Link every case candidate to metadata without excluding unfavorable annotations."""
import csv
import json
from collections import defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results/experimental_structures')
KEY=['family','gene_a','gene_b','pfam_accession','screen']


def read(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def main():
    sources={}
    def checked(directory,name):
        root=BASE/directory;rp=root/'receipt.json';receipt=json.loads(rp.read_text());path=root/name
        assert sha(path)==receipt['artifacts'][name]
        sources[str(path)]=sha(path);sources[str(rp)]=sha(rp)
        return read(path)
    coverage='whole-domain-case-domain-coverage-20260927-v1'
    review='whole-domain-case-metadata-review-20260927-v1'
    candidates=checked(coverage,'shared_entity_screens.tsv')
    cases=checked(coverage,'case_summary.tsv')
    entry_rows=checked(review,'entries.tsv');entity_rows=checked(review,'entities.tsv')
    entries={r['entry_id']:r for r in entry_rows};entities={r['entity_id']:r for r in entity_rows}
    assert len(entries)==len(entry_rows) and len(entities)==len(entity_rows)
    linked=[];groups=defaultdict(list);seen=set()
    for row in candidates:
        key=tuple(row[k] for k in KEY);entity=entities[row['entity_id']];entry=entries[entity['entry_id']]
        assert (*key,row['entity_id']) not in seen;seen.add((*key,row['entity_id']))
        joined=dict(row,**{k:v for k,v in entry.items() if k!='entry_id'},entry_id=entry['entry_id'],
                reported_mutation_count=entity['mutation_count'],reported_nonstandard_monomer_count=entity['nonstandard_monomer_count'],
                mutation_annotation=entity['mutation_annotation'],fragment_annotation=entity['fragment_annotation'],
                source_organisms_json=entity['source_organisms_json'])
        linked.append(joined);groups[key].append(joined)
    summary=[]
    for case in cases:
        key=tuple(case[k] for k in KEY);group=groups[key]
        for region,field,total_field in [('domain','shared_domain_pass','shared_domain_entities'),('domain_and_outside','shared_domain_and_outside_pass','shared_domain_and_outside_entities')]:
            chosen=[r for r in group if r[field]=='1'];assert len(chosen)==int(case[total_field])
            experimental=[r for r in chosen if r['methodology']=='experimental']
            flagged=[r for r in experimental if r['alphafold_mentioned_in_starting_annotations']=='1']
            mutated=[r for r in experimental if r['reported_mutation_count']!='unknown' and int(r['reported_mutation_count'])>0]
            missing=[r for r in experimental if r['starting_model_category']=='starting_model_annotation_missing']
            summary.append({**{k:case[k] for k in KEY},'species_name':case['species_name'],'coverage_region':region,
                            'all_candidate_entities':len(chosen),'experimental_entities':len(experimental),
                            'nonexperimental_entities':len(chosen)-len(experimental),
                            'experimental_entries':len({r['entry_id'] for r in experimental}),
                            'alphafold_starting_annotation_entities':len(flagged),
                            'alphafold_starting_annotation_entries':len({r['entry_id'] for r in flagged}),
                            'alphafold_starting_entry_ids_json':json.dumps(sorted({r['entry_id'] for r in flagged})),
                            'mutated_experimental_entities':len(mutated),
                            'missing_starting_annotation_entries':len({r['entry_id'] for r in missing})})
    assert len(summary)==156
    out=BASE/'whole-domain-case-metadata-links-20260927-v1';out.mkdir(exist_ok=False)
    artifacts={}
    for name,rows in [('annotated_candidate_entities.tsv',linked),('case_metadata_summary.tsv',summary)]:
        path=out/name
        with path.open('w') as handle:
            writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
        assert read(path)==[{k:str(v) for k,v in r.items()} for r in rows];artifacts[name]=sha(path)
    for path,digest in sources.items():assert sha(path)==digest
    receipt=dict(status='complete_case_experimental_metadata_links',source_hashes=sources,script_sha256=sha(__file__),
                 candidate_rows=len(linked),case_screen_region_rows=len(summary),artifacts=artifacts,
                 scope='All candidates joined by exact entity and entry IDs. Counts stratify alignment coverage only, not resolved-coordinate coverage or quality. Entry-level starting-model flags are not target-chain assignments. Multiple entities, entries and cases may share experiments; no independence or accuracy claim. Zero-candidate case/screens retained. Mutation flags are reported metadata, not inferred effects.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    for row in summary:
        if row['screen']=='n50_c70' and row['coverage_region']=='domain_and_outside':print(row['family'],row['experimental_entities'],row['experimental_entries'],row['alphafold_starting_annotation_entries'],row['mutated_experimental_entities'])


if __name__=='__main__':main()
