#!/usr/bin/env python3
"""Retain method, construct and starting-model annotations for every case candidate."""
import csv
import json
from collections import Counter
from pathlib import Path

from screen_duplication_domain_alignment_coverage import sha
from review_experimental_model_provenance import classify

BASE = Path('results/experimental_structures')
ROOT = BASE / 'whole-domain-case-metadata-20260927-v1'
AUDIT = BASE / 'whole-domain-case-subject-readback-20260927-v1/receipt.json'
OUTPUT = BASE / 'whole-domain-case-metadata-review-20260927-v1'


def packed(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def main():
    audit = json.loads(AUDIT.read_text())
    assert audit['status'] == 'complete_case_subject_metadata_readback'
    rp = ROOT / 'receipt.json'
    assert audit['source_hashes'][str(rp)] == sha(rp)
    receipt = json.loads(rp.read_text())
    sources = {str(rp): sha(rp), str(AUDIT): sha(AUDIT),
               'scripts/review_experimental_model_provenance.py': sha('scripts/review_experimental_model_provenance.py')}
    entries, entities, annotations = [], [], []
    raw = {}
    for item in receipt['responses']:
        path = ROOT / item['kind'] / (item['identifier'] + '.json')
        assert sha(path) == item['response_sha256']
        sources[str(path)] = sha(path)
        data = json.loads(path.read_text())
        raw[item['kind'], item['identifier']] = data
        if item['kind'] == 'entry':
            info = data['rcsb_entry_info']
            starting = data.get('pdbx_initial_refinement_model', [])
            entry = dict(entry_id=item['identifier'],
                         methodology=info.get('structure_determination_methodology', 'unknown'),
                         methods_json=packed([r.get('method') for r in data.get('exptl', [])]),
                         resolutions_json=packed(info.get('resolution_combined', [])),
                         initial_release_date=data.get('rcsb_accession_info', {}).get('initial_release_date', 'unknown'),
                         starting_model_category=classify(starting),
                         alphafold_mentioned_in_starting_annotations=int('alphafold' in packed(starting).lower()),
                         alphafold_mentioned_in_software_annotations=int('alphafold' in packed(data.get('software', [])).lower()),
                         target_chain_attribution='unresolved_entry_level_annotation',
                         training_independence='unresolved',
                         primary_citation_pubmed_id=data.get('rcsb_primary_citation', {}).get('pdbx_database_id_PubMed', ''),
                         citation_doi=data.get('rcsb_primary_citation', {}).get('pdbx_database_id_DOI', ''))
            entries.append(entry)
            keys = [k for k in data if k in ('refine','em_3d_reconstruction','pdbx_initial_refinement_model','software','rcsb_primary_citation','exptl') or k.startswith('pdbx_vrpt_summary')]
            annotations.append(dict(entry_id=item['identifier'],fields={k: data[k] for k in sorted(keys)}))
        elif item['kind'] == 'polymer_entity':
            polymer = data['entity_poly']
            summary = data.get('rcsb_polymer_entity', {})
            entities.append(dict(entity_id=item['identifier'],entry_id=item['identifier'].rsplit('_',1)[0],
                                 canonical_length=len(''.join(polymer['pdbx_seq_one_letter_code_can'].split())),
                                 mutation_count=polymer.get('rcsb_mutation_count','unknown'),
                                 nonstandard_monomer_count=polymer.get('rcsb_non_std_monomer_count','unknown'),
                                 mutation_annotation=summary.get('pdbx_mutation','unknown'),
                                 fragment_annotation=summary.get('pdbx_fragment','unknown'),
                                 description=summary.get('pdbx_description','unknown'),
                                 source_organisms_json=packed(data.get('rcsb_entity_source_organism',[])),
                                 host_organisms_json=packed(data.get('rcsb_entity_host_organism',[])),
                                 source_construct_fields_json=packed({k:v for k,v in data.items() if k in ('entity_src_gen','entity_src_nat','pdbx_entity_src_syn')})))
    assert len(entries)==receipt['entry_count']==707 and len(entities)==receipt['entity_count']==1205
    assert {r['entry_id'] for r in entities}=={r['entry_id'] for r in entries}
    OUTPUT.mkdir(exist_ok=False)
    artifacts={}
    for name, rows in [('entries.tsv', entries), ('entities.tsv', entities)]:
        path=OUTPUT/name
        with path.open('w') as f:
            writer=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
        with path.open() as f:
            readback=list(csv.DictReader(f,delimiter='\t'))
        assert readback==[{k:str(v) for k,v in r.items()} for r in rows]
        artifacts[name]=sha(path)
    path=OUTPUT/'method_specific_annotations.json'
    path.write_text(json.dumps(annotations,indent=2)+'\n')
    assert json.loads(path.read_text())==annotations
    artifacts[path.name]=sha(path)
    # Re-read source objects to check every exported method/construct field directly.
    for row in entries:
        data=json.loads((ROOT/'entry'/(row['entry_id']+'.json')).read_text())
        assert json.loads(row['methods_json'])==[r.get('method') for r in data.get('exptl',[])]
        assert json.loads(row['resolutions_json'])==data['rcsb_entry_info'].get('resolution_combined',[])
    for row in entities:
        data=json.loads((ROOT/'polymer_entity'/(row['entity_id']+'.json')).read_text())
        assert json.loads(row['source_organisms_json'])==data.get('rcsb_entity_source_organism',[])
        assert json.loads(row['source_construct_fields_json'])=={k:v for k,v in data.items() if k in ('entity_src_gen','entity_src_nat','pdbx_entity_src_syn')}
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='complete_case_experimental_metadata_annotation_review',source_hashes=sources,
                script_sha256=sha(__file__),entries=len(entries),entities=len(entities),
                methodologies=dict(Counter(r['methodology'] for r in entries)),
                starting_model_categories=dict(Counter(r['starting_model_category'] for r in entries)),
                entries_with_alphafold_starting_annotation=sum(r['alphafold_mentioned_in_starting_annotations'] for r in entries),
                entries_with_alphafold_software_annotation=sum(r['alphafold_mentioned_in_software_annotations'] for r in entries),
                entities_with_positive_reported_mutation_count=sum(isinstance(r['mutation_count'],int) and r['mutation_count']>0 for r in entities),
                entities_with_unknown_mutation_count=sum(r['mutation_count']=='unknown' for r in entities),
                artifacts=artifacts,
                scope='All candidates retained, including integrative entry. Metadata annotations only; missing annotation is not absence of modification or prediction-assisted refinement. Resolution and validation metrics remain method-specific. Entry starting models are not attributed to the target chain. Release dates are not training/template cutoffs; independence unresolved for every entry. No quality-based ranking or geometry acceptance.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','artifacts')},indent=2))


if __name__=='__main__':main()
