#!/usr/bin/env python3
"""Retain all case screens and compare reference-dependent signs across fit variants."""
import csv,json
from collections import defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results');CASE=['family','gene_a','gene_b','pfam_accession']
ID=['triad_id','domain_triad','entity_id','experimental_model','label_asym_id','mask','context_indices_json']
UNIT=['triad_id','entity_id','experimental_model','label_asym_id','context_indices_json']
METRICS=['whole','domain','outside_independent','outside_domain_anchored']
OUTPUT=BASE/'experimental_structures/whole-domain-case-reference-robustness-20260927-v1'


def direction(values,margin):
    if not values:return 'unavailable'
    if min(values)>margin+1e-8:return 'positive'
    if max(values)<-margin-1e-8:return 'negative'
    return 'variable_or_within_margin'


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name;assert sha(p)==r['artifacts'][name];sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return read(p)
    experimental=BASE/'experimental_structures';root=experimental/'whole-domain-case-quartet-contrasts-20260927-v1'
    ap=experimental/'whole-domain-case-quartet-contrast-readback-20260927-v1/receipt.json';audit=json.loads(ap.read_text());assert audit['status']=='complete_full_case_quartet_contrast_readback' and audit['source_hashes'][str(root/'receipt.json')]==sha(root/'receipt.json');sources[str(ap)]=sha(ap)
    cases=BASE/'structural_comparisons/whole-domain-case-dossiers-20260927-v1'
    dossiers=checked(cases,'case_dossiers.tsv');links=checked(cases,'domain_reference_links.tsv');whole=checked(cases,'whole_reference_links.tsv')
    fields=['guide','family','gene_a','gene_b','reference_gene'];wl={tuple(r[k] for k in fields):r['triad_id'] for r in whole};case_triads=defaultdict(set);variants=defaultdict(set)
    for r in links:
        ck=tuple(r[k] for k in CASE);tid=wl[tuple(r[k] for k in fields)];case_triads[ck].add(tid)
        for mask in ['full','plddt70']:variants[ck,tid].add((r['triad_key'],mask))
    triad_cases=defaultdict(set)
    for ck,tids in case_triads.items():
        for tid in tids:triad_cases[tid].add(ck)
    cover={}
    for r in checked(root,'coverage_screens.tsv'):
        key=tuple(r[k] for k in ID)+(r['screen'],);assert key not in cover;cover[key]=r
    grouped=defaultdict(list)
    for r in checked(root,'paired_reference_contrasts.tsv'):
        for ck in triad_cases[r['triad_id']]:grouped[ck,tuple(r[k] for k in UNIT),r['metric']].append(r)
    dependence={r['entity_id']:r for r in checked(experimental/'whole-domain-case-experimental-dependence-20260927-v1','entity_components.tsv')}
    metadata=experimental/'whole-domain-case-metadata-review-20260927-v1'
    entries={r['entry_id']:r for r in checked(metadata,'entries.tsv')};entities={r['entity_id']:r for r in checked(metadata,'entities.tsv')}
    OUTPUT.mkdir(exist_ok=False);path=OUTPUT/'reference_unit_sensitivity.tsv';aggregate=defaultdict(list);nrows=0
    with path.open('w') as f:
        writer=None
        for (ck,unit,metric),records in sorted(grouped.items()):
            assert {(r['domain_triad'],r['mask']) for r in records}==variants[ck,unit[0]] and len(records)==len(variants[ck,unit[0]])
            requirements={r['coverage_requirement'] for r in records};assert len(requirements)==1;requirement=next(iter(requirements))
            numeric=all(r['fungal_reference_status']==r['experimental_reference_status']=='computed' for r in records)
            fvalues=[float(r['fungal_reference_ar_minus_br']) for r in records if r['fungal_reference_status']=='computed'];evalues=[float(r['experimental_reference_ae_minus_be']) for r in records if r['experimental_reference_status']=='computed']
            entity=unit[1];entry=entity.rsplit('_',1)[0];dep=dependence[entity]
            for n in [30,50]:
                for c in [50,70,90]:
                    screen=f'n{n}_c{c}'
                    coverage=all(cover[tuple(r[k] for k in ID)+(screen,)][requirement+'_pass']=='1' for r in records)
                    for margin in [0.,.01,.1]:
                        qualified=numeric and coverage
                        fd=direction(fvalues,margin) if qualified else 'not_qualified';ed=direction(evalues,margin) if qualified else 'not_qualified'
                        state='not_qualified'
                        if qualified:state='same_direction' if fd==ed and fd in ['positive','negative'] else 'opposite_direction' if {fd,ed}=={'positive','negative'} else 'variable_or_within_margin'
                        row=dict(zip(CASE,ck),**dict(zip(UNIT,unit)),metric=metric,screen=screen,margin_angstrom=margin,required_variants=len(records),all_coverage_pass=int(coverage),all_numeric_pass=int(numeric),qualified=int(qualified),fungal_min=min(fvalues) if fvalues else '',fungal_max=max(fvalues) if fvalues else '',experimental_min=min(evalues) if evalues else '',experimental_max=max(evalues) if evalues else '',fungal_direction=fd,experimental_direction=ed,reference_agreement=state,sequence_sha256=dep['sequence_sha256'],dependency_component=dep['component_id'],reported_mutation_count=entities[entity]['mutation_count'],starting_model_category=entries[entry]['starting_model_category'])
                        if writer is None:writer=csv.DictWriter(f,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                        writer.writerow(row);nrows+=1;aggregate[ck,metric,screen,str(margin)].append(row)
    summaries=[]
    for case in dossiers:
        ck=tuple(case[k] for k in CASE)
        for metric in METRICS:
            for n in [30,50]:
                for c in [50,70,90]:
                    screen=f'n{n}_c{c}'
                    for margin in [0.,.01,.1]:
                        data=aggregate[ck,metric,screen,str(margin)];eligible=[r for r in data if r['qualified']]
                        row=dict(zip(CASE,ck),species_name=case['species_name'],metric=metric,screen=screen,margin_angstrom=margin,candidate_units=len(data),qualified_units=len(eligible),qualified_entities=len({r['entity_id'] for r in eligible}),qualified_entries=len({r['entity_id'].rsplit('_',1)[0] for r in eligible}),qualified_sequences=len({r['sequence_sha256'] for r in eligible}),qualified_dependency_components=len({r['dependency_component'] for r in eligible}))
                        for state in ['same_direction','opposite_direction','variable_or_within_margin']:
                            selected=[r for r in eligible if r['reference_agreement']==state];row[state+'_units']=len(selected);row[state+'_entities']=len({r['entity_id'] for r in selected})
                        assert sum(row[s+'_units'] for s in ['same_direction','opposite_direction','variable_or_within_margin'])==row['qualified_units'];summaries.append(row)
    assert len(summaries)==13*4*6*3
    p=OUTPUT/'case_summary.tsv'
    with p.open('w') as f:w=csv.DictWriter(f,list(summaries[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summaries)
    assert read(p)==[{k:str(v) for k,v in r.items()} for r in summaries]
    assert len(read(path))==nrows
    for p,digest in sources.items():assert sha(p)==digest
    result=dict(status='complete_descriptive_case_reference_robustness',source_hashes=sources,script_sha256=sha(__file__),unit_sensitivity_rows=nrows,case_summary_rows=len(summaries),artifacts={p.name:sha(p) for p in OUTPUT.iterdir()},scope='All case/metric/screen/margin combinations, including zero cases. Each unit fixes fungal triplet, experimental entity/model/chain and alignment contexts; qualification requires every domain-boundary/full-plddt70 variant. Direction requires all signed contrasts beyond margin plus 1e-8 numeric tolerance. Units and entity-level categories are not independent; entity categories may overlap across chains. Descriptive, selected-case reference sensitivity, not mechanism, ancestral polarity, experimental validation or significance.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in summaries:
        if r['screen']=='n50_c70' and r['margin_angstrom']==.1 and r['qualified_units']:print(r['family'],r['metric'],r['qualified_units'],r['same_direction_units'],r['opposite_direction_units'],r['variable_or_within_margin_units'])


if __name__=='__main__':main()
