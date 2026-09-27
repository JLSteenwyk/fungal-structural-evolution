#!/usr/bin/env python3
"""Inventory frozen model coverage for all audited reported duplication events."""
import argparse,csv,gzip,json,sqlite3
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

HEADER=['Orthogroup','Species Tree Node','Gene Tree Node','Support','Type','Genes 1','Genes 2']


def event_coverage(row, models):
    if len(row)!=7:raise ValueError('Malformed event')
    sides=[s.split(', ') if s else [] for s in row[5:]]
    if not sides[0] or len(set(sides[0]+sides[1]))!=sum(map(len,sides)):raise ValueError('Empty or repeated event genes')
    covered=[[g for g in side if g in models] for side in sides]
    simple=row[4]=='Terminal' and all(len(s)==1 for s in sides)
    both=all(covered)
    return sides,covered,simple,both


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text())
    def verify():
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin '+path)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    audit=json.loads(Path(plan['event_audit']).read_text());ap=json.loads(Path(plan['event_plan']).read_text())
    if audit['status']!='passed_full_duplication_event_identity_audit' or audit['plan_sha256']!=sha(plan['event_plan']):raise ValueError('Event audit mismatch')
    bridge=Path(plan['bridge']);br=json.loads(Path(plan['bridge_receipt']).read_text())
    readback=json.loads(Path(plan['bridge_readback']).read_text())
    if readback['status']!='passed_independent_full_structure_family_coverage_readback' or readback['producer_receipt_sha256']!=sha(plan['bridge_receipt']):raise ValueError('Bridge audit mismatch')
    if sha(bridge)!=br['artifacts'][bridge.name]:raise ValueError('Bridge database mismatch')
    db=sqlite3.connect('file:'+str(bridge.resolve())+'?mode=ro',uri=True)
    models={}
    for taxon,protein,seq,model in db.execute('SELECT taxon_id,protein_id,sequence_sha256,model_id FROM structures'):
        key=taxon+'_'+protein
        if key in models:raise ValueError('Ambiguous model label')
        models[key]=(seq,model)
    if len(models)!=br['protein_links']:raise ValueError('Bridge protein count mismatch')
    db.close();out.mkdir(parents=True);summaries=[]
    for entry in ap['guides']:
        guide=entry['guide'];source=Path(entry['result'])/'Gene_Duplication_Events/Duplications.tsv'
        if sha(source)!=ap['pins'][str(source)]:raise ValueError('Changed audited event source')
        counts=Counter();families=set();taxa=set()
        with source.open() as handle,gzip.open(out/(guide+'_event_coverage.tsv.gz'),'wt') as full,(out/(guide+'_terminal_two_model_candidates.tsv')).open('w') as candidates:
            reader=csv.reader(handle,delimiter='\t');assert next(reader)==HEADER
            w=csv.writer(full,delimiter='\t',lineterminator='\n');cw=csv.writer(candidates,delimiter='\t',lineterminator='\n')
            w.writerow(['family','species_node','gene_node','support','type','left_genes','right_genes','left_models','right_models','terminal_singleton_sides','both_sides_covered'])
            cw.writerow(['family','taxon_id','gene_node','support','gene_a','gene_b','model_a','model_b','same_sequence','same_model'])
            for row in reader:
                sides,covered,simple,both=event_coverage(row,models);counts['events']+=1;counts['genes_across_events']+=sum(map(len,sides));counts['covered_genes_across_events']+=sum(map(len,covered));counts['both_sides_covered']+=both;counts['empty_second_side']+=not sides[1];counts['terminal_singleton_sides']+=simple
                counts['events_with_any_model']+=any(covered)
                w.writerow(row[:5]+[len(sides[0]),len(sides[1]),len(covered[0]),len(covered[1]),int(simple),int(both)])
                if simple and both:
                    ga,gb=sorted(sides[0]+sides[1]);sa,ma=models[ga];sb,mb=models[gb]
                    if not ga.startswith(row[1]+'_') or not gb.startswith(row[1]+'_'):raise ValueError('Terminal taxon mismatch')
                    cw.writerow([row[0],row[1],row[2],row[3],ga,gb,ma,mb,int(sa==sb),int(ma==mb)])
                    counts['terminal_two_model_candidates']+=1;counts['candidate_same_sequence']+=sa==sb;counts['candidate_same_model']+=ma==mb;families.add(row[0]);taxa.add(row[1])
        expected=next(r for r in audit['guides'] if r['guide']==guide)
        if counts['events']!=expected['counts']['events'] or counts['genes_across_events']!=expected['counts']['gene_assignments_across_events']:raise ValueError('Audited event universe differs')
        summaries.append(dict(guide=guide,counts=dict(counts),candidate_families=len(families),candidate_taxa=len(taxa)));print(json.dumps(summaries[-1]),flush=True)
    verify();r=dict(status='complete_duplication_frozen_structure_coverage',plan_sha256=sha(a.plan),script_sha256=sha(__file__),model_protein_links=len(models),guides=summaries,artifacts={p.name:sha(p) for p in out.iterdir()},scope='All reported events in both audited guides; frozen AlphaFold bridge coverage only. Nested gene/event counts are not independent observations. Terminal singleton-side candidates require gene-tree membership validation, model quality and alignment checks, and phylogenetic controls before duplication effects or asymmetry tests. Missing models are not biological absences; identical sequences/models remain explicitly flagged.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')

if __name__=='__main__':
    csv.field_size_limit(32*1024*1024);main()
