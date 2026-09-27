#!/usr/bin/env python3
"""Check all taxonomy annotations and sampling summaries using independent CSV sets/counters."""
import argparse,csv,hashlib,json
from collections import defaultdict,Counter
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def fresh():return dict(n=0,pairs=set(),families=set(),pfams=set(),taxa=set())
def add(b,r):
    b['n']+=1;b['pairs'].add((r['family'],r['gene_a'],r['gene_b']));b['families'].add(r['family']);b['pfams'].add(r['pfam_accession']);b['taxa'].add(r['taxon'])

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();p=json.loads(a.plan.read_text());root=Path(p['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_domain_candidate_sampling_pending_readback' and r['plan_sha256']==sha(a.plan)
    for path,h in p['pins'].items():assert sha(path)==h,path
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    taxa={};lineage_taxa=defaultdict(set)
    for x in csv.DictReader((Path(p['membership'])/'taxon_coverage_with_membership.tsv').open(),delimiter='\t'):
        if x['in_reconciliation']!='1':continue
        values={f:x[f] for f in ['species_name','study_role','lineage']};values['lineage_group']=x['lineage'].split(';')[0]
        assert x['taxon'] not in taxa or taxa[x['taxon']]==values
        taxa[x['taxon']]=values;lineage_taxa[values['study_role'],values['lineage_group']].add(x['taxon'])
    keys=['family','gene_a','gene_b','pfam_accession','screen','margin'];original={}
    for x in csv.DictReader((Path(p['comparison'])/'guide_comparison.tsv').open(),delimiter='\t'):
        k=tuple(x[f] for f in keys);assert k not in original;original[k]=x
    lineage=defaultdict(fresh);family=defaultdict(fresh);overview=defaultdict(fresh);roles=defaultdict(Counter);lineage_overview=defaultdict(set);seen=set();screens=set()
    categories=['structural_direction_unresolved','guide_specific_domain','stable_structure_sequence_unresolved','stable_concordant','stable_discordant','stable_structure_sequence_guide_sensitive']
    for x in csv.DictReader((root/'annotated_comparisons.tsv').open(),delimiter='\t'):
        k=tuple(x[f] for f in keys);assert k not in seen;seen.add(k);src=original[k]
        for f,v in src.items():assert x[f]==v,f
        taxon=src['mafft_taxon_id'] or src['profile_taxon_id'];assert x['taxon']==taxon
        for f,v in taxa[taxon].items():assert x[f]==v
        if src['structural_guide_agreement']=='same_stable_direction':
            category={'concordant_in_both_guides':'stable_concordant','discordant_in_both_guides':'stable_discordant','guide_sensitive_relationship':'stable_structure_sequence_guide_sensitive','unresolved_or_missing':'stable_structure_sequence_unresolved'}[src['joint_sequence_structure_relation']]
        elif src['presence']!='both_guides':category='guide_specific_domain'
        else:category='structural_direction_unresolved'
        assert x['candidate_class']==category and json.loads(x['gene_pair_key'])==[x['family'],x['gene_a'],x['gene_b']]
        screen=x['screen'],x['margin'];screens.add(screen)
        add(lineage[(*screen,category,x['study_role'],x['lineage_group'])],x);add(family[(*screen,category,x['family'])],x)
        if screen==('n30_c70','direction_margin_0_1') and category.startswith('stable_'):
            add(overview[category],x);roles[category][x['study_role']]+=1;lineage_overview[category].add(x['lineage_group'])
    assert seen==set(original) and len(seen)==r['annotated_rows']
    expected_lineage={(s,m,c,role,lin) for s,m in screens for c in categories for role,lin in lineage_taxa};seen_lineage=set()
    for x in csv.DictReader((root/'lineage_summary.tsv').open(),delimiter='\t'):
        k=tuple(x[f] for f in ['screen','margin','candidate_class','study_role','lineage_group']);assert k not in seen_lineage;seen_lineage.add(k);b=lineage.get(k,fresh())
        expected=dict(reconciled_taxa=len(lineage_taxa[k[3:]]),event_domain_combinations=b['n'],distinct_gene_pairs=len(b['pairs']),families=len(b['families']),pfams=len(b['pfams']),taxa_represented=len(b['taxa']))
        for f,v in expected.items():assert int(x[f])==v,(k,f)
    assert seen_lineage==expected_lineage and len(seen_lineage)==r['lineage_summary_rows']
    seen_family=set()
    for x in csv.DictReader((root/'family_summary.tsv').open(),delimiter='\t'):
        k=tuple(x[f] for f in ['screen','margin','candidate_class','family']);assert k not in seen_family and k in family;seen_family.add(k);b=family[k]
        for f,v in dict(event_domain_combinations=b['n'],distinct_gene_pairs=len(b['pairs']),taxa_represented=len(b['taxa']),pfams=len(b['pfams'])).items():assert int(x[f])==v
    assert seen_family==set(family) and len(seen_family)==r['family_summary_rows'] and len(taxa)==r['reconciled_taxa']
    expected_overview={c:dict(candidate_class=c,event_domain_combinations=b['n'],distinct_gene_pairs=len(b['pairs']),families=len(b['families']),pfams=len(b['pfams']),taxa=len(b['taxa']),manifest_lineage_groups=len(lineage_overview[c]),study_roles=dict(roles[c])) for c,b in overview.items()}
    assert {x['candidate_class']:x for x in r['n30_c70_margin_0_1_overview']}==expected_overview
    selected_lineages=[]
    for (screen,margin,c,role,lin),b in lineage.items():
        if screen=='n30_c70' and margin=='direction_margin_0_1' and c.startswith('stable_'):selected_lineages.append(dict(candidate_class=c,study_role=role,lineage_group=lin,event_domain_combinations=b['n'],gene_pairs=len(b['pairs']),taxa=len(b['taxa'])))
    result=dict(status='passed_full_domain_candidate_sampling_readback',producer_receipt_sha256=sha(rp),plan_sha256=sha(a.plan),checker_sha256=sha(__file__),annotated_rows_checked=len(seen),lineage_rows_checked=len(seen_lineage),family_rows_checked=len(seen_family),overview=list(expected_overview.values()),selected_lineages=sorted(selected_lineages,key=lambda x:(x['candidate_class'],-x['event_domain_combinations'],x['lineage_group'])),scope='Every copied source field and taxonomy assignment checked; exact event-domain/gene-pair/family/Pfam/taxon sets and all zero lineage rows independently reconstructed. Descriptive concentration only; no enrichment or phylogenetic correction.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
