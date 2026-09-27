#!/usr/bin/env python3
"""Quantify frozen model ascertainment for all reported terminal singleton-side duplications."""
import argparse,csv,gzip,itertools,json,sqlite3
from collections import Counter,defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

COLUMNS=['terminal_singleton_events','neither_model','one_model','both_models','both_same_model']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed coverage plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();cp=json.loads(Path(plan['coverage_plan']).read_text());coverage=Path(cp['output']);cr=json.loads((coverage/'receipt.json').read_text());ep=json.loads(Path(cp['event_plan']).read_text())
    if cr['status']!='complete_duplication_frozen_structure_coverage' or cr['plan_sha256']!=sha(plan['coverage_plan']):raise ValueError('Coverage binding differs')
    with Path(plan['sampling']).open() as f:taxa={r['taxon_id']:r for r in csv.DictReader(f,delimiter='\t')}
    db=sqlite3.connect('file:'+str(Path(cp['bridge']).resolve())+'?mode=ro',uri=True)
    models={t+'_'+p:m for t,p,m in db.execute('SELECT taxon_id,protein_id,model_id FROM structures')};db.close()
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);reports=[]
    taxon_counts={};family_counts={}
    with gzip.open(out/'terminal_event_coverage.tsv.gz','wt') as output:
        w=csv.writer(output,delimiter='\t',lineterminator='\n');w.writerow(['guide','family','taxon_id','gene_node','gene_a','gene_b','model_a','model_b','coverage_class','same_model'])
        for entry in ep['guides']:
            guide=entry['guide'];path=Path(entry['result'])/'Gene_Duplication_Events/Duplications.tsv';export=coverage/(guide+'_event_coverage.tsv.gz')
            if sha(path)!=ep['pins'][str(path)] or sha(export)!=cr['artifacts'][export.name]:raise ValueError('Changed event source/export')
            totals=Counter();tc={t:Counter() for t in taxa};fc=defaultdict(Counter);rows=0
            with path.open() as f,gzip.open(export,'rt') as g:
                original=csv.DictReader(f,delimiter='\t');exported=csv.DictReader(g,delimiter='\t')
                for event,old in itertools.zip_longest(original,exported):
                    if event is None or old is None:raise ValueError('Event export row count differs')
                    rows+=1
                    for x,y in [('Orthogroup','family'),('Species Tree Node','species_node'),('Gene Tree Node','gene_node'),('Support','support'),('Type','type')]:
                        if event[x]!=old[y]:raise ValueError('Event identity differs')
                    a_gene,b_gene=event['Genes 1'],event['Genes 2']
                    simple=event['Type']=='Terminal' and bool(a_gene) and bool(b_gene) and ', ' not in a_gene and ', ' not in b_gene
                    if int(old['terminal_singleton_sides'])!=int(simple):raise ValueError('Terminal-singleton classification differs')
                    if not simple:continue
                    taxon=event['Species Tree Node'];family=event['Orthogroup']
                    if taxon not in taxa or not all(g.startswith(taxon+'_') for g in [a_gene,b_gene]) or a_gene==b_gene:raise ValueError('Terminal taxon/gene identity differs')
                    left,right=models.get(a_gene,''),models.get(b_gene,'');n=bool(left)+bool(right);label=['neither_model','one_model','both_models'][n];same=int(bool(left) and left==right)
                    if [int(old[k]) for k in ['left_genes','right_genes','left_models','right_models','both_sides_covered']]!=[1,1,int(bool(left)),int(bool(right)),int(n==2)]:raise ValueError('Terminal model coverage differs')
                    for counts in [totals,tc[taxon],fc[family]]:counts['terminal_singleton_events']+=1;counts[label]+=1;counts['both_same_model']+=same
                    w.writerow([guide,family,taxon,event['Gene Tree Node'],a_gene,b_gene,left,right,label,same])
            expected=next(x for x in cr['guides'] if x['guide']==guide)
            if rows!=expected['counts']['events'] or totals['terminal_singleton_events']!=expected['counts']['terminal_singleton_sides'] or totals['both_models']!=expected['counts']['terminal_two_model_candidates'] or totals['both_same_model']!=expected['counts']['candidate_same_model']:raise ValueError('Coverage summary differs')
            taxon_counts[guide]=tc;family_counts[guide]=fc
            report=dict(guide=guide,all_source_events=rows,counts={k:totals[k] for k in COLUMNS},taxa_with_terminal_singletons=sum(c['terminal_singleton_events']>0 for c in tc.values()),taxa_with_both_models=sum(c['both_models']>0 for c in tc.values()),families_with_terminal_singletons=len(fc),families_with_both_models=sum(c['both_models']>0 for c in fc.values()))
            reports.append(report);print(json.dumps(report),flush=True)
    for level,mapping in [('taxon',taxon_counts),('family',family_counts)]:
        with (out/(level+'_coverage.tsv')).open('w') as f:
            fields=['guide',level]+(['species_name','study_role','lineage'] if level=='taxon' else [])+COLUMNS+['both_model_fraction']
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            for guide,groups in sorted(mapping.items()):
                for key,counts in sorted(groups.items()):
                    total=counts['terminal_singleton_events'];row=dict(guide=guide,**{level:key},**{k:counts[k] for k in COLUMNS},both_model_fraction=counts['both_models']/total if total else '')
                    if level=='taxon':row.update({k:taxa[key][k] for k in ['species_name','study_role','lineage']})
                    writer.writerow(row)
    verify();r=dict(status='complete_terminal_duplication_sampling_coverage',plan_sha256=ph,guides=reports,sampled_taxa=len(taxa),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All native reported terminal singleton-side events independently classified and joined to the frozen structural bridge, with existing export identities and singleton counts checked row by row. Taxon/family denominators include zero-model events; all sampled taxa retained. Not a complete duplication-event universe, biological event validation, missing-at-random assumption, causal effect or corrected evolutionary test. Frozen coverage excludes subsequent retrievals and ESMFold structures.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')


if __name__=='__main__':
    csv.field_size_limit(32*1024*1024);main()
