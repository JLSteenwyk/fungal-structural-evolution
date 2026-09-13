#!/usr/bin/env python3
"""Screen genus-label/code groups for codon coverage without inferring selection."""
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT,sha,read_table


def coverage(sequences):
    if not sequences:return {'taxa':0,'codon_columns':0,'fully_called_columns':0,'columns_at_least_80pct_called':0,'variable_columns_at_least_80pct_called':0}
    lengths={len(x) for x in sequences}
    if len(lengths)!=1 or next(iter(lengths))%3:raise ValueError('Nonrectangular codon alignment')
    full=covered=variable=0;n=len(sequences)
    for start in range(0,len(sequences[0]),3):
        called=[s[start:start+3] for s in sequences if set(s[start:start+3])<=set('ACGT')]
        full+=len(called)==n
        if len(called)>=math.ceil(.8*n):covered+=1;variable+=len(set(called))>1
    return {'taxa':n,'codon_columns':len(sequences[0])//3,'fully_called_columns':full,'columns_at_least_80pct_called':covered,'variable_columns_at_least_80pct_called':variable}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable coverage assessment')
    base=ROOT/'results/cds/marker-codon-alignments-v1';readback=ROOT/'results/cds/marker-codon-readback-v1/receipt.json'
    r=json.loads((base/'receipt.json').read_text());vr=json.loads(readback.read_text())
    if vr['status']!='complete_full_marker_codon_alignment_readback' or vr['source_receipt_sha256']!=sha(base/'receipt.json'):raise ValueError('Verified full codon projection required')
    for name,h in r['artifacts'].items():
        if sha(base/name)!=h:raise ValueError('Changed codon projection')
    mp=ROOT/'metadata/analysis_manifest.tsv';manifest={x['taxon_id']:x for x in read_table(mp)}
    lp=ROOT/'metadata/taxon_label_review.tsv';label_review={x['taxon_id'] for x in read_table(lp)}
    rows=read_table(base/'sequence_audit.tsv');accepted={(x['marker'],x['taxon_id']):x for x in rows if x['status']=='translation_verified_codon_alignment'}
    # Only fungal genus labels with at least four input entries are screened.
    groups=defaultdict(list)
    for t,m in manifest.items():
        if m['study_role']=='ingroup':groups[m['species_name'].split()[0]].append(t)
    groups={g:ts for g,ts in groups.items() if len(ts)>=4};membership=[];results=[]
    for g,ts in groups.items():
        for t in ts:membership.append({'genus_label':g,'taxon_id':t,'species_name':manifest[t]['species_name'],'lineage':manifest[t]['lineage'],'label_review_flag':t in label_review})
    for marker_record in r['markers']:
        marker=marker_record['marker'];sequences={x.id:str(x.seq) for x in SeqIO.parse(base/(marker+'.fna'),'fasta')}
        for genus,ts in groups.items():
            codes=sorted({accepted[marker,t]['translation_table'] for t in ts if (marker,t) in accepted})
            # Emit zero-coverage groups explicitly when all sequences were excluded.
            for code in codes or ['unavailable']:
                compatible=[t for t in ts if (marker,t) in accepted and accepted[marker,t]['translation_table']==code]
                for policy in ['all_aligned','exclude_recorded_annotation_gene_label_flags']:
                    retained=compatible if policy=='all_aligned' else [t for t in compatible if t not in label_review and accepted[marker,t]['gene_mapping_status']=='unique_gene' and accepted[marker,t]['representative_decision']=='longest_per_gene_lexical_tiebreak' and not accepted[marker,t]['annotation_flags']]
                    stats=coverage([sequences[t] for t in retained])
                    results.append({'marker':marker,'genus_label':genus,'translation_table':code,'policy':policy,'input_genus_entries':len(ts),'aligned_code_compatible_entries':len(compatible),**stats,'passes_coverage_screen':len(retained)>=4 and stats['columns_at_least_80pct_called']>=100,'retained_taxa_json':json.dumps(sorted(retained))})
        print(marker,len(results),flush=True)
    a.output.mkdir(parents=True)
    for name,data in [('group_codon_coverage.tsv',results),('genus_label_membership.tsv',membership)]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    summary=[]
    for policy in ['all_aligned','exclude_recorded_annotation_gene_label_flags']:
        for genus in sorted(groups):
            selected=[x for x in results if x['policy']==policy and x['genus_label']==genus]
            passing=[x for x in selected if x['passes_coverage_screen']]
            summary.append({'genus_label':genus,'policy':policy,'input_entries':len(groups[genus]),'marker_code_groups_screened':len(selected),'marker_code_groups_passing_coverage':len(passing),'distinct_markers_passing_coverage':len({x['marker'] for x in passing})})
    with (a.output/'group_summary.tsv').open('w') as f:
        w=csv.DictWriter(f,list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
    result={'status':'complete_full_marker_genus_code_coverage_screen','genus_labels':len(groups),'input_entries_in_groups':sum(map(len,groups.values())),'markers':len(r['markers']),'group_policy_rows':len(results),'passing_by_policy':dict(Counter(x['policy'] for x in results if x['passes_coverage_screen'])),'source_codon_receipt_sha256':sha(base/'receipt.json'),'source_readback_sha256':sha(readback),'manifest_sha256':sha(mp),'label_review_sha256':sha(lp),'script_sha256':sha(Path(__file__)),'thresholds':{'minimum_entries':4,'minimum_columns':100,'minimum_called_fraction':0.8},'interpretation':'Coverage screening, not selection eligibility or evidence of monophyly. Groups share a fungal genus label and translation code; entries are not assumed unique species or independent ecological transitions. The stricter policy removes only recorded annotation/gene/label flags, not all possible quality problems. Original protein masks are unchanged; canonical complete codons define called data. Tree support/reconciliation, divergence and saturation, alignment sensitivity, taxonomy and model adequacy still require assessment. No group is called a replicated transition.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))

if __name__=='__main__':main()
