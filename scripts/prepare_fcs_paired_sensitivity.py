#!/usr/bin/env python3
"""Prepare all affected paired markers after omitting FCS EXCLUDE/FIX/TRIM observations."""
import argparse,csv,json,hashlib
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from prepare_paired_phylogenetic_inputs import site_counts


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['inputs','mapping','audit','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable output')
    source=json.loads((a.inputs/'receipt.json').read_text());mapping=json.loads((a.mapping/'receipt.json').read_text());audit=json.loads((a.audit/'receipt.json').read_text())
    if audit['status']!='passed_full_region_cds_intersection_and_protein_marker_join_audit' or audit['mapping_receipt_sha256']!=sha(a.mapping/'receipt.json'):raise ValueError('Missing matching overlap audit')
    for root,receipt in [(a.inputs,source),(a.mapping,mapping)]:
        for name,digest in receipt['artifacts'].items():
            if sha(root/name)!=digest:raise ValueError('Changed source artifact')
    flagged=defaultdict(dict)
    for row in table(a.mapping/'marker_overlap_review.tsv'):
        if {'EXCLUDE','FIX','TRIM'} & set(row['fcs_actions'].split(';')):flagged[row['marker']][row['taxon_id']]=row
    original_summary=table(a.inputs/'marker_summary.tsv');summaries=[];dispositions=[];omissions=[];changes=[];a.output.mkdir(parents=True)
    for summary in original_summary:
        marker=summary['marker'];new=dict(summary);before=int(summary['eligible_taxa']);cols=int(summary['retained_columns']);remove=set();keptcols=cols;after=before
        if summary['status']!='ready_for_inference':status='baseline_not_ready_no_refit'
        else:
            sequences={label:{r.id:str(r.seq) for r in SeqIO.parse(a.inputs/marker/(label+'.faa'),'fasta')} for label in ['aa','3di']}
            taxa=set(sequences['aa'])
            if taxa!=set(sequences['3di']) or len(taxa)!=before:raise ValueError('Baseline taxon grid differs')
            remove=taxa & set(flagged[marker]);retained=taxa-remove;after=len(retained)
            if not remove:status='unchanged_reuse_baseline_fit'
            else:
                for t in sorted(remove):omissions.append(dict(flagged[marker][t],baseline_input_receipt_sha256=sha(a.inputs/'receipt.json')))
                keep=[i for i in range(cols) if any(sequences['aa'][t][i]!='?' for t in retained)];keptcols=len(keep)
                if after<4 or not keep:status='insufficient_taxa_or_columns_after_omission'
                else:
                    status='ready_for_inference';folder=a.output/marker;folder.mkdir()
                    for label in ['aa','3di']:
                        projected={t:''.join(sequences[label][t][i] for i in keep) for t in sorted(retained)}
                        with (folder/(label+'.faa')).open('w') as f:
                            for t,seq in projected.items():f.write('>'+t+'\n'+seq+'\n')
                        loaded={r.id:str(r.seq) for r in SeqIO.parse(folder/(label+'.faa'),'fasta')}
                        if loaded!=projected:raise ValueError('Projected FASTA readback differs')
                        variable,informative=site_counts(list(loaded.values()));new[label+'_variable_columns']=variable;new[label+'_parsimony_informative_columns']=informative
                    aa={r.id:str(r.seq) for r in SeqIO.parse(folder/'aa.faa','fasta')};ss={r.id:str(r.seq) for r in SeqIO.parse(folder/'3di.faa','fasta')}
                    if set(aa)!=retained or set(ss)!=retained:raise ValueError('Omission taxon grid differs')
                    for t in retained:
                        if [c=='?' for c in aa[t]]!=[c=='?' for c in ss[t]]:raise ValueError('Paired masks differ')
                        if sum(c!='?' for c in aa[t])!=sum(c!='?' for c in sequences['aa'][t]):raise ValueError('Retained observations changed')
                    columns=table(a.inputs/marker/'columns.tsv')
                    if len(columns)!=cols or [int(x['paired_column_1based']) for x in columns]!=list(range(1,cols+1)):raise ValueError('Baseline column grid differs')
                    projected=[]
                    for j,i in enumerate(keep):
                        row=dict(columns[i]);row['baseline_paired_column_1based']=row['paired_column_1based'];row['paired_column_1based']=j+1;projected.append(row)
                    write(folder/'columns.tsv',projected)
                    for x in table(folder/'columns.tsv'):
                        old=columns[int(x['baseline_paired_column_1based'])-1]
                        if any(x[k]!=v for k,v in old.items() if k!='paired_column_1based'):raise ValueError('Source coordinate mapping differs')
                    changes.append({'marker':marker,'retained_taxa':after,'retained_columns':keptcols,'removed_observations':len(remove)})
        new.update(eligible_taxa=after,retained_columns=keptcols,status=status);summaries.append(new)
        dispositions.append({'marker':marker,'baseline_status':summary['status'],'sensitivity_status':status,'baseline_eligible_taxa':before,'sensitivity_eligible_taxa':after,'omitted_taxa':';'.join(sorted(remove)),'baseline_columns':cols,'sensitivity_columns':keptcols,'all_missing_columns_removed':cols-keptcols})
    if not omissions:raise ValueError('No affected observations in this dataset')
    write(a.output/'marker_summary.tsv',summaries);write(a.output/'marker_disposition.tsv',dispositions);write(a.output/'omitted_marker_observations.tsv',omissions)
    result={'status':'complete_affected_marker_fcs_omission_sensitivity_inputs','ready_markers':len(changes),'baseline_ready_markers':source['ready_markers'],'unchanged_ready_markers_reuse_baseline':sum(x['sensitivity_status']=='unchanged_reuse_baseline_fit' for x in dispositions),'affected_markers_below_eligibility':sum(x['sensitivity_status']=='insufficient_taxa_or_columns_after_omission' for x in dispositions),'omitted_marker_observations':len(omissions),'source_input_receipt_sha256':sha(a.inputs/'receipt.json'),'source_mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'source_audit_receipt_sha256':sha(a.audit/'receipt.json'),'script_sha256':sha(Path(__file__)),'site_statistics_helper_sha256':sha(Path(__file__).with_name('prepare_paired_phylogenetic_inputs.py')),'changes':changes,'interpretation':'Sensitivity refit inputs for every changed baseline-ready marker; unchanged estimates remain in the baseline dataset and are not rerun. Omit entire marker/taxon observations with CDS overlap to FCS EXCLUDE/FIX/TRIM, retain REVIEW-only observations. Other taxa and observed AA/3Di characters unchanged; drop only now-all-missing columns with exact coordinate projection. Existing confidence masks and taxon coverage eligibility inherited; at least four remaining taxa. Not confirmed contamination, whole-species deletion or completed sensitivity inference.','artifacts':{str(p.relative_to(a.output)):sha(p) for p in a.output.rglob('*') if p.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['artifacts','changes']},indent=2))


if __name__=='__main__':main()
