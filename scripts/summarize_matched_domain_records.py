#!/usr/bin/env python3
"""Describe all matched domain contrasts with equal record, family and taxon weighting."""
import argparse,csv,json,itertools
from pathlib import Path
import duckdb
from screen_duplication_domain_alignment_coverage import sha

METRICS={
 'rmsd_difference':'rmsd_target_minus_background',
 'identity_difference':'sequence_identity_target_minus_background',
 'target_rmsd':'target_rmsd_recomputed','background_rmsd':'background_rmsd_recomputed',
 'target_identity':'target_sequence_identity_exact','background_identity':'background_sequence_identity_exact',
 'original_coverage_difference':'least(target_original_coverage_a,target_original_coverage_b)-least(background_original_coverage_a,background_original_coverage_b)',
 'log_aligned_length_ratio':'ln(target_aligned_length/background_aligned_length)',
 'confidence_fraction_difference':'target_joint_plddt70_fraction-background_joint_plddt70_fraction'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();bindings={}
    for name in ['measurements','inventory','coverage']:
        root=Path(plan[name]);r=json.loads((root/'receipt.json').read_text());proof=json.loads(Path(plan[name+'_audit']).read_text())
        assert proof['status']==plan[name+'_audit_status'] and proof['source_receipt_sha256']==sha(root/'receipt.json')
        for f,h in r['artifacts'].items():assert sha(root/f)==h
        bindings[name]=dict(receipt_sha256=sha(root/'receipt.json'),audit_sha256=sha(plan[name+'_audit']))
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);(out/'configuration_means').mkdir()
    db=duckdb.connect(str(out/'summaries.duckdb'));db.execute("SET threads=1");db.execute("SET memory_limit='12GB'");db.execute("SET max_temp_directory_size='20GB'")
    db.execute('CREATE TABLE nodes AS SELECT node_id,guide,family,taxon_id FROM read_json_auto(?)',[plan['target_nodes']])
    db.execute("CREATE TABLE selections AS SELECT s.*,n.guide,n.family,n.taxon_id FROM read_csv(?,delim='\t',header=true,all_varchar=true) s JOIN nodes n ON s.target_id=n.node_id",[str(Path(plan['inventory'])/'selection_domain_links.tsv.gz')])
    assert db.execute('SELECT count(*) FROM selections').fetchone()[0]==2786912
    # Parse numeric columns explicitly; blanks represent audited exclusions.
    numeric=['rmsd_target_minus_background','sequence_identity_target_minus_background','target_rmsd_recomputed','background_rmsd_recomputed','target_sequence_identity_exact','background_sequence_identity_exact','target_original_coverage_a','target_original_coverage_b','background_original_coverage_a','background_original_coverage_b','target_aligned_length','background_aligned_length','target_joint_plddt70_fraction','background_joint_plddt70_fraction']
    fixed=['domain_config_id','pfam_accession','boundary','mask','target_order','background_order','contrast_status']
    flags=[s+'_'+suffix for s in plan['screens'] for suffix in ['mask_pass','both_masks_pass']]
    sql='CREATE TABLE measurements AS SELECT '+','.join(fixed+flags+["CAST(NULLIF("+x+",'') AS DOUBLE) AS "+x for x in numeric])+" FROM read_csv(?,delim='\t',header=true,all_varchar=true)"
    db.execute(sql,[str(Path(plan['measurements'])/'matched_domain_contrasts.tsv.gz')]);assert db.execute('SELECT count(*) FROM measurements').fetchone()[0]==1200640
    db.execute('CREATE TABLE strata(guide VARCHAR,policy VARCHAR,scenario_id VARCHAR)')
    scenarios=[x['scenario_id'] for x in json.loads(Path(plan['scenarios']).read_text())];grid=list(itertools.product(['profile','mafft'],plan['policies'],scenarios));assert len(grid)==432
    db.executemany('INSERT INTO strata VALUES (?,?,?)',grid)
    coverage={}
    for x in csv.DictReader((Path(plan['coverage'])/'selected_domain_coverage.tsv').open(),delimiter='\t'):
        coverage[tuple(x[k] for k in ['guide','policy','scenario_id','boundary','mask_cohort','screen'])]=x
    keys='guide,policy,scenario_id';parts=[];total_rows=0
    fields=['guide','policy','scenario_id','boundary','mask','cohort','screen','target_order','background_order','selected_records','matched_records','families','taxa','backgrounds','mean_eligible_domains']
    for metric in METRICS:fields += [metric+'_record_mean',metric+'_family_equal_mean',metric+'_taxon_equal_mean']
    conditions=list(itertools.product(['alignment','envelope'],['full','plddt70'],['same_mask','both_masks'],plan['screens'],['0','1'],['0','1']));assert len(conditions)==192
    with (out/'record_summary.tsv').open('w') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for ix,(boundary,mask,cohort,screen,to,bo) in enumerate(conditions):
            flag=screen+('_mask_pass' if cohort=='same_mask' else '_both_masks_pass')
            query='SELECT domain_config_id,count(*) eligible_domains,'+','.join('avg('+expression+') AS '+metric for metric,expression in METRICS.items())+' FROM measurements WHERE boundary=? AND mask=? AND target_order=? AND background_order=? AND '+flag+"='1' GROUP BY domain_config_id"
            db.execute('CREATE OR REPLACE TABLE config_means AS '+query,[boundary,mask,to,bo])
            nc=db.execute('SELECT count(*) FROM config_means').fetchone()[0]
            parquet=out/'configuration_means'/f'{ix:03d}.parquet'
            # Paths originate from the pinned local plan; SQL parameters protect quoting.
            db.execute('COPY config_means TO ? (FORMAT PARQUET)',[str(parquet)])
            db.execute('CREATE OR REPLACE VIEW records AS SELECT s.*,c.* EXCLUDE(domain_config_id) FROM selections s JOIN config_means c USING(domain_config_id)')
            db.execute('CREATE OR REPLACE TABLE record_stats AS SELECT '+keys+',count(*) matched_records,count(DISTINCT family) families,count(DISTINCT taxon_id) taxa,count(DISTINCT background_id) backgrounds,avg(eligible_domains) mean_eligible_domains,'+','.join('avg('+m+') AS '+m+'_record_mean' for m in METRICS)+' FROM records GROUP BY '+keys)
            for grouping,label in [('family','family'),('taxon_id','taxon')]:
                db.execute('CREATE OR REPLACE TABLE '+label+'_stats AS SELECT '+keys+','+','.join('avg('+m+') AS '+m+'_'+label+'_equal_mean' for m in METRICS)+' FROM (SELECT '+keys+','+grouping+','+','.join('avg('+m+') AS '+m for m in METRICS)+' FROM records GROUP BY '+keys+','+grouping+') GROUP BY '+keys)
            query='SELECT s.*,r.* EXCLUDE(guide,policy,scenario_id),f.* EXCLUDE(guide,policy,scenario_id),t.* EXCLUDE(guide,policy,scenario_id) FROM strata s LEFT JOIN record_stats r USING('+keys+') LEFT JOIN family_stats f USING('+keys+') LEFT JOIN taxon_stats t USING('+keys+') ORDER BY s.guide,s.policy,s.scenario_id'
            cursor=db.execute(query);names=[x[0] for x in cursor.description];data=cursor.fetchall();assert len(data)==432
            for values in data:
                row=dict(zip(names,values));ck=(row['guide'],row['policy'],row['scenario_id'],boundary,mask if cohort=='same_mask' else 'both',screen);cv=coverage[ck]
                for count in ['matched_records','families','taxa','backgrounds']:row[count]=row[count] or 0
                assert row['matched_records']==int(cv['usable_selected_records']) and row['families']==int(cv['usable_families']) and row['taxa']==int(cv['usable_taxa']) and row['backgrounds']==int(cv['usable_backgrounds'])
                row.update(boundary=boundary,mask=mask,cohort=cohort,screen=screen,target_order=to,background_order=bo,selected_records=int(cv['selected_records']))
                writer.writerow(row);total_rows+=1
            parts.append(dict(index=ix,boundary=boundary,mask=mask,cohort=cohort,screen=screen,target_order=to,background_order=bo,configurations=nc,path=str(parquet.relative_to(out)),sha256=sha(parquet)))
            print('Summarized record setting',ix+1,'/ 192',flush=True)
    assert total_rows==82944;db.close();verify()
    (out/'partition_manifest.json').write_text(json.dumps(parts,indent=2)+'\n')
    result=dict(status='complete_matched_domain_record_summaries_pending_independent_readback',plan_sha256=ph,summary_rows=total_rows,settings=len(parts),metrics=METRICS,source_bindings=bindings,artifacts={n:sha(out/n) for n in ['record_summary.tsv','partition_manifest.json']},duckdb_version=duckdb.__version__,scope='Every selected-record sensitivity setting; equal averaging across eligible Pfam domains within each record, preserving all input orders/masks/boundaries. Record, family-equal and taxon-equal descriptive means shown separately. Every denominator/diversity count checked against full audited coverage. Empty strata retain blank means. No uncertainty interval, significance test, phylogenetic adjustment or causal duplication effect claimed.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
