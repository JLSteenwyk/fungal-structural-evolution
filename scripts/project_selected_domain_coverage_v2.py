#!/usr/bin/env python3
"""Project all qualified domain configurations to selected biological records."""
import argparse,csv,json,time,itertools,subprocess
from pathlib import Path
import psutil,duckdb
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();dep=plan['producer']
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines());assert state=={'ActiveState':'inactive','ExecMainStatus':'0'}
    verify();bindings={}
    for name in ['inventory','qualification']:
        root=Path(plan[name]);r=json.loads((root/'receipt.json').read_text());audit=json.loads(Path(plan[name+'_audit']).read_text())
        assert audit['status']==plan[name+'_audit_status'] and audit['source_receipt_sha256']==sha(root/'receipt.json')
        for f,h in r['artifacts'].items():assert sha(root/f)==h
        bindings[name]=dict(receipt_sha256=sha(root/'receipt.json'),audit_sha256=sha(plan[name+'_audit']))
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    db=duckdb.connect(str(out/'projection.duckdb'));db.execute("SET threads=1");db.execute("SET memory_limit='12GB'");db.execute("SET max_temp_directory_size='20GB'")
    db.execute('CREATE TABLE nodes AS SELECT node_id,guide,family,taxon_id FROM read_json_auto(?)',[plan['target_nodes']])
    db.execute("CREATE TABLE selections AS SELECT s.*,n.guide,n.family,n.taxon_id FROM read_csv(?,delim='\t',header=true,all_varchar=true) s JOIN nodes n ON s.target_id=n.node_id",[str(Path(plan['inventory'])/'selection_domain_links.tsv.gz')])
    assert db.execute('SELECT count(*) FROM selections').fetchone()[0]==2786912
    assert db.execute('SELECT count(*) FROM (SELECT target_id,policy,scenario_id,count(*) n FROM selections GROUP BY ALL HAVING n<>1)').fetchone()[0]==0
    db.execute("CREATE TABLE qualification AS SELECT *,CAST(eligible_domain_count AS BIGINT) AS eligible FROM read_csv(?,delim='\t',header=true,all_varchar=true)",[str(Path(plan['qualification'])/'configuration_eligibility.tsv.gz')])
    assert db.execute('SELECT count(*) FROM qualification').fetchone()[0]==3777048
    db.execute('CREATE VIEW expanded AS SELECT s.*,q.boundary,q.mask_cohort,q.screen,q.eligible,q.eligibility_status FROM selections s JOIN qualification q USING(domain_config_id)')
    keys='guide,policy,scenario_id,boundary,mask_cohort,screen'
    statuses=['all_shared_domains_pass','some_shared_domains_pass','no_shared_domains_pass','no_shared_domain_comparison','identical_model_requires_separate_handling']
    status_sql=','.join("count(*) FILTER (WHERE eligibility_status='"+s+"') AS "+s for s in statuses)
    coverage_query='SELECT '+keys+''',count(*) selected_records,
        count(*) FILTER (WHERE eligible>0) usable_selected_records,
        sum(eligible) eligible_domain_occurrences,
        count(DISTINCT taxon_id) selected_taxa,
        count(DISTINCT family) selected_families,
        count(DISTINCT background_id) selected_backgrounds,
        count(DISTINCT taxon_id) FILTER (WHERE eligible>0) usable_taxa,
        count(DISTINCT family) FILTER (WHERE eligible>0) usable_families,
        count(DISTINCT background_id) FILTER (WHERE eligible>0) usable_backgrounds,
        '''+status_sql+' FROM expanded WHERE boundary=? AND mask_cohort=? AND screen=? GROUP BY '+keys
    reuse_query='SELECT '+keys+',max(n) maximum_usable_background_reuse FROM (SELECT '+keys+',background_id,count(*) n FROM expanded WHERE eligible>0 AND boundary=? AND mask_cohort=? AND screen=? GROUP BY '+keys+',background_id) GROUP BY '+keys
    for ix,condition in enumerate(itertools.product(['alignment','envelope'],['full','plddt70','both'],plan['screens'])):
        db.execute(('CREATE TABLE coverage AS ' if ix==0 else 'INSERT INTO coverage ')+coverage_query,list(condition))
        db.execute(('CREATE TABLE reuse_max AS ' if ix==0 else 'INSERT INTO reuse_max ')+reuse_query,list(condition))
        print('Projected setting',ix+1,'/ 36',condition,flush=True)
    db.execute('CREATE TABLE strata(guide VARCHAR,policy VARCHAR,scenario_id VARCHAR,boundary VARCHAR,mask_cohort VARCHAR,screen VARCHAR)')
    scenarios=[x['scenario_id'] for x in json.loads(Path(plan['scenarios']).read_text())];assert len(set(scenarios))==54
    grid=list(itertools.product(['profile','mafft'],plan['policies'],scenarios,['alignment','envelope'],['full','plddt70','both'],plan['screens']));assert len(grid)==15552
    db.executemany('INSERT INTO strata VALUES (?,?,?,?,?,?)',grid)
    metrics=['selected_records','usable_selected_records','eligible_domain_occurrences','selected_taxa','selected_families','selected_backgrounds','usable_taxa','usable_families','usable_backgrounds']+statuses
    query='SELECT '+','.join('s.'+k for k in keys.split(','))+','+','.join('coalesce(c.'+m+',0) AS '+m for m in metrics)+',coalesce(r.maximum_usable_background_reuse,0) AS maximum_usable_background_reuse FROM strata s LEFT JOIN coverage c USING('+keys+') LEFT JOIN reuse_max r USING('+keys+') ORDER BY '+','.join('s.'+k for k in keys.split(','))
    data=db.execute(query);fields=[x[0] for x in data.description];rows=data.fetchall();assert len(rows)==15552
    with (out/'selected_domain_coverage.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields);w.writerows(rows)
    assert sum(x[fields.index('selected_records')] for x in rows)==2786912*36
    for values in rows:
        row=dict(zip(fields,values));assert sum(row[s] for s in statuses)==row['selected_records'];assert row['usable_selected_records']==row['all_shared_domains_pass']+row['some_shared_domains_pass']
        assert row['usable_taxa']<=row['selected_taxa'] and row['usable_families']<=row['selected_families'] and row['usable_backgrounds']<=row['selected_backgrounds']
    db.close();verify()
    result=dict(status='complete_selected_domain_coverage_projection_pending_independent_readback',plan_sha256=ph,selected_records=2786912,expanded_record_cells=2786912*36,summary_cells=len(rows),source_bindings=bindings,duckdb_version=duckdb.__version__,artifacts={'selected_domain_coverage.tsv':sha(out/'selected_domain_coverage.tsv')},working_database='projection.duckdb',scope='All selections projected across boundaries, masks and six screens with explicit empty strata. Counts retain reused controls and guide/scenario alternatives. Taxa/families/backgrounds and maximum reuse shown before effect inference. Denominator is metadata-selected controls, not all original targets; unmatched targets remain in upstream selection inventory. No structural effect or causal claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
