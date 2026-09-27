#!/usr/bin/env python3
"""Build a normalized, provenance-bound domain/sequence triad dataset with explicit missing matches."""
import argparse,csv,hashlib,json,sqlite3
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed source: '+p)
    verify()
    for source in plan['sources']:
        receipt=json.loads(Path(source['receipt']).read_text())
        if receipt['status']!=source['status']:raise ValueError('Wrong source status')
        for path in source['files']:
            if sha(path)!=receipt['artifacts'][Path(path).name]:raise ValueError('Unbound source artifact')
    audit=json.loads(Path(plan['coverage_readback']).read_text())
    if audit['status']!='passed_complete_domain_coverage_screen_readback' or audit['producer_receipt_sha256']!=sha(plan['sources'][2]['receipt']):raise ValueError('Unbound coverage readback')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);db=sqlite3.connect(out/'domain_triads.sqlite')
    db.execute('PRAGMA journal_mode=DELETE');db.execute('PRAGMA temp_store=FILE')
    counts={}
    def import_table(table,paths,add_guide=False):
        count=0;fields=None
        for path in paths:
            with Path(path).open() as f:
                reader=csv.DictReader(f,delimiter='\t');cols=(['guide'] if add_guide else [])+reader.fieldnames
                if fields is None:
                    fields=cols;db.execute('CREATE TABLE '+table+' ('+','.join('"'+c+'" TEXT NOT NULL' for c in cols)+')')
                if fields!=cols:raise ValueError('Source columns differ')
                command='INSERT INTO '+table+' VALUES ('+','.join('?' for _ in fields)+')';batch=[]
                for row in reader:
                    values=([Path(path).name.split('_')[0]] if add_guide else [])+[row[c] for c in reader.fieldnames]
                    batch.append(values);count+=1
                    if len(batch)>=2000:db.executemany(command,batch);batch=[]
                if batch:db.executemany(command,batch)
        counts[table]=count;db.commit()
    import_table('triads',[plan['sources'][0]['files'][0]])
    import_table('domain_links',[plan['sources'][1]['files'][0]])
    import_table('coverage',[plan['sources'][2]['files'][0]])
    import_table('sequence_covariates',plan['sources'][3]['files'],True)
    event='guide,family,gene_node,gene_a,gene_b'
    db.executescript(f'''
    CREATE UNIQUE INDEX triads_key ON triads({event},reference_gene,policy);
    CREATE UNIQUE INDEX links_key ON domain_links(pair_key,policy,boundary,pfam_accession);
    CREATE UNIQUE INDEX coverage_key ON coverage(pair_key,mask);
    CREATE UNIQUE INDEX sequence_key ON sequence_covariates({event});
    CREATE TABLE boundaries(boundary TEXT PRIMARY KEY);
    INSERT INTO boundaries VALUES ('alignment'),('envelope');
    CREATE TABLE masks(mask TEXT PRIMARY KEY);
    INSERT INTO masks VALUES ('full'),('plddt70');
    ''')
    # A domain must be the SAME interval at a shared vertex of all three edges.
    def interval(alias,model,version):
        return f"CASE WHEN {alias}.model_left={model} AND {alias}.version_left={version} THEN {alias}.interval_left WHEN {alias}.model_right={model} AND {alias}.version_right={version} THEN {alias}.interval_right ELSE NULL END"
    da=interval('d','t.model_a','t.version_a');dbb=interval('d','t.model_b','t.version_b')
    aa=interval('a','t.model_a','t.version_a');ar=interval('a','t.reference_model','t.reference_version')
    bb=interval('b','t.model_b','t.version_b');br=interval('b','t.reference_model','t.reference_version')
    db.execute(f'''CREATE TABLE matched_domains AS
        SELECT t.*,d.boundary,d.pfam_accession,d.domain_pair_key AS duplicate_domain_pair_key,
        a.domain_pair_key AS a_reference_domain_pair_key,b.domain_pair_key AS b_reference_domain_pair_key,
        {da} AS interval_a,{dbb} AS interval_b,{ar} AS interval_reference
        FROM triads t JOIN domain_links d ON d.pair_key=t.duplicate_pair_key AND d.policy=t.policy
        JOIN domain_links a ON a.pair_key=t.a_reference_pair_key AND a.policy=t.policy AND a.boundary=d.boundary AND a.pfam_accession=d.pfam_accession
        JOIN domain_links b ON b.pair_key=t.b_reference_pair_key AND b.policy=t.policy AND b.boundary=d.boundary AND b.pfam_accession=d.pfam_accession
        WHERE {da}={aa} AND {dbb}={bb} AND {ar}={br}''')
    db.execute(f'CREATE UNIQUE INDEX matched_key ON matched_domains({event},reference_gene,policy,boundary,pfam_accession)')
    # Preserve a complete triad × policy × boundary disposition universe.
    on=' AND '.join('t.'+k+'=m.'+k for k in (event+',reference_gene,policy').split(','))
    db.execute(f'''CREATE VIEW domain_availability AS SELECT t.*,z.boundary,
        count(m.pfam_accession) AS common_domain_count,
        CASE WHEN count(m.pfam_accession)>0 THEN 'matched_single_copy_domain'
        WHEN t.triad_class='identical_model_in_triad' THEN 'identical_model_in_triad'
        ELSE 'no_common_single_copy_domain_candidate' END AS domain_status
        FROM triads t CROSS JOIN boundaries z LEFT JOIN matched_domains m ON {on} AND m.boundary=z.boundary
        GROUP BY {','.join('t.'+k for k in (event+',reference_gene,policy').split(','))},z.boundary''')
    son=' AND '.join('m.'+k+'=s.'+k for k in event.split(','))
    passed=','.join(f"(d.{sid}_pass='1' AND a.{sid}_pass='1' AND b.{sid}_pass='1') AS all_three_{sid}_pass" for sid in plan['screens'])
    seqfields=['duplicate_pair_sequence_distance','sequence_tip_a','sequence_tip_b','sequence_signed_difference','sequence_normalized_contrast','sequence_direction','sequence_covariate_status']
    seqselect=','.join('s.'+x for x in seqfields)
    db.execute(f'''CREATE VIEW triad_domain_coverage AS SELECT m.*,d.mask,{passed},
        {seqselect},s.chosen_reference_gene AS sequence_covariate_reference,
        CASE WHEN s.chosen_reference_gene=m.reference_gene THEN 'exact_reference' ELSE 'different_tied_reference' END AS sequence_reference_status,
        CASE WHEN s.chosen_reference_gene=m.reference_gene THEN s.distance_a_to_reference ELSE NULL END AS distance_a_to_this_reference,
        CASE WHEN s.chosen_reference_gene=m.reference_gene THEN s.distance_b_to_reference ELSE NULL END AS distance_b_to_this_reference
        FROM matched_domains m
        JOIN coverage d ON d.pair_key=m.duplicate_domain_pair_key
        JOIN coverage a ON a.pair_key=m.a_reference_domain_pair_key AND a.mask=d.mask
        JOIN coverage b ON b.pair_key=m.b_reference_domain_pair_key AND b.mask=d.mask
        JOIN sequence_covariates s ON {son}''')
    # An inner sequence join may not discard a matched event unnoticed.
    matched=db.execute('SELECT count(*) FROM matched_domains').fetchone()[0]
    unchecked=db.execute('''SELECT count(*) FROM triads t
        JOIN domain_links d ON d.pair_key=t.duplicate_pair_key AND d.policy=t.policy
        JOIN domain_links a ON a.pair_key=t.a_reference_pair_key AND a.policy=t.policy AND a.boundary=d.boundary AND a.pfam_accession=d.pfam_accession
        JOIN domain_links b ON b.pair_key=t.b_reference_pair_key AND b.policy=t.policy AND b.boundary=d.boundary AND b.pfam_accession=d.pfam_accession''').fetchone()[0]
    if unchecked!=matched:raise ValueError('Inconsistent interval or model at shared triad vertex')
    joined=db.execute('SELECT count(*) FROM triad_domain_coverage').fetchone()[0]
    if joined!=2*matched:raise ValueError('Missing sequence covariate or coverage source')
    availability=db.execute('SELECT count(*) FROM domain_availability').fetchone()[0]
    if availability!=2*counts['triads']:raise ValueError('Incomplete availability universe')
    summaries=[]
    for sid in plan['screens']:
        for guide,mask,total,passed_count in db.execute(f'SELECT guide,mask,count(*),sum(all_three_{sid}_pass) FROM triad_domain_coverage GROUP BY guide,mask'):
            summaries.append(dict(screen=sid,guide=guide,mask=mask,domain_triad_rows=total,all_three_pass=passed_count))
    counts.update(matched_domains=matched,domain_availability=availability,triad_domain_coverage=joined)
    statuses=[dict(guide=g,domain_status=s,rows=n) for g,s,n in db.execute('SELECT guide,domain_status,count(*) FROM domain_availability GROUP BY guide,domain_status')]
    db.execute('CREATE TABLE provenance(path TEXT PRIMARY KEY,sha256 TEXT NOT NULL)');db.executemany('INSERT INTO provenance VALUES (?,?)',plan['pins'].items())
    db.commit()
    if db.execute('PRAGMA integrity_check').fetchone()[0]!='ok':raise ValueError('Database integrity failed')
    db.close();verify()
    result=dict(status='complete_domain_triad_integration_pending_independent_readback',plan_sha256=ph,counts=counts,availability_statuses=statuses,screen_counts=summaries,artifacts={'domain_triads.sqlite':sha(out/'domain_triads.sqlite')},biological_inference_eligible=False,scope='Normalized joins, shared-vertex interval consistency, every triad/policy/boundary retained in availability view; matched domain/mask rows link all three comparisons and exact gene-oriented sequence covariates. Alternative guides/policies/boundaries/references are dependent. No biological effect test or ancestral/causal interpretation.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
