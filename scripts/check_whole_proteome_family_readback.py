#!/usr/bin/env python3
"""Hand-calculated coverage fixture and semantic corruption checks."""
import csv
import json
from pathlib import Path
import sqlite3
import tempfile
from readback_whole_proteome_family_coverage import reconstruct


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        bridge, observed = root/'source.sqlite', root/'observed.sqlite'
        with sqlite3.connect(bridge) as c:
            c.execute('CREATE TABLE proteins(native_gene_id TEXT PRIMARY KEY,taxon_id TEXT,protein_id TEXT,sequence_id TEXT)')
            c.executemany('INSERT INTO proteins VALUES(?,?,?,?)', [('g1','t1','p1','Saaa'),('g2','t1','p2','Saaa'),('g3','t2','p3','Sbbb'),('g4','t3','p4','Sccc')])
            c.execute('CREATE TABLE assignments(guide TEXT,native_gene_id TEXT,family TEXT)')
            c.execute('CREATE INDEX assignments_family ON assignments(guide,family)')
            c.executemany('INSERT INTO assignments VALUES(?,?,?)', [('profile',g,f) for g,f in [('g1','f1'),('g2','f1'),('g3','f1'),('g4','f2')]])
        models = [('g1','t1','p1','aaa','m1',1,'/m1'),('g2','t1','p2','aaa','m1',1,'/m1'),('g3','t2','p3','bbb','m2',2,'/m2')]
        with sqlite3.connect(observed) as c:
            c.execute('CREATE TABLE structures(native_gene_id TEXT,taxon_id TEXT,protein_id TEXT,sequence_sha256 TEXT,model_id TEXT,version INTEGER,model_path TEXT)')
            c.executemany('INSERT INTO structures VALUES(?,?,?,?,?,?,?)',models)
        catalog = root/'links.tsv'
        with catalog.open('w') as h:
            w=csv.writer(h,delimiter='\t');w.writerow(['taxon_id','protein_id','sequence_sha256','model_id','version','model_path']);w.writerows(r[1:] for r in models)
        coverage = root/'coverage.tsv'
        text='guide\tfamily\tproteins\ttaxa\tstructure_proteins\tstructure_taxa\tstructure_sequences\tstructure_models\nprofile\tf1\t3\t2\t3\t2\t2\t2\nprofile\tf2\t1\t1\t0\t0\t0\t0\n'
        coverage.write_text(text)
        def run():return reconstruct(catalog,bridge,observed,coverage,['profile'])
        result=run()
        assert result['protein_links']==3 and result['source_proteins']==4
        assert result['guides'][0]['families_at_minimum_modeled_taxa']=={'2':1,'4':0,'10':0,'25':0,'50':0,'100':0}
        checks=[]
        def reject(name,expected):
            try:run()
            except ValueError as error:
                assert expected in str(error),str(error)
                checks.append(name)
            else:raise AssertionError(name+' was accepted')
        coverage.write_text(text.replace('3\t2\t3\t2\t2\t2','3\t2\t3\t3\t2\t2'))
        reject('paralog_count_as_taxa','Family coverage mismatch')
        coverage.write_text(text+text.splitlines(True)[-1]);reject('extra_family','Extra coverage rows')
        coverage.write_text('\n'.join(text.splitlines()[:-1])+'\n');reject('omitted_unmodeled_family','Family coverage mismatch')
        coverage.write_text(text)
        with sqlite3.connect(observed) as c:c.execute("UPDATE structures SET model_id='wrong' WHERE native_gene_id='g2'")
        reject('wrong_model_identity','Output model identity mismatch')
        with sqlite3.connect(observed) as c:c.execute("UPDATE structures SET model_id='m1' WHERE native_gene_id='g2'")
        with sqlite3.connect(bridge) as c:c.execute("UPDATE proteins SET sequence_id='Sbad' WHERE native_gene_id='g3'")
        reject('wrong_source_sequence','Source sequence mismatch')
        print(json.dumps({'status':'passed','rejected_corruptions':checks,'scope':'Hand-calculated source membership and model identity fixture; not full-data verification.'}))


if __name__=='__main__':main()
