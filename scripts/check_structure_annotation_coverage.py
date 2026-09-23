#!/usr/bin/env python3
"""Full-CLI fixtures for annotation-stratified availability accounting."""
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile


def main():
    script=Path(__file__).with_name('summarize_structure_annotation_coverage.py').resolve()
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);ann=root/'annotations.sqlite';bridge=root/'bridge.sqlite'
        c=sqlite3.connect(ann);c.executescript("CREATE TABLE proteins(taxon_id,protein_id,sequence_id); CREATE TABLE queries(sequence_id,raw_hits); INSERT INTO proteins VALUES ('T','p1','Saaa'),('T','p2','Sbbb'),('T','p3','Sccc'); INSERT INTO queries VALUES ('Saaa',1),('Sbbb',0),('Sccc',0);");c.commit();c.close()
        c=sqlite3.connect(bridge);c.executescript("CREATE TABLE structures(native_gene_id,taxon_id,protein_id); INSERT INTO structures VALUES ('0_1','T','p1'),('0_2','T','p2');");c.commit();c.close()
        links=root/'links.tsv';links.write_text('taxon_id\tprotein_id\tsequence_sha256\nT\tp1\taaa\nT\tp2\tbbb\n')
        coverage=root/'coverage.tsv';coverage.write_text('taxon_id\trepresentative_proteins\tproteins_with_model\nT\t3\t2\n')
        plan={'architectures':str(ann),'structure_bridge':str(bridge),'links':str(links),'taxon_coverage':str(coverage),'output':str(root/'out'),'expected_proteins':3,'expected_links':2,'pins':{}}
        config=root/'plan.json';config.write_text(json.dumps(plan))
        result=subprocess.run([sys.executable,str(script),'--plan',str(config)],capture_output=True,text=True)
        assert result.returncode==0,result.stderr
        r=json.loads((root/'out'/'receipt.json').read_text());assert r['pooled']['with_pfam_hit']['modeled_fraction']==1 and r['pooled']['without_pfam_hit']['modeled_fraction']==0.5
        links.write_text('taxon_id\tprotein_id\tsequence_sha256\nT\tp1\twrong\nT\tp2\tbbb\n');plan['output']=str(root/'bad');config.write_text(json.dumps(plan))
        result=subprocess.run([sys.executable,str(script),'--plan',str(config)],capture_output=True,text=True)
        assert result.returncode!=0 and 'sequence mismatch' in result.stderr
    print('Passed full-CLI strata counts, independent aggregation and mismatched-sequence rejection')


if __name__=='__main__':main()
