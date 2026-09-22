#!/usr/bin/env python3
"""Check independent variation metrics and complete-table corruption detection."""
import csv
import gzip
import json
from pathlib import Path
import random
import sqlite3
import subprocess
import sys
import tempfile
from check_family_architecture_variation import row
from summarize_family_architectures import summarize, POLICIES
from readback_family_architecture_variation import reconstruct, check_family, sha


def fixture():
    a, b = ['PF1', 'Domain'], ['PF2', 'Repeat']
    data = [row('F1','s1',[a,b]),row('F2','s2',[b,a]),row('F2','s3',[a,a,b]),
            row('F3','s4',[]),row('F4','s5',[a,b],partial=1),row('F5','s6',[a,b],agree=0)]
    result = reconstruct(data)[0]
    for key, value in dict(proteins=6, observed_ordered_signatures=3,
                           observed_multisets=2, observed_model_type_sets=1,
                           observed_multisets_with_order_variation=1,
                           observed_model_type_sets_with_multiplicity_variation=1,
                           conservative_proteins=3, no_ga_hit_proteins=1).items():
        assert result[key] == value, key
    rng = random.Random(21926)
    for case in range(100):
        records = []
        for i in range(rng.randrange(1,40)):
            r = list(row('T'+str(rng.randrange(6)), 's'+str(i),
                         [rng.choice([a,b,['PF1','Family']]) for _ in range(rng.randrange(6))],
                         agree=rng.randrange(2), partial=rng.randrange(2)))
            payload = json.loads(r[4])
            for p in POLICIES:
                for key in ['primary_rank_tie_pairs','unresolved_overlap_pairs','candidate_nested_pairs']:
                    payload['policies'][p][key] = rng.randrange(2)
            payload['alternatives'][0]['alignment_overlap_pairs'] = rng.randrange(2)
            r[4] = json.dumps(payload)
            records.append(tuple(r))
        assert reconstruct(records) == summarize(iter(records)), case
    expected = [dict(guide='profile',family='OG0',**r) for r in summarize(data)]
    valid = [{k:str(v) for k,v in r.items()} for r in expected]
    for variant in ['valid','changed','missing','repeated']:
        table = [dict(r) for r in valid]
        if variant == 'changed': table[0]['observed_multisets'] = '999'
        if variant == 'missing': table.pop()
        if variant == 'repeated': table[1] = table[0]
        try: check_family(iter(table),'profile','OG0',data)
        except ValueError: assert variant != 'valid'
        else: assert variant == 'valid'
    return data


def complete_fixture(data):
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp); bridge_db=root/'bridge.sqlite'; arch_db=root/'architecture.sqlite'
        db=sqlite3.connect(bridge_db)
        db.executescript('CREATE TABLE assignments(guide TEXT,native_gene_id TEXT,family TEXT); CREATE INDEX assignments_family ON assignments(guide,family); CREATE TABLE proteins(native_gene_id TEXT PRIMARY KEY,taxon_id TEXT,sequence_id TEXT);')
        arch=sqlite3.connect(arch_db)
        arch.execute('CREATE TABLE queries(sequence_id TEXT PRIMARY KEY,raw_hits INT,retained_sets_agree INT,candidate_architectures_json TEXT)')
        for i,(taxon,seq,hits,agree,payload) in enumerate(data):
            db.execute('INSERT INTO proteins VALUES (?,?,?)',(str(i),taxon,seq))
            for guide in ['profile','mafft']:db.execute('INSERT INTO assignments VALUES (?,?,?)',(guide,str(i),'OG'+str(i//3)))
            arch.execute('INSERT INTO queries VALUES (?,?,?,?)',(seq,hits,agree,payload))
        db.commit();db.close();arch.commit();arch.close()
        def save(name,d):
            path=root/name;path.write_text(json.dumps(d));return str(path)
        bridge=save('bridge.json',dict(proteins=6,guides=[dict(guide=g,families=2) for g in ['profile','mafft']],bridge_sha256=sha(bridge_db),architecture_database_sha256=sha(arch_db)))
        audit=save('bridge_audit.json',dict(status='passed_complete_independent_family_domain_bridge_readback',producer_receipt_sha256=sha(bridge)))
        folder=root/'producer';folder.mkdir()
        for guide in ['profile','mafft']:
            rr=[dict(guide=guide,family='OG'+str(i),**r) for i in range(2) for r in summarize(data[i*3:i*3+3])]
            with gzip.open(folder/(guide+'_family_architectures.tsv.gz'),'wt') as f:
                w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
        source=save('producer_plan.json',dict(output=str(folder),bridge_receipt=bridge,bridge_readback=audit,bridge_database=str(bridge_db),architecture_database=str(arch_db)))
        receipt=dict(status='complete_family_architecture_variation_inventory',input_hashes={source:sha(source)},artifacts={p.name:sha(p) for p in folder.iterdir()},guides=[dict(guide=g,families=2,proteins=6,rows=8) for g in ['profile','mafft']])
        (folder/'receipt.json').write_text(json.dumps(receipt))
        plan=save('audit_plan.json',dict(output=str(root/'readback'),producer_pid=0,producer_start_ticks=0,producer_plan=source,pins={source:sha(source)},resources=dict(minimum_free_disk_gib=0)))
        script=Path(__file__).with_name('readback_family_architecture_variation.py')
        subprocess.run([sys.executable,str(script),'--plan',plan],check=True)
        r=json.loads((root/'readback/receipt.json').read_text())
        assert r['status']=='passed_full_family_architecture_variation_readback'
        assert sum(x['rows'] for x in r['guides'])==16


if __name__ == '__main__':
    complete_fixture(fixture())
    print('PASS: known order/multiplicity counts, 100 varied cases, changed/missing/repeated row rejection, full two-guide fixture')
