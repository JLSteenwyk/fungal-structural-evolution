#!/usr/bin/env python3
"""Exercise cross-taxon ID collisions and reject incomplete domain mappings."""
import json
from pathlib import Path
import sqlite3
import tempfile
from build_family_domain_bridge import build, sha


def check(missing=False):
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p/'species').write_text('0: F1.faa\n1: F2.faa\n')
        (p/'sequences').write_text('0_0: repeated_id\n1_0: repeated_id\n')
        (p/'clusters').write_text('(mclmatrix\nbegin\n0 0_0 1_0 $\n)\n')
        (p/'members').write_text('0_0\n1_0\n')
        (p/'crosswalk').write_text('new_family\tsource_type\tsource_clade\tsource_family\tproteins\tmembership_sha256\nOG0000000\tretained\t\tOG0000000\t2\t'+sha(p/'members')+'\n')
        db = sqlite3.connect(p/'domains.sqlite')
        db.execute('CREATE TABLE proteins(taxon_id,protein_id,sequence_id,query_source)')
        db.execute("INSERT INTO proteins VALUES ('F1','repeated_id','sequence_a','fixture')")
        if not missing:
            db.execute("INSERT INTO proteins VALUES ('F2','repeated_id','sequence_b','fixture')")
        db.commit(); db.close()
        (p/'receipt').write_text(json.dumps({'artifacts':{'candidate_architectures.sqlite':sha(p/'domains.sqlite')}}))
        (p/'audit').write_text(json.dumps({'status':'passed_complete_independent_candidate_architecture_readback','production_receipt_sha256':sha(p/'receipt')}))
        (p/'merge').write_text('{}')
        (p/'merge_audit').write_text(json.dumps({'status':'passed_complete_guide_discovery_merge_readback','merge_receipt_sha256':sha(p/'merge')}))
        plan = dict(architecture_readback=str(p/'audit'), architecture_receipt=str(p/'receipt'),
                    architecture_database=str(p/'domains.sqlite'), merge_readback=str(p/'merge_audit'),
                    merge_receipt=str(p/'merge'), species_ids=str(p/'species'), sequence_ids=str(p/'sequences'),
                    output=str(p/'out'), taxa=2, proteins=2, resources={'minimum_free_disk_gib':0},
                    guides=[dict(name=g,clusters=str(p/'clusters'),crosswalk=str(p/'crosswalk'),families=1) for g in ['profile','mafft']],
                    pins={str(x):sha(x) for x in p.iterdir() if x.is_file()})
        (p/'plan').write_text(json.dumps(plan))
        try:
            build(p/'plan')
        except ValueError as e:
            if missing and str(e)=='Native and domain protein universes differ':
                return
            raise
        if missing:
            raise AssertionError('Accepted incomplete domain mapping')
        db = sqlite3.connect(p/'out/family_domain_bridge.sqlite')
        assert db.execute('SELECT native_gene_id,sequence_id FROM proteins ORDER BY native_gene_id').fetchall()==[('0_0','sequence_a'),('1_0','sequence_b')]
        assert db.execute('SELECT COUNT(*) FROM assignments').fetchone()[0]==4
        db.close()


if __name__=='__main__':
    check()
    check(missing=True)
    print('PASS: same protein ID in different taxa stays distinct; missing domain mapping rejected')
