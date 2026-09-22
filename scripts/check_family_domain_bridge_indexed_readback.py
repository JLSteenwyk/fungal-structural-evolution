#!/usr/bin/env python3
"""Check independent readback and rejection of changed links and memberships."""
import json
from pathlib import Path
import sqlite3
import tempfile
import check_family_domain_bridge as fixture
from readback_family_domain_bridge_indexed import audit, partitions, sha


original_build = fixture.build


def audited_build(plan_path):
    original_build(plan_path)
    source = json.loads(plan_path.read_text())
    root = plan_path.parent
    database = Path(source['output'])/'family_domain_bridge.sqlite'
    receipt_path = Path(source['output'])/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    for variant in ['valid','sequence','membership']:
        db = sqlite3.connect(database)
        if variant == 'sequence':
            db.execute("UPDATE proteins SET sequence_id='wrong' WHERE native_gene_id='0_0'")
        if variant == 'membership':
            db.execute("UPDATE proteins SET sequence_id='sequence_a' WHERE native_gene_id='0_0'")
            db.execute("DELETE FROM assignments WHERE guide='profile' AND native_gene_id='0_0'")
        db.commit(); db.close()
        receipt['bridge_sha256'] = sha(database)
        receipt_path.write_text(json.dumps(receipt))
        plan = dict(output=str(root/('audit_'+variant)),producer_pid=0,producer_start_ticks=0,
                    producer_plan=str(plan_path),pins={str(plan_path):sha(plan_path)})
        path = root/('audit_plan_'+variant)
        path.write_text(json.dumps(plan))
        try:
            audit(path)
        except ValueError as e:
            expected = {'sequence':'Changed architecture sequence/source link','membership':'Family membership differs'}
            assert variant in expected and str(e)==expected[variant], str(e)
        else:
            assert variant=='valid','Corruption was accepted'
    # Restore the fixture for its own post-build checks.
    db=sqlite3.connect(database)
    db.execute("INSERT INTO assignments VALUES ('profile','0_0','OG0000000')")
    db.commit(); db.close()


if __name__=='__main__':
    fixture.build=audited_build
    fixture.check()
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'clusters'
        for text in ['(mclmatrix\nbegin\n0 0_0 0_0 $\n)\n',
                     '(mclmatrix\nbegin\n0 0_0 $\n',
                     '(mclmatrix\nbegin\n1 0_0 $\n)\n']:
            path.write_text(text)
            try:
                list(partitions(path))
            except ValueError:
                pass
            else:
                raise AssertionError('Malformed partition accepted')
    print('PASS: full fixture readback; changed sequence and membership rejected; three malformed partitions rejected')
