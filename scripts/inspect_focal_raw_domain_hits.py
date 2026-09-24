#!/usr/bin/env python3
"""Check retained focal architectures against all stored precompetition hits."""
import argparse
import json
from pathlib import Path
import sqlite3
from inspect_focal_domain_architectures import sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--architectures', type=Path, required=True)
    ap.add_argument('--database', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    receipt = a.database.parent / 'receipt.json'
    paths = [a.architectures, a.database, receipt, Path(__file__),
             Path(__file__).with_name('inspect_focal_domain_architectures.py')]
    pins = {str(p): sha(p) for p in paths}
    assert json.loads(receipt.read_text())['artifacts'][a.database.name] == pins[str(a.database)]
    source = json.loads(a.architectures.read_text())
    assert source['status'] == 'complete_focal_and_sister_candidate_architecture_export'
    con = sqlite3.connect(a.database.resolve().as_uri() + '?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    records = []
    for protein in source['records']:
        sid = protein['sequence_id']
        assert con.execute('SELECT sequence_id FROM proteins WHERE taxon_id=? AND protein_id=?',
                           (protein['taxon_id'], protein['protein_id'])).fetchone()[0] == sid
        hits = [dict(r) for r in con.execute('SELECT * FROM hits WHERE sequence_id=? ORDER BY hit_number', (sid,))]
        assert len(hits) == protein['raw_hits']
        lookup = {h['hit_id']: h for h in hits}
        retained = set()
        for alternative in protein['candidate_architectures']['alternatives']:
            for hit in alternative['annotations']:
                assert all(str(v) == str(lookup[hit['hit_id']][k]) for k, v in hit.items())
                retained.add(hit['hit_id'])
        records.append(dict(tree_label=protein['tree_label'], sequence_id=sid,
                            raw_hits=hits, retained_in_any_policy=sorted(retained),
                            excluded_in_all_policies=sorted(set(lookup) - retained)))
    con.close()
    assert all(sha(p) == h for p, h in pins.items())
    out = dict(status='passed_focal_raw_hit_retained_annotation_join', source_hashes=pins,
               records=records, scope='All stored precompetition hits for the focal proteins. '
               'These hits already passed search reporting thresholds; missing hits do not prove '
               'domain absence. Annotation descriptions do not establish biological function.')
    with a.output.open('x') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    print(json.dumps({'proteins': len(records), 'raw_hits': sum(len(r['raw_hits']) for r in records)}))


if __name__ == '__main__':
    main()
