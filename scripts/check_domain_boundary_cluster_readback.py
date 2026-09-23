#!/usr/bin/env python3
"""Exercise independent full-table boundary reconstruction and corruptions."""
import csv
import gzip
import json
from pathlib import Path
import tempfile
from readback_domain_boundary_clusters import audit_tables, sha


def table(path, rows):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def main():
    rejected = []
    with tempfile.TemporaryDirectory(prefix='domain-boundary-readback-') as tmp:
        root = Path(tmp); lookup = root/'domains.lookup'; members = root/'members.tsv'; output = root/'result.tsv.gz'
        intervals = [dict(interval_id=i, model_key='m') for i in ['a','b','c','x','y']]
        pairs = [('a','a','identical_interval'), ('a','b','distinct_intervals_same_cluster'),
                 ('a','c','distinct_intervals_different_clusters'), ('x','a','alignment_unclustered'),
                 ('a','x','envelope_unclustered'), ('x','y','both_unclustered')]
        assignments = {'a':'a','b':'a','c':'c'}
        links = [dict(model_key='m', hit_id=str(i), boundary=kind, interval_id=interval)
                 for i,(a,e,_) in enumerate(pairs) for kind,interval in [('alignment',a),('envelope',e)]]
        rows = [dict(model_key='m',hit_id=str(i),alignment_interval=a,envelope_interval=e,
                     alignment_cluster=assignments.get(a,''),envelope_cluster=assignments.get(e,''),disposition=d)
                for i,(a,e,d) in enumerate(pairs)]
        def reset():
            table(root/'intervals.tsv',intervals);table(root/'boundary_links.tsv.gz',links);table(output,rows)
            lookup.write_text('0\ta\t0\n1\tb\t0\n2\tc\t0\n');members.write_text('a\ta\na\tb\nc\tc\n')
        def check():return audit_tables(root,lookup,members,output)
        def reject(label, mutate):
            reset();mutate()
            try:check()
            except (ValueError,KeyError):rejected.append(label)
            else:raise AssertionError('Accepted corruption: '+label)
        reset();counts=check()
        assert counts==dict(candidate_pairs=6,intervals=5,clustered_intervals=3,unclustered_intervals=2,boundary_links=12,dispositions={p[2]:1 for p in pairs})
        for field in rows[0]:
            reject('altered_output_'+field, lambda field=field:table(output,[dict(rows[0],**{field:'incorrect'})]+rows[1:]))
        reject('missing_output_pair',lambda:table(output,rows[:-1]))
        reject('duplicate_output_pair',lambda:table(output,rows+rows[:1]))
        reject('missing_source_boundary',lambda:table(root/'boundary_links.tsv.gz',links[:-1]))
        reject('duplicate_source_boundary',lambda:table(root/'boundary_links.tsv.gz',links+links[:1]))
        reject('wrong_source_model',lambda:table(root/'intervals.tsv',[dict(intervals[0],model_key='other')]+intervals[1:]))
        reject('unknown_source_boundary',lambda:table(root/'boundary_links.tsv.gz',[dict(links[0],boundary='unknown')]+links[1:]))
        reject('duplicate_lookup',lambda:lookup.write_text('0\ta\t0\n1\ta\t0\n2\tc\t0\n'))
        reject('missing_partition_member',lambda:members.write_text('a\ta\nc\tc\n'))
        reject('nonself_representative',lambda:members.write_text('b\ta\na\tb\nc\tc\n'))
        reject('duplicate_partition_member',lambda:members.write_text('a\ta\na\ta\na\tb\nc\tc\n'))
    result=dict(status='passed_independent_boundary_cluster_table_fixtures',categories=6,
                rejected_cases=len(rejected),rejected=rejected,script_sha256=sha(__file__),
                readback_script_sha256=sha('scripts/readback_domain_boundary_clusters.py'),
                scope='Synthetic complete-table reconstruction and malformed source/output tests. Production output readback remains pending; not a controller lifecycle or biological validation test.')
    Path('metadata/domain_boundary_cluster_readback_fixture_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
