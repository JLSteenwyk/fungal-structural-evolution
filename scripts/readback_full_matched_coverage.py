#!/usr/bin/env python3
"""Rebuild every fixed selection, unmatched policy and attrition cell with independent SQL."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import sqlite3
import tempfile
from collections import Counter
from pathlib import Path
from full_matched_coverage_sources import load, MASKS, ATTRITION_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); source, bindings = load(plan, plan_path); original_bindings = dict(bindings)
    out = Path(plan['output']); rp = out / 'receipt.json'; receipt = json.loads(rp.read_text()); rh = sha(rp)
    assert receipt['status'] == 'complete_full_matched_coverage_pending_independent_readback' and receipt['plan_sha256'] == sha(plan_path) and receipt['scientific_eligibility'] is False
    assert receipt['source_hashes'] == original_bindings; bind(bindings, rp, rh)
    for name, digest in receipt['artifacts'].items(): bind(bindings, out / name, digest)
    verify(bindings)
    sid_ids = {s['scenario_id']: i for i, s in enumerate(source['scenarios'])}
    with tempfile.TemporaryDirectory(prefix='independent-matched-coverage-', dir=out) as directory:
        db = sqlite3.connect(Path(directory) / 'source.sqlite')
        db.execute('CREATE TABLE quality(kind TEXT,pair TEXT,mask TEXT,record TEXT,PRIMARY KEY(kind,pair,mask))')
        db.execute('CREATE TABLE nodes(kind TEXT,id TEXT,guide TEXT,full INTEGER,plddt70 INTEGER,both INTEGER,disposition TEXT,record TEXT,PRIMARY KEY(kind,id))')
        db.execute('CREATE TABLE policies(tid TEXT,policy TEXT,guide TEXT,full INTEGER,plddt70 INTEGER,both INTEGER,matched INTEGER,PRIMARY KEY(tid,policy))')
        db.execute('CREATE TABLE chosen(tid TEXT,policy TEXT,sid TEXT,bid TEXT,guide TEXT,tf INTEGER,tp INTEGER,tb INTEGER,bf INTEGER,bp INTEGER,bb INTEGER,PRIMARY KEY(tid,policy,sid))')
        for kind, key, count in [('target', 'target_quality', 'target_pairs'), ('background', 'background_quality', 'background_pairs')]:
            seen = set()
            with Path(source[key]).open() as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    assert row['mask'] in MASKS[:2]
                    for spec in plan['screens']: assert row[spec['id'] + '_pass'] in ['0', '1'] and (row[spec['id'] + '_pass'] == '1') == (row[spec['id'] + '_exclusions'] == '')
                    db.execute('INSERT INTO quality VALUES(?,?,?,?)', (kind, row['pair_key'], row['mask'], json.dumps(row))); seen.add((row['pair_key'], row['mask']))
            assert len(seen) == 2 * plan['expected'][count] and len({k[0] for k in seen}) == plan['expected'][count]
            assert not db.execute('SELECT pair FROM quality WHERE kind=? GROUP BY pair HAVING COUNT(*)<>2', (kind,)).fetchall()
        db.commit()
        for kind, key in [('target', 'targets'), ('background', 'backgrounds')]:
            for ident, node in source[key].items():
                endpoints = [(node['model_id_' + s], node['version_' + s]) for s in ['a', 'b']]
                assert node['pair_key'] == hashlib.sha256(json.dumps(sorted(endpoints), separators=(',', ':')).encode()).hexdigest()
                same = endpoints[0] == endpoints[1]; assert bool(node['same_model']) == same
                values = []
                if same: values = [0, 0]
                else:
                    for mask in MASKS[:2]:
                        record = db.execute('SELECT record FROM quality WHERE kind=? AND pair=? AND mask=?', (kind, node['pair_key'], mask)).fetchone(); assert record is not None; row = json.loads(record[0])
                        lengths = {(row['model_' + s], int(row['version_' + s])): int(row['length_' + s]) for s in ['a', 'b']}
                        assert set(lengths) == set(endpoints) and all(lengths[end] == node['length_' + s] for end, s in zip(endpoints, ['a', 'b']))
                        number = 0
                        for j, spec in enumerate(plan['screens']):
                            if int(row[spec['id'] + '_pass']): number |= 2**j
                        values.append(number)
                values.append(values[0] & values[1]); disposition = 'identical_model_no_alignment' if same else 'distinct_model_pair'
                db.execute('INSERT INTO nodes VALUES(?,?,?,?,?,?,?,?)', (kind, ident, node['guide'], *values, disposition, json.dumps(node)))
        db.commit()
        dispositions = {kind + ':' + guide + ':' + disposition: n for kind, guide, disposition, n in db.execute('SELECT kind,guide,disposition,COUNT(*) FROM nodes GROUP BY kind,guide,disposition')}
        statuses = 0; unmatched = 0; matched_expected = 0
        with (source['selection'] / 'target_policy_selection_status.tsv').open() as original, gzip.open(out / 'target_policy_coverage_status.tsv.gz', 'rt') as exported:
            for raw, actual in itertools.zip_longest(csv.DictReader(original, delimiter='\t'), csv.DictReader(exported, delimiter='\t')):
                assert raw is not None and actual is not None and raw['policy'] in plan['policies']
                result = db.execute('SELECT guide,full,plddt70,both,disposition,record FROM nodes WHERE kind=? AND id=?', ('target', raw['target_id'])).fetchone(); assert result is not None
                guide, full, masked, both, disposition, text = result; node = json.loads(text); assert guide in plan['guides']
                matched = raw['matched_scenarios'].split(',') if raw['matched_scenarios'] else []; absent = raw['unmatched_scenarios'].split(',') if raw['unmatched_scenarios'] else []
                assert len(matched) == len(set(matched)) and len(absent) == len(set(absent)) and not set(matched) & set(absent) and set(matched + absent) == set(sid_ids)
                bits = sum(2**sid_ids[sid] for sid in matched); unmatched += len(absent); matched_expected += len(matched)
                db.execute('INSERT INTO policies VALUES(?,?,?,?,?,?,?)', (raw['target_id'], raw['policy'], guide, full, masked, both, bits))
                expected = dict(raw, guide=guide, family=node['family'], focal_taxon=node['taxon_id'], gene_node=node['gene_node'], target_pair_key=node['pair_key'], target_comparison_disposition=disposition,
                                target_full_pass_bits=str(full), target_plddt70_pass_bits=str(masked), target_both_pass_bits=str(both))
                assert actual == expected; statuses += 1
                if statuses % 100000 == 0: print('Independent full matched policy states', statuses, '/', plan['expected']['target_policy_records'], flush=True)
        db.commit(); assert statuses == plan['expected']['target_policy_records'] == len(source['targets']) * len(plan['policies'])
        # Unique known target/policy keys of this cardinality cover their complete Cartesian grid.
        selected = 0
        with gzip.open(source['selection'] / 'selections.tsv.gz', 'rt') as original, gzip.open(out / 'selected_pair_coverage.tsv.gz', 'rt') as exported:
            for raw, actual in itertools.zip_longest(csv.DictReader(original, delimiter='\t'), csv.DictReader(exported, delimiter='\t')):
                assert raw is not None and actual is not None and raw['scenario_id'] in sid_ids and raw['endpoint_order'] in ['0', '1']
                key = raw['target_id'], raw['policy']; policy = db.execute('SELECT matched FROM policies WHERE tid=? AND policy=?', key).fetchone(); assert policy and policy[0] & 2**sid_ids[raw['scenario_id']]
                records = [db.execute('SELECT guide,full,plddt70,both,disposition,record FROM nodes WHERE kind=? AND id=?', (kind, ident)).fetchone() for kind, ident in [('target', raw['target_id']), ('background', raw['background_id'])]]
                assert all(records) and records[0][0] == records[1][0]; t, b = [json.loads(row[5]) for row in records]
                tf, tp, tb = records[0][1:4]; bf, bp, bb = records[1][1:4]
                expected = dict(raw, guide=t['guide'], family=t['family'], focal_taxon=t['taxon_id'], gene_node=t['gene_node'], target_pair_key=t['pair_key'], background_pair_key=b['pair_key'], target_comparison_disposition=records[0][4], background_comparison_disposition=records[1][4], target_sequence_distance=str(t['sequence_distance']), background_sequence_distance=str(b['sequence_distance']))
                for mask, target, control in zip(MASKS, [tf, tp, tb], [bf, bp, bb]): expected.update({f'target_{mask}_pass_bits':str(target), f'control_{mask}_pass_bits':str(control), f'joint_{mask}_pass_bits':str(target & control)})
                assert actual == expected
                db.execute('INSERT INTO chosen VALUES(?,?,?,?,?,?,?,?,?,?,?)', (*key, raw['scenario_id'], raw['background_id'], t['guide'], tf, tp, tb, bf, bp, bb)); selected += 1
                if selected % 100000 == 0: print('Independent full matched selection states', selected, '/', plan['expected']['selected_records'], flush=True)
        db.commit(); assert selected == matched_expected == plan['expected']['selected_records'] and unmatched == plan['expected']['unmatched_decisions']
        aggregates, baseline = {}, {}
        for mask, tcol, bcol in [('full','tf','bf'),('plddt70','tp','bp'),('both','tb','bb')]:
            for j, spec in enumerate(plan['screens']):
                query = f'SELECT guide,policy,sid,COUNT(*),SUM(({tcol}>>{j})&1),SUM(({bcol}>>{j})&1),SUM((({tcol}>>{j})&1)*(({bcol}>>{j})&1)) FROM chosen GROUP BY guide,policy,sid'
                for guide, policy, sid, n, tpass, bpass, joint in db.execute(query): aggregates[guide,policy,sid,mask,spec['id']] = n,tpass,bpass,joint
                query = f'SELECT guide,policy,COUNT(*),SUM(({mask}>>{j})&1) FROM policies GROUP BY guide,policy'
                for guide, policy, n, passed in db.execute(query): baseline[guide,policy,mask,spec['id']] = n,passed
        expected_rows = {}
        for guide, policy, sid, mask, spec in itertools.product(plan['guides'],plan['policies'],sid_ids,MASKS,plan['screens']):
            key = guide,policy,sid,mask,spec['id']; n, allpass = baseline.get((guide,policy,mask,spec['id']),(0,0)); selected_n,tpass,bpass,joint = aggregates.get(key,(0,0,0,0))
            prefix='|'.join([guide,policy,sid]); counts=source['matching_receipt']['counts'];assert n==counts.get(prefix+'|targets',0) and selected_n==counts.get(prefix+'|matched',0) and n-selected_n==counts.get(prefix+'|unmatched',0)
            assert 0<=allpass-tpass<=n-selected_n and min(joint,tpass-joint,bpass-joint,selected_n-tpass-bpass+joint)>=0
            row=dict(guide=guide,policy=policy,scenario_id=sid,mask=mask,screen=spec['id'],all_target_records=n,matched_records=selected_n,unmatched_records=n-selected_n,target_pass_all_records=allpass,target_pass_matched_records=tpass,target_pass_unmatched_records=allpass-tpass,control_pass_matched_records=bpass,joint_pass_matched_records=joint,target_only_pass_matched_records=tpass-joint,control_only_pass_matched_records=bpass-joint,neither_pass_matched_records=selected_n-tpass-bpass+joint)
            expected_rows[key]={k:str(v) for k,v in row.items()}
        attrition_rows=len(expected_rows)
        with (out/'matched_attrition_counts.tsv').open() as handle:
            reader=csv.DictReader(handle,delimiter='\t');assert reader.fieldnames==ATTRITION_FIELDS
            for actual in reader:
                key=tuple(actual[k] for k in ['guide','policy','scenario_id','mask','screen']);assert expected_rows.pop(key)==actual
        assert not expected_rows; db.close()
    summary=dict(target_nodes=len(source['targets']),background_nodes=len(source['backgrounds']),target_policy_records=statuses,scenarios=len(sid_ids),scenario_decisions=statuses*len(sid_ids),selected_records=selected,unmatched_decisions=unmatched,selection_screen_cells=selected*len(MASKS)*len(plan['screens']),full_scenario_screen_cells=statuses*len(sid_ids)*len(MASKS)*len(plan['screens']),attrition_rows=attrition_rows,screens=plan['screens'],masks=MASKS,guides=plan['guides'],policies=plan['policies'],node_dispositions=dispositions)
    assert all(receipt[key]==summary[key] for key in SUMMARY_FIELDS);verify(bindings)
    result=dict(status='passed_full_matched_coverage_sql_readback',plan_sha256=sha(plan_path),producer_receipt_sha256=rh,**summary,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2));return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();run(args.plan,args.output)
