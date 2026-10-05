#!/usr/bin/env python3
"""Close complete domain reconstruction and retain all taxon denominators."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--taxon-table', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists() and not args.taxon_table.exists()
    plan_path = Path('metadata/completed_afdb_domain_registry_plan_20261005_v1.json')
    plan = json.loads(plan_path.read_text()); root = Path(plan['output'])
    producer_path = Path('metadata/completed_afdb_domain_registry_20261005_v1.json')
    reader_path = Path('metadata/completed_afdb_domain_registry_readback_20261005_v1.json')
    producer, reader = [json.loads(p.read_text()) for p in [producer_path, reader_path]]
    raw_path = root/'receipt.json'; raw = json.loads(raw_path.read_text())
    assert producer['status'] == 'completed_full_refreshed_afdb_domain_registry_pending_independent_readback'
    assert reader['status'] == 'passed_full_structure_domain_registry_readback'
    assert reader['producer_receipt_sha256'] == producer['raw_registry_receipt_sha256'] == sha(raw_path)
    assert reader['database_sha256'] == producer['database_sha256'] == raw['database_sha256']
    assert producer['counts'] == reader['counts'] == raw['counts']
    assert producer['candidate_domain_intervals_by_policy'] == reader['candidate_domain_intervals_by_policy'] == raw['candidate_domain_intervals_by_policy']
    assert producer['counts']['models'] == 2935733 and producer['counts']['protein_links'] == 2994868
    pins = {}; transports = []
    for prefix, validation, session in [('completed_afdb_domain_registry',producer_path,93892),
                                         ('completed_afdb_domain_registry_readback',reader_path,49632)]:
        ep = Path(f'metadata/{prefix}_execution_20261005_v1.json')
        tp = Path(f'metadata/{prefix}_transport_20261005_v1.json')
        pp = Path(f'metadata/{prefix}_original_tool_payloads_20261005_v1.json')
        lp = Path(f'metadata/{prefix}_launch_20261005_v1.json')
        e, t, tool, launch = [json.loads(p.read_text()) for p in [ep,tp,pp,lp]]
        assert e['exit_code'] == 0 and not e['timed_out']
        assert e['receipt_sha256'] == t['validation_sha256'] == sha(validation)
        assert t['original_tool_session_id'] == tool['original_tool_session_id'] == launch['original_tool_session_id'] == session
        assert tool['initial']['session_id'] == session and tool['terminal']['exit_code'] == 0
        assert e['wrapper'] == dict(pid=launch['pid'],created=launch['created'],cmdline=launch['cmdline'])
        assert e['invocation_id'] == launch['invocation_id'] == t['invocation_id']
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records'] == t['manager_completion_records'] == 1
        rows = [json.loads(line) for line in subprocess.check_output(
            ['journalctl','--user','-u',t['unit'],'--all','-o','json','--no-pager'],text=True).splitlines()]
        inv = e['invocation_id']; rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
        exact = [r for r in rows if r.get('_PID') == str(e['wrapper']['pid']) and r.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
        assert len(exact) == 2
        assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'],invocation_id=inv)
        assert json.loads(exact[1]['MESSAGE']) == {k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
        jp = Path(f'metadata/{prefix}_original_whole_journal_20261005_v1.jsonl')
        with jp.open('x') as handle:
            for row in rows: handle.write(json.dumps(row,sort_keys=True)+'\n')
        for mapping in [e['source_hashes'],e['artifacts'],t['source_hashes']]:
            for path,digest in mapping.items(): bind(pins,path,digest)
        for path in [ep,tp,pp,lp,jp]: bind(pins,path)
        transports.append(dict(unit=t['unit'],invocation_id=inv,original_tool_session_id=session,
                               actual_terminal_exit_code=0,whole_original_payloads_verified=True))
    for mapping in [producer['source_hashes'],reader['source_hashes'],plan['pins']]:
        for path,digest in mapping.items(): bind(pins,path,digest)
    coverage_path = Path('metadata/completed_afdb_catalog_taxon_coverage_20261005_v1.tsv')
    with coverage_path.open() as handle: coverage = list(csv.DictReader(handle,delimiter='\t'))
    assert len(coverage) == len({r['taxon_id'] for r in coverage}) == 526
    dbpath = root/'structure_domains.sqlite'; bind(pins,dbpath,raw['database_sha256'])
    db = sqlite3.connect('file:'+str(dbpath.resolve())+'?mode=ro',uri=True)
    by_taxon = {taxon:(n,none) for taxon,n,none in db.execute(
        'SELECT p.taxon_id,count(*),sum(m.raw_hits=0) FROM protein_links p JOIN models m USING(model_key) GROUP BY p.taxon_id')}
    policies = sorted(raw['candidate_domain_intervals_by_policy'])
    assert policies == ['alignment_bitscore','alignment_evalue','envelope_bitscore','envelope_evalue']
    candidates = {(taxon,policy):(intervals,proteins) for taxon,policy,intervals,proteins in db.execute(
        'SELECT p.taxon_id,h.policy,count(*),count(distinct p.protein_id) FROM protein_links p JOIN policy_hits h USING(model_key) WHERE h.domain_interval_candidate=1 GROUP BY p.taxon_id,h.policy')}
    summaries = []; totals = Counter()
    for r in coverage:
        n = int(r['representative_proteins']); k = int(r['proteins_with_model'])
        linked,without_hits = by_taxon.get(r['taxon_id'],(0,0)); assert linked == k
        assert 0 <= without_hits <= linked <= n
        row = dict(taxon_id=r['taxon_id'],species_name=r['species_name'],study_role=r['study_role'],
            representative_proteins=n,afdb_linked_proteins=linked,afdb_unlinked_proteins=n-linked,
            linked_proteins_with_raw_pfam_hits=linked-without_hits,linked_proteins_without_raw_pfam_hits=without_hits)
        for policy in policies:
            intervals,proteins = candidates.get((r['taxon_id'],policy),(0,0))
            assert 0 <= proteins <= linked and intervals >= proteins
            row[policy+'_candidate_interval_occurrences'] = intervals
            row[policy+'_proteins_with_candidate_interval'] = proteins
        summaries.append(row); totals.update({k:v for k,v in row.items() if isinstance(v,int)})
    assert totals['representative_proteins'] == 5815847 and totals['afdb_linked_proteins'] == 2994868
    assert totals['afdb_unlinked_proteins'] == 2820979
    for policy in policies:
        observed = totals[policy+'_candidate_interval_occurrences']
        expected = db.execute('SELECT count(*) FROM protein_links p JOIN policy_hits h USING(model_key) WHERE h.policy=? AND h.domain_interval_candidate=1',(policy,)).fetchone()[0]
        assert observed == expected
    db.close()
    for path in [coverage_path,plan_path,producer_path,reader_path,raw_path,Path(__file__)]: bind(pins,path)
    verify(pins)
    with args.taxon_table.open('x') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(summaries[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(summaries)
    bind(pins,args.taxon_table)
    result = dict(status='complete_verified_full_refreshed_afdb_domain_registry',
        checked_utc=datetime.now(timezone.utc).isoformat(),taxa=526,representative_proteins=5815847,
        counts=raw['counts'],candidate_domain_intervals_by_policy=raw['candidate_domain_intervals_by_policy'],
        public_taxon_table=str(args.taxon_table),protein_occurrence_totals=dict(totals),
        original_transports=transports,complete_bound_files=len(pins),source_hashes=pins,
        scientific_eligibility=False,gpu=False,new_predictions=0,all_eight_aims_incomplete=True,
        scope='Complete2,935,733model/2,994,868protein domain registry and every independent interval/policy '
              'and protein-link reconstruction close against actual original tool waits/whole journals/pins. '
              'All526taxon denominators and source-unlinked proteins retained. Taxon table counts candidate '
              'interval occurrences per linked protein; identical shared models can recur across proteins. '
              'Four policies are alternatives, not independent replicates or counts to sum. Raw Pfam '
              'absence is not lack of function. No residue-confidence/PAE, structural boundary, homology '
              'or evolutionary domain event acceptance; ESMFold integration remains separate.')
    with args.receipt.open('x') as handle: json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__ == '__main__': main()
