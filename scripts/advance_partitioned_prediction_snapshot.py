#!/usr/bin/env python3
"""Audit disjoint prediction partitions and convert them without merging provenance."""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from Bio import SeqIO
from audit_busco_gene_copies import ROOT, sha


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def link_counts(rows):
    return Counter(tuple(sorted(row.items())) for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    config_sha = sha(args.config)

    def verify():
        if sha(args.config) != config_sha:
            raise ValueError('Configuration changed')
        for name, expected in config['pinned_files'].items():
            if sha(ROOT / name) != expected:
                raise ValueError('Pinned file changed: ' + name)

    verify()
    output = ROOT / config['output']
    if output.exists():
        raise FileExistsError('Inspect previous execution; use immutable output paths')
    for cohort in config['cohorts']:
        if any((ROOT / cohort[key]).exists() for key in ['audit', 'conversion']):
            raise FileExistsError('Audit or conversion already exists')
    if shutil.disk_usage(ROOT).free < config['minimum_free_disk_bytes']:
        raise ValueError('Insufficient storage headroom')
    records = list(SeqIO.parse(ROOT / config['original_fasta'], 'fasta'))
    expected_ids = {r.id for r in records}
    if len(records) != len(expected_ids) or len(records) != config['expected_predictions']:
        raise ValueError('Original prediction universe differs')
    observed = set()
    for cohort in config['cohorts']:
        pred = ROOT / cohort['predictions']
        chunk = json.loads((pred / 'last_chunk.json').read_text())
        ids = {p.stem for p in pred.glob('S*.json')}
        if observed & ids or len(ids) != cohort['expected_predictions']:
            raise ValueError('Overlapping partitions or incorrect partition size')
        if (chunk['status'] != 'production_chunk_finished'
                or chunk['config_sha256'] != sha(pred / 'config.json')
                or chunk['oom_deferred'] != 0
                or chunk['new_predictions'] + chunk['cached_predictions'] != len(ids)):
            raise ValueError('Invalid partition completion receipt')
        if not cohort['partial_source'] and (chunk['interrupted'] or chunk['remaining_eligible']):
            raise ValueError('Split queue is incomplete')
        observed.update(ids)
    if observed != expected_ids:
        raise ValueError('Partitions do not cover exactly the original queue')
    output.mkdir(parents=True)
    (output / 'launch.json').write_text(json.dumps({
        'config_sha256': config_sha, 'started_unix': time.time(),
        'resource_plan': config['resource_plan']}, indent=2) + '\n')
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
    audited_ids = set()
    links = []
    receipts = []
    for cohort in config['cohorts']:
        verify()
        audit_command = [sys.executable, 'scripts/audit_local_predictions.py',
                         '--predictions', cohort['predictions'], '--inputs', cohort['inputs'],
                         '--output', cohort['audit']]
        if cohort['partial_source']:
            audit_command.append('--snapshot-live')
        conversion_command = [sys.executable, 'scripts/convert_esmfold_snapshot.py',
                              '--predictions', cohort['predictions'], '--audit', cohort['audit'],
                              '--output', cohort['conversion']]
        for stage, command in [('audit', audit_command), ('conversion', conversion_command)]:
            verify()
            print(cohort['name'], stage, flush=True)
            with (output / (cohort['name'] + '-' + stage + '.log')).open('w') as log:
                subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                               stderr=subprocess.STDOUT, check=True)
        audit_path = ROOT / cohort['audit'] / 'receipt.json'
        conversion_path = ROOT / cohort['conversion'] / 'receipt.json'
        audit = json.loads(audit_path.read_text())
        conversion = json.loads(conversion_path.read_text())
        rows = table(ROOT / cohort['audit'] / 'predictions.tsv')
        ids = {r['sequence_id'] for r in rows}
        if (len(ids) != len(rows) or audited_ids & ids
                or audit['predictions'] != cohort['expected_predictions']
                or len(ids) != cohort['expected_predictions']
                or conversion['status'] != 'complete_audited_local_model_conversion'
                or conversion['models'] != len(ids)
                or conversion['source_audit_receipt_sha256'] != sha(audit_path)):
            raise ValueError('Audit/conversion count or provenance mismatch')
        audited_ids.update(ids)
        links.extend(table(ROOT / cohort['audit'] / 'taxon_links.tsv'))
        receipts.append({'name': cohort['name'], 'predictions': len(ids),
                         'audit': cohort['audit'], 'audit_receipt_sha256': sha(audit_path),
                         'conversion': cohort['conversion'],
                         'conversion_receipt_sha256': sha(conversion_path)})
    original_links = [r for r in table(ROOT / config['original_links'])
                      if r['sequence_id'] in expected_ids]
    if audited_ids != expected_ids or link_counts(links) != link_counts(original_links):
        raise ValueError('Audited partition union or taxon-marker link multiplicity differs')
    verify()
    result = {'status': 'complete_partitioned_prediction_snapshot_handoff',
              'config_sha256': config_sha, 'predictions': len(audited_ids),
              'taxon_marker_links': len(links), 'cohorts': receipts,
              'finished_unix': time.time(),
              'interpretation': 'Exact original candidate union audited and converted; source configs remain separate. Previously reused models are outside this prediction queue. Native features and evolutionary analyses remain pending.'}
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
