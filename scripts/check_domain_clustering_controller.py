#!/usr/bin/env python3
"""Exercise the full domain controller CLI with a native two-domain database."""
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
from catalog_whole_proteome_structures import sha
from build_whole_proteome_foldseek_database import readback


def main():
    production = json.loads(Path('metadata/domain_clustering_plan.json').read_text())
    foldseek = production['foldseek']
    root = Path('results/domains/domain-coordinates-20260923-v1/shards')
    source = root / 'shard_00000.tar'
    records = {}
    with (root / 'shard_00000.jsonl').open() as f:
        for line in f:
            row = json.loads(line)
            if row['status'] == 'exported':
                records[row['member']] = row
    with tempfile.TemporaryDirectory() as temporary:
        tmp = Path(temporary); db = tmp / 'database'; db.mkdir()
        archive = tmp / 'fixture.tar'; models = []
        with tarfile.open(source) as src, tarfile.open(archive, 'w') as dst:
            for member in src.getmembers()[:2]:
                data = src.extractfile(member).read(); row = records[member.name]
                import hashlib
                assert hashlib.sha256(data).hexdigest() == row['pdb_sha256']
                dst.addfile(member, io.BytesIO(data))
                models.append({'path': member.name, 'length': row['residues'], 'sequence_sha256': row['fragment_sequence_sha256']})
        paths = tmp / 'archives.tsv'; paths.write_text(str(archive)+'\n')
        prefix = db / 'domains'
        subprocess.run([foldseek, 'createdb', str(paths), str(prefix), '--threads', '1', '--gpu', '0', '--mask-bfactor-threshold', '70', '--coord-store-mode', '1'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=dict(os.environ, CUDA_VISIBLE_DEVICES=''))
        checked = readback(prefix, models)
        database_plan = tmp / 'database-plan.json'; database_plan.write_text('{"fixture":true}\n')
        # This receipt is deliberately fixture-only: exercise gate wiring, not production certification.
        receipt = {'status': 'complete_domain_foldseek_database_with_full_sequence_and_coordinate_readback',
                   **checked, 'plan_sha256': sha(database_plan), 'excluded_rejected_intervals': 3,
                   'exported_intervals_with_missing_backbone': 1,
                   'artifacts': {p.name: sha(p) for p in db.iterdir() if p.is_file()}}
        rp = db / 'receipt.json'; rp.write_text(json.dumps(receipt))
        plan = dict(production)
        plan.update(database=str(db), database_plan=str(database_plan), output=str(tmp/'valid'), predecessor={'pid':os.getpid(), 'create_time':0}, pins={})
        plan['resources'] = dict(production['resources'], cpu=1, minimum_free_disk_gib=0, minimum_available_memory_gib=0, emergency_free_disk_gib=0)
        arguments = list(production['cluster_arguments'])
        arguments[arguments.index('--threads')+1] = '1'
        arguments[arguments.index('--split-memory-limit')+1] = '64G'
        plan['cluster_arguments'] = arguments
        pp = tmp / 'plan.json'
        def run():
            pp.write_text(json.dumps(plan))
            return subprocess.run([sys.executable, 'scripts/advance_domain_clustering.py', '--plan', str(pp)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=180)
        completed = run()
        if completed.returncode:
            logs = '\n'.join(p.read_text()[-8000:] for p in (tmp/'valid').glob('stage_*.log'))
            raise RuntimeError(completed.stderr + logs)
        actual = json.loads((tmp/'valid/receipt.json').read_text())
        assert actual['status'] == 'complete_domain_candidate_partition_membership_readback'
        assert actual['intervals'] == 2 and actual['excluded_rejected_intervals'] == 3 and actual['exported_intervals_with_missing_backbone'] == 1
        failures = []
        for label, altered in [('incomplete_status', dict(receipt,status='incomplete')), ('wrong_plan', dict(receipt,plan_sha256='0'*64)), ('changed_artifact', dict(receipt,artifacts={**receipt['artifacts'],'domains.lookup':'0'*64}))]:
            rp.write_text(json.dumps(altered)); plan['output'] = str(tmp/label)
            rejected = run()
            if rejected.returncode == 0 or (tmp/label/'stage_0.log').exists():
                raise AssertionError('Invalid completion gate reached native execution: '+label)
            failures.append(label)
    result = {'status':'passed_full_domain_clustering_cli_and_rejection_gates', 'native_fixture_intervals':2,
              'native_fixture_clusters':actual['clusters'], 'rejected_before_native_execution':failures,
              'controller_sha256':sha('scripts/advance_domain_clustering.py'), 'script_sha256':sha(__file__),
              'foldseek_sha256':sha(foldseek), 'scope':'Real native createdb, cluster and createtsv on two source-checked exported domains; full controller CLI including domains prefix, partition and exclusion-count propagation. Fixture-only database completion receipt. Rejects incomplete status, wrong plan binding and changed artifact before clustering. Does not validate production completion, every alignment threshold, homology or evolutionary effects.'}
    Path('metadata/domain_clustering_cli_fixture_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
