#!/usr/bin/env python3
"""Bind complete retained parallel software artifacts, including private failures."""
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    output=Path('metadata/weighted_parallel_retained_software_artifacts_20261004_v2.json')
    assert not output.exists();pins={};counts={};negatives=0;captures=0
    prep=Path('metadata/weighted_parallel_fixtures_software_preparation_20261004_v4.json')
    prepared=json.loads(prep.read_text());verify(prepared['source_hashes']);bind(pins,prep)
    for grid in prepared['grids']:
        for version in ['serial','parallel']:
            plan=Path(grid[version]);root=Path(json.loads(plan.read_text())['output'])
            for rp in [root/'receipt.json',root/'readback.json']:
                record=json.loads(rp.read_text());verify(record['source_hashes']);bind(pins,rp)
                for name,digest in record.get('artifacts',{}).items():bind(pins,root/name,digest)
                for name,digest in record.get('reader_checkpoint_artifacts',{}).items():bind(pins,root/name,digest)
                for path,digest in record['source_hashes'].items():bind(pins,path,digest)
    for item in prepared['negatives']:
        root=Path(json.loads(Path(item['plan']).read_text())['output'])
        assert not (root/'readback.json').exists()
        archive=root/'original-private-serialized-fixture';assert len(list(archive.iterdir()))==5
        failures=list((root/'failures').glob('*/failure.json'))
        for path in failures:
            assert not json.loads(path.read_text())['scientific_eligibility']
            assert (path.parent/'original_numeric_inputs.npz').exists();captures+=1
        files=[p for p in root.rglob('*') if p.is_file()]
        for path in files:bind(pins,path)
        counts[item['case']]=len(files);negatives+=1
    assert negatives==19 and captures==13
    for stem,version in [('weighted_parallel_qualification',2),('weighted_parallel_source_guards',1),
            ('weighted_parallel_fit_source_adapter',1)]:
        for kind in ['validation','execution','transport']:
            path=Path('metadata/'+stem+'_software_'+kind+'_20261004_v'+str(version)+'.json')
            record=json.loads(path.read_text());verify(record['source_hashes']);verify(record.get('artifacts',{}))
            bind(pins,path)
            for p,d in record['source_hashes'].items():bind(pins,p,d)
    bind(pins,Path(__file__));verify(pins)
    result=dict(status='verified_complete_retained_parallel_software_inputs_outputs_and_private_failures_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),private_negative_namespaces=negatives,
        private_original_fixture_archives=19,original_numeric_failure_captures=captures,
        private_case_file_counts=counts,bound_file_count=len(pins),source_hashes=pins,
        scientific_eligibility=False,
        scope='Every successful serial/parallel software producer/readback artifact and '
              'parallel reader checkpoint, every private negative namespace and five-file '
              'original fixture archive, thirteen exact numeric failure captures, and all '
              'three original software transport proofs rehashed. No production completion, '
              'native memory-corruption repair or biological acceptance.')
    with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','private_case_file_counts']},indent=2))


if __name__=='__main__':main()
