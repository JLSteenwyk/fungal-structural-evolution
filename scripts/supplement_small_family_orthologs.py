#!/usr/bin/env python3
"""Write missing small-family pairs separately from immutable native output."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from assess_small_family_output_exposure import groups, sha

ROOT = Path(__file__).resolve().parents[1]


def expected_pairs(source):
    sizes = {}
    small = {}
    for family, genes in groups(source/'clusters_OrthoFinder.txt_id_pairs.txt'):
        sizes[family] = len(genes)
        if len(genes) in (2, 3):
            small[family] = genes
    expected = {(family, left, right) for family, genes in small.items()
                for left in genes for right in genes
                if left.split('_')[0] != right.split('_')[0]}
    needed = {gene for genes in small.values() for gene in genes}
    names = {}
    with (source/'SequenceIDs.txt').open() as f:
        for line in f:
            native, label = line.rstrip('\n').split(': ', 1)
            # Fail closed on unsupported label normalization/fallback modes.
            if label != label.strip() or len(label.split()) != 1 or any(c in label for c in ',:()'):
                raise ValueError('Sequence label requires a separately tested native naming mode')
            if native in needed:
                if native in names:
                    raise ValueError('Duplicate native identifier')
                names[native] = label
    if names.keys() != needed:
        raise ValueError('Missing sequence labels')
    species = {}
    with (source/'SpeciesIDs.txt').open() as f:
        for line in f:
            native, label = line.rstrip('\n').split(': ', 1)
            if native in species:
                raise ValueError('Duplicate species identifier')
            species[native] = label.rsplit('.', 1)[0]
    if len(set(species.values())) != len(species):
        raise ValueError('Ambiguous species labels')
    reverse = {}
    for native, label in names.items():
        key = (species[native.split('_')[0]], label)
        if key in reverse:
            raise ValueError('Ambiguous protein label within species')
        reverse[key] = native
    return expected, sizes, names, species, reverse


def scan_native(result, sizes, species, reverse, expected):
    observed = set()
    hashes = {}
    rows_scanned = 0
    for taxon in species.values():
        path = result/'Orthologues'/(taxon+'.tsv')
        stat = path.stat()
        digest = hashlib.sha256()
        with path.open('rb') as f:
            header = f.readline(); digest.update(header)
            if header.decode().rstrip('\r\n').split('\t') != ['Orthogroup','Species',taxon,'Orthologs']:
                raise ValueError('Unexpected grouped ortholog header: '+str(path))
            for raw in f:
                digest.update(raw); rows_scanned += 1
                family = raw.split(b'\t',1)[0].decode()
                if family not in sizes:
                    raise ValueError('Unknown output family')
                if sizes[family] > 3:
                    continue  # No expansion or validation of large-family pairs here.
                row = next(csv.reader([raw.decode()], delimiter='\t'))
                if len(row) != 4:
                    raise ValueError('Malformed small-family ortholog row')
                _, target, lefts, rights = row
                for left in lefts.split(', '):
                    for right in rights.split(', '):
                        pair = (family, reverse.get((taxon,left)), reverse.get((target,right)))
                        if pair not in expected or pair in observed:
                            raise ValueError('Unexpected or duplicate small-family pair')
                        observed.add(pair)
        after = path.stat()
        if (stat.st_size,stat.st_mtime_ns,stat.st_ino) != (after.st_size,after.st_mtime_ns,after.st_ino):
            raise ValueError('Native table changed during read')
        hashes[str(path)] = digest.hexdigest()
    return observed, hashes, rows_scanned


def supplement(source, result, output, completion_snapshot):
    output.mkdir(parents=True, exist_ok=False)
    inputs = {str(source/name):sha(source/name) for name in
              ['clusters_OrthoFinder.txt_id_pairs.txt','SequenceIDs.txt','SpeciesIDs.txt']}
    expected, sizes, names, species, reverse = expected_pairs(source)
    observed, native_hashes, scanned = scan_native(result, sizes, species, reverse, expected)
    missing = sorted(expected-observed)
    path = output/'supplemental_directed_pairs.tsv'
    with path.open('w') as f:
        writer = csv.writer(f,delimiter='\t',lineterminator='\n')
        writer.writerow(['family','source_native_gene','target_native_gene','source_taxon',
                         'target_taxon','source_protein','target_protein','provenance'])
        for family,left,right in missing:
            writer.writerow([family,left,right,species[left.split('_')[0]],species[right.split('_')[0]],
                             names[left],names[right],'supplement_small_family_cross_species_rule'])
    for filename, digest in inputs.items():
        if sha(filename) != digest:
            raise ValueError('Input changed during supplement')
    (output/'native_completion_snapshot.json').write_text(json.dumps(completion_snapshot,indent=2)+'\n')
    receipt = dict(status='complete_separate_small_family_ortholog_supplement',
        expected_directed_pairs=len(expected),native_directed_pairs=len(observed),
        supplemental_directed_pairs=len(missing),supplemented_families=len({x[0] for x in missing}),
        native_rows_scanned=scanned,input_hashes=inputs,native_ortholog_hashes=native_hashes,
        artifacts={path.name:sha(path),'native_completion_snapshot.json':sha(output/'native_completion_snapshot.json')},
        script_sha256=sha(__file__),cluster_parser_sha256=sha(Path(__file__).with_name('assess_small_family_output_exposure.py')),
        scope='Missing directed cross-species pairs in two-/three-gene families under the intended native small-family rule. Raw native outputs and statistics are unchanged. Larger families, biological orthology, HOGs, duplications and final combined statistics are not validated here.')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True,help='Completed native WorkingDirectory')
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--completion',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--execution-plan',type=Path,help='Required to bind production source inputs')
    a=p.parse_args()
    completion=json.loads(a.completion.read_text())
    if completion['status']=='completed_synthetic_small_family_output_observation':
        if a.results.resolve()!=a.completion.resolve().parent/'Results_small_fixture':
            raise ValueError('Fixture result binding mismatch')
        if a.source.resolve()!=a.completion.resolve().parent/'Source/WorkingDirectory':
            raise ValueError('Fixture source binding mismatch')
        for relative,digest in completion['source_hashes'].items():
            if sha(a.source.parent/relative)!=digest:raise ValueError('Fixture source changed')
        for relative,digest in completion['result_hashes'].items():
            if sha(a.results/relative)!=digest:raise ValueError('Fixture outputs changed')
    elif completion['status']=='both_native_executions_complete_pending_full_output_readback':
        if a.execution_plan is None or sha(a.execution_plan)!=completion['plan_sha256']:
            raise ValueError('Production execution plan is required and must match')
        plan=json.loads(a.execution_plan.read_text())
        if (ROOT/plan['output']).resolve()!=a.completion.resolve().parent:
            raise ValueError('Completion location differs from execution plan')
        inputs=ROOT/plan['inputs']
        if sha(inputs/'receipt.json')!=plan['pins'][str((inputs/'receipt.json').relative_to(ROOT))]:
            raise ValueError('Changed production input receipt')
        matches=[s for s in completion['stages'] if (ROOT/s['result']).resolve()==a.results.resolve()]
        if len(matches)!=1 or not matches[0]['source_inputs_unchanged']:
            raise ValueError('No completed matching native stage')
        if a.source.resolve()!=a.results.resolve().parent/'Source/WorkingDirectory':
            raise ValueError('Production source binding mismatch')
        input_receipt=json.loads((inputs/'receipt.json').read_text())
        guide=matches[0]['guide']
        manifest=inputs/guide/'copied_files.tsv'
        expected_manifest=next(g['manifest_sha256'] for g in input_receipt['guides'] if g['guide']==guide)
        if sha(manifest)!=expected_manifest:raise ValueError('Changed input manifest')
        records={row['relative_path']:row['sha256'] for row in csv.DictReader(manifest.open(),delimiter='\t')}
        for name in ['clusters_OrthoFinder.txt_id_pairs.txt','SpeciesIDs.txt','SequenceIDs.txt']:
            if sha(a.source/name)!=records['Source/WorkingDirectory/'+name]:
                raise ValueError('Production source input changed')
        for relative,digest in matches[0]['mandatory_artifacts'].items():
            if sha(a.results/relative)!=digest:raise ValueError('Native mandatory output changed')
    else:
        raise ValueError('Native execution is not complete')
    receipt=supplement(a.source,a.results,a.output,completion)
    print(json.dumps({k:receipt[k] for k in ['status','expected_directed_pairs','native_directed_pairs','supplemental_directed_pairs']},indent=2))


if __name__=='__main__':main()
