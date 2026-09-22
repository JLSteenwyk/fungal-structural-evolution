#!/usr/bin/env python3
"""Observe installed OrthoFinder output coverage on synthetic small families."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot(path):
    return {str(p.relative_to(path)): sha(p) for p in path.rglob('*') if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--early-singleton', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = out/'Source'
    wd = source/'WorkingDirectory'
    trees = wd/'Trees_ids'
    trees.mkdir(parents=True)
    # Test both size ordering and a singleton preceding small non-singletons.
    species_patterns = [(0,1,2,3), (0,1,2), (0,0,2), (1,1,1), (0,2), (3,3), (2,)]
    if args.early_singleton:
        species_patterns = [species_patterns[0], species_patterns[-1], *species_patterns[1:-1]]
    used = [0]*4
    families = {}
    for index, pattern in enumerate(species_patterns):
        genes = []
        for taxon in pattern:
            genes.append(f'{taxon}_{used[taxon]}')
            used[taxon] += 1
        families[f'OG{index:07d}'] = genes
    (wd/'SpeciesIDs.txt').write_text(''.join(f'{i}: Taxon{i}.faa\n' for i in range(4)))
    (wd/'SequenceIDs.txt').write_text(''.join(f'{i}_{g}: protein_{i}_{g}\n' for i in range(4) for g in range(used[i])))
    for i in range(4):
        (wd/f'Species{i}.fa').write_text(''.join(f'>{i}_{g}\n'+ 'ACDEFGHIKLMNPQRSTVWY'*3+'\n' for g in range(used[i])))
    clusters = wd/'clusters_OrthoFinder.txt_id_pairs.txt'
    clusters.write_text('(mclmatrixbegin\n'+''.join(f'{i} '+ ' '.join(genes)+' $\n' for i,genes in enumerate(families.values()))+')\n')
    (trees/'OG0000000.txt').write_text('((0_0:0.1,1_0:0.1)1:0.1,(2_0:0.1,3_0:0.1)1:0.1);\n')
    (wd/'SpeciesTree_unrooted_ids.txt').write_text('((0:0.1,1:0.1):0.1,(2:0.1,3:0.1):0.1);\n')
    (source/'Log.txt').write_text('Synthetic small-family output contract.\n'+f'WorkingDirectory_Base: {wd}/\nFN_Orthogroups: {clusters}\nWorkingDirectory_Trees: {wd}/\n')
    user_tree=out/'species_tree.nwk'
    user_tree.write_text('((Taxon0:0.1,Taxon1:0.1):0.1,(Taxon2:0.1,Taxon3:0.1):0.1);\n')
    before=snapshot(source)
    exe=ROOT/'.cache/envs/orthofinder/bin/orthofinder'
    command=[str(exe),'--from-trees',str(source),'-s',str(user_tree),'-n','small_fixture','-t','2','-a','1','-M','msa','-S','diamond','-A','famsa','-T','fasttree','--no-fix-files','--save-space']
    started=time.time()
    with (out/'stdout.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=out,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
        try:
            rc=proc.wait(timeout=120)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            raise RuntimeError('Synthetic fixture exceeded 120 seconds')
    if rc or snapshot(source)!=before:
        raise RuntimeError('Native execution failed or changed staged source')
    result=out/'Results_small_fixture'
    pairs=[]
    for i in range(4):
        with (result/f'Orthologues/Taxon{i}.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                for left in row[f'Taxon{i}'].split(', '):
                    for right in row['Orthologs'].split(', '):
                        pairs.append((row['Orthogroup'],left,right))
    expected={(og,'protein_'+left,'protein_'+right) for og,genes in families.items()
              for left in genes for right in genes if left.split('_')[0]!=right.split('_')[0]}
    actual=set(pairs)
    if len(actual)!=len(pairs) or actual-expected:
        raise ValueError('Duplicate or unexpected directed pair')
    missing=sorted(expected-actual)
    package=ROOT/'.cache/envs/orthofinder/lib/python3.12/site-packages/orthofinder'
    receipt=dict(status='completed_synthetic_small_family_output_observation',genes=sum(used),families=families,early_singleton=args.early_singleton,
        expected_directed_pairs=len(expected),observed_directed_pairs=len(actual),missing_directed_pairs=missing,
        source_unchanged=True,source_hashes=before,command=command,returncode=rc,
        script_sha256=sha(__file__),executable_sha256=sha(exe),
        native_orthologues_source_sha256=sha(package/'comparative_genomics/orthologues.py'),
        result_hashes=snapshot(result),elapsed_seconds=time.time()-started,
        scope='Synthetic output coverage observation, not a passing full reconciliation audit or a biological validation. Missing expected pairs remain explicit.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','genes','expected_directed_pairs','observed_directed_pairs','missing_directed_pairs']},indent=2))


if __name__=='__main__':
    main()
