#!/usr/bin/env python3
"""Describe ten retained per-fit probability extremes; no causal attribution."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np
from Bio import SeqIO

ROOT = Path('results/ancestral')


def main():
    hashes = {}

    def bind(path):
        path = Path(path)
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    def read_json(path):
        return json.loads(bind(path).read_text())

    def rows(path):
        with bind(path).open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    completion = read_json('metadata/whole_optimization_probability_comparison_full_completed_20260927.json')
    source = bind(completion['completed_receipt_path'])
    assert hashes[str(source)] == completion['completed_receipt_sha256']
    plan = read_json('metadata/ancestral_whole_multistart_plan_20260927.json')
    jobs = {j['job_id']: j for j in plan['jobs']}
    mapping = rows(ROOT/'case-local-trees-20260927-v1/ancestral_node_mapping.tsv')
    fits = {}
    for label, folder in [('baseline', 'whole-protein-model-readback-20260927-v1'),
                          ('alternate', 'whole-multistart-readback-20260927-v1'),
                          ('refined', 'whole-refinement-readback-20260927-v1')]:
        fits[label] = {r['job_id']: r for r in rows(ROOT/folder/'fit_readback.tsv')}

    folders = {'baseline': 'whole-protein-ancestors-20260927-v1',
               'alternate': 'alternate-whole-ancestors-20260927-v1',
               'refined': 'refined-whole-ancestors-20260927-v2'}
    receipts = {k: read_json(ROOT/v/'receipt.json') for k, v in folders.items()}

    def vector(label, job, level, col):
        path = bind(ROOT/folders[label]/(job+'.npz'))
        assert hashes[str(path)] == receipts[label]['artifacts'][path.name]
        with np.load(path) as data:
            index = list(data['levels']).index(level)
            aa = data['amino_acids'].tolist()
            vec = data['posterior'][index, col].copy()
        assert np.isfinite(vec).all() and abs(float(vec.sum())-1) < 1e-8
        return aa, vec

    results = []
    for extreme in completion['largest_alternate_site_shifts']:
        job = jobs[extreme['alternate_job']]
        base = extreme['baseline_job']
        col = extreme['column']-1
        level = extreme['level']
        path = bind(job['alignment'])
        assert hashes[str(path)] == job['alignment_sha256']
        records = {r.id: str(r.seq) for r in SeqIO.parse(path, 'fasta')}
        assert len(records) == job['proteins']
        mappings = [r for r in mapping if r['family'] == job['family']
                    and r['dataset'] == 'whole' and int(r['level']) == level]
        sets = {tuple(sorted(json.loads(r['retained_set_json']))) for r in mappings}
        assert len(sets) == 1
        descendants = set(next(iter(sets)))
        assert descendants <= set(records)
        observations = {}
        for group, names in [('all', set(records)), ('descendants', descendants),
                             ('outside', set(records)-descendants)]:
            counts = Counter(records[name][col] for name in names)
            observations[group] = dict(tips=len(names), characters=dict(sorted(counts.items())),
                                       canonical_amino_acids=sum(n for a,n in counts.items() if a in 'ARNDCQEGHILKMFPSTWYV'))
        ids = {'baseline': base, 'alternate': extreme['alternate_job'],
               'refined': base+'-refine-amin'+str(job['alpha_min'])}
        values = {}; vectors = {}
        for label, jid in ids.items():
            aa, vec = vector(label, jid, level, col)
            vectors[label] = vec
            fit = fits[label][jid]
            values[label] = dict(job_id=jid, log_likelihood=float(fit['checkpoint_log_likelihood']),
                                gamma_shape=float(fit['gamma_shape']),
                                edges_below_1e_5=int(fit['edges_below_1e_5']),
                                probabilities=dict(zip(aa, vec.tolist())),
                                most_probable_amino_acid=aa[int(vec.argmax())])
        tv = float(np.abs(vectors['alternate']-vectors['baseline']).sum()/2)
        assert abs(tv-extreme['total_variation']) < 1e-12
        results.append(dict(**extreme, observations=observations, fits=values,
                            alternate_minus_baseline_log_likelihood=values['alternate']['log_likelihood']-values['baseline']['log_likelihood'],
                            refined_versus_baseline_total_variation=float(np.abs(vectors['refined']-vectors['baseline']).sum()/2)))
    output = ROOT/'whole-ancestral-extreme-diagnostics-20260927-v1'
    output.mkdir(exist_ok=False)
    artifact = output/'diagnostics.json'
    artifact.write_text(json.dumps(results, indent=2)+'\n')
    bind(__file__)
    receipt = dict(status='complete_descriptive_diagnostics', retained_per_fit_extremes=len(results),
                   source_hashes=hashes, artifacts={artifact.name: hashlib.sha256(artifact.read_bytes()).hexdigest()},
                   scope='Ten previously retained per-fit maxima, not the ten largest individual sites globally. Coverage describes observed alignment characters, not alignment correctness. Likelihood differences and parameter associations do not establish a causal mechanism. Probabilities are conditional on residue existence; no indel or biological validation is implied.')
    (output/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'output':str(output),'largest':results[0]}))


if __name__ == '__main__':
    main()
