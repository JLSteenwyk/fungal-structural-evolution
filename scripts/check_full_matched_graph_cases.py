#!/usr/bin/env python3
"""Regression cases for separated graph/covariate and matching proof contracts.

Proofs are explicitly synthetic. This checks refusal of wrong provenance, not
production journals, physical structures or scientific model validity.
"""
import argparse
import json
import tempfile
from pathlib import Path
from check_full_matched_coverage_cases import setup, write
from full_matched_graph_sources import load_graph
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run():
    labels = ['wrong_matching_closure_used_as_graph_proof', 'missing_graph_receipt_binding',
              'missing_original_journal', 'wrong_covariate_graph_lineage',
              'wrong_covariate_readback_lineage', 'wrong_fixed_matching_covariate_lineage',
              'wrong_full_graph_count', 'changed_node_bytes_without_updated_proof']
    with tempfile.TemporaryDirectory(prefix='full-matched-graph-contract-') as temp:
        root = Path(temp); plan = json.loads(setup(root).read_text())
        matching = json.loads((Path(plan['selection']) / 'receipt.json').read_text())
        bindings = {}; load_graph(plan, bindings, matching); verify(bindings)
        assert str(Path(plan['graph']) / 'receipt.json') not in json.loads(Path(plan['matching_completion']).read_text())['source_hashes']
        paths = [p for p in root.rglob('*') if p.is_file()]; saved = {p: p.read_bytes() for p in paths}
        for label in labels:
            for p, raw in saved.items(): p.write_bytes(raw)
            config = dict(plan); mr = dict(matching); gc = Path(plan['graph_completion']); closed = json.loads(gc.read_text())
            if label == labels[0]: config['graph_completion'] = plan['matching_completion']
            elif label == labels[1]: del closed['source_hashes'][str(Path(plan['graph']) / 'receipt.json')]; write(gc, closed)
            elif label == labels[2]: closed['services'].pop(); write(gc, closed)
            elif label == labels[3]:
                cp = Path(plan['covariates']) / 'receipt.json'; cr = json.loads(cp.read_text()); cr['graph_receipt_sha256'] = '0'*64; write(cp, cr)
                closed['source_hashes'][str(cp)] = sha(cp); write(gc, closed)
            elif label == labels[4]:
                cp = Path(plan['covariate_readback']); cr = json.loads(cp.read_text()); cr['producer_receipt_sha256'] = '0'*64; write(cp, cr)
                closed['source_hashes'][str(cp)] = sha(cp); write(gc, closed)
            elif label == labels[5]: mr['source_covariate_receipt_sha256'] = '0'*64
            elif label == labels[6]: closed['summary']['duplicate_target_links'] += 1; write(gc, closed)
            else:
                p = Path(plan['graph']) / 'target_nodes.jsonl'; p.write_text(p.read_text() + '\n')
            try:
                bindings = {}; load_graph(config, bindings, mr); verify(bindings)
            except (AssertionError, KeyError): pass
            else: raise AssertionError('False graph proof accepted: ' + label)
    return dict(status='passed_separate_full_matched_graph_source_contract_cases',
                matching_closure_has_no_graph_receipt=True, rejected_false_proofs=labels,
                synthetic_prior_proofs=True, production_journals_tested=False, scope=__doc__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    result = run(); result['script_sha256'] = sha(__file__)
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
