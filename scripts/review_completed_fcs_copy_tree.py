#!/usr/bin/env python3
"""Review candidate and selected-copy neighborhoods in a completed unrooted tree."""
import argparse
import csv
import json
import math
from pathlib import Path
from Bio import Phylo, SeqIO
from audit_busco_gene_copies import sha


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['trees', 'inputs', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    p.add_argument('--marker', required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    folder = a.trees / a.marker
    r = json.loads((folder / 'receipt.json').read_text())
    if r['status'] != 'inferred' or r['marker'] != a.marker:
        raise ValueError('Completed marker receipt required')
    for name, digest in r['artifacts'].items():
        if sha(folder / name) != digest:
            raise ValueError('Tree artifact changed')
    config = json.loads((folder / 'run_config.json').read_text())
    common = json.loads((a.trees / 'run_config.json').read_text())
    if config['common_config_sha256'] != sha(a.trees / 'run_config.json') or common['input_receipt_sha256'] != sha(a.inputs / 'receipt.json'):
        raise ValueError('Input configuration lineage differs')
    alignment = a.inputs / a.marker / 'input.faa'
    if sha(alignment) != config['input_sha256']:
        raise ValueError('Alignment changed')
    command = config['command']
    if command[command.index('--alrt')+1] != '1000' or '-B' in command:
        raise ValueError('Unexpected support method')
    records = list(SeqIO.parse(alignment, 'fasta'))
    tree = Phylo.read(folder / 'tree.treefile', 'newick')
    tips = {x.name: x for x in tree.get_terminals()}
    if len(tips) != len(records) or set(tips) != {x.id for x in records} or len(tips) != r['tips']:
        raise ValueError('Tip grid differs')
    if any(len(x.seq) != config['columns'] for x in records):
        raise ValueError('Alignment dimensions differ')
    taxa_path = Path('metadata/analysis_manifest.tsv')
    taxa = {row['taxon_id']: row for row in table(taxa_path)}
    candidates = sorted(set(tips)-set(taxa))
    if not candidates or any(len(x.split('|')) != 3 or x.split('|')[1] != a.marker or x.split('|')[0] not in taxa for x in candidates):
        raise ValueError('Invalid candidate identities')
    focal = sorted(set(candidates) | {x.split('|')[0] for x in candidates if x.split('|')[0] in tips})
    graph = {node: [] for node in tree.find_clades()}
    for node in graph:
        for child in node.clades:
            length = child.branch_length
            if length is None or not math.isfinite(length) or length < 0:
                raise ValueError('Invalid branch length')
            if child.confidence is not None and not 0 <= child.confidence <= 100:
                raise ValueError('Invalid support')
            graph[node].append((child, length)); graph[child].append((node, length))
    def info(tip):
        taxon = tip.split('|')[0]
        return {'tip': tip, 'taxon_id': taxon, 'species_name': taxa[taxon]['species_name'],
                'lineage': taxa[taxon]['lineage'], 'kind': 'candidate' if tip in candidates else 'selected'}
    neighborhoods, checked, maximum = [], 0, 0.0
    for name in focal:
        start = tips[name]; stack = [(start, None, 0.0)]; distances = {}
        while stack:
            node, parent, distance = stack.pop()
            if node.is_terminal():
                distances[node.name] = distance
            stack.extend((n, node, distance+length) for n,length in graph[node] if n is not parent)
        for other, distance in distances.items():
            error = abs(distance-tree.distance(name,other)); maximum = max(maximum,error); checked += 1
            if error > 1e-10:
                raise ValueError('Independent graph path differs')
        nearest = sorted((d,t) for t,d in distances.items() if t != name)[:10]
        sides = []
        for node in tree.find_clades():
            if node is tree.root:
                continue
            side = {tip.name for tip in node.get_terminals()}
            group = side if name in side else set(tips)-side
            if 2 <= len(group) <= 8:
                sides.append({'tips': sorted(group), 'reported_sh_alrt': node.confidence,
                              'separating_branch_length': node.branch_length})
        neighborhoods.append({'focal': info(name), 'nearest_tree_path_tips': [dict(info(t),path_length=d) for d,t in nearest],
                              'all_unrooted_split_sides_of_size_2_to_8': sorted(sides,key=lambda x:(len(x['tips']),x['tips']))})
    fcs_path = Path('results/qc/fcs-cds-overlap-v1/marker_overlap_review.tsv')
    fr = json.loads((fcs_path.parent/'receipt.json').read_text())
    if sha(fcs_path) != fr['artifacts'][fcs_path.name]:
        raise ValueError('FCS mapping changed')
    flags = [row for row in table(fcs_path) if row['marker'] == a.marker and row['taxon_id'] in focal]
    a.output.mkdir(parents=True)
    result = {'status': 'complete_exploratory_candidate_copy_tree_neighborhood_review',
              'marker': a.marker, 'tips': len(tips), 'columns': config['columns'],
              'candidate_tips': candidates, 'focal_neighborhoods': neighborhoods, 'selected_copy_fcs_rows': flags,
              'independent_graph_path_checks': checked, 'maximum_path_difference': maximum,
              'source_receipt_sha256': sha(folder/'receipt.json'), 'input_receipt_sha256': sha(a.inputs/'receipt.json'),
              'tree_sha256': sha(folder/'tree.treefile'), 'manifest_sha256': sha(taxa_path),
              'fcs_mapping_receipt_sha256': sha(fcs_path.parent/'receipt.json'), 'script_sha256': sha(Path(__file__)),
              'interpretation': 'Exploratory unrooted neighborhoods for all candidate tips and available same-taxon selected tips. All focal-to-tip paths independently checked by graph traversal. Small split sides are exhaustively reported for sizes 2–8, not asserted rooted clades. Reported support is SH-aLRT, not bootstrap probability. No candidate is accepted as an ortholog, source contamination is not proven, and full paralog/root/reconciliation and wider copy-tree review remain necessary.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'marker':a.marker,'tips':len(tips),'candidate_tips':candidates,'independent_graph_path_checks':checked,'maximum_path_difference':maximum},indent=2))


if __name__ == '__main__':
    main()
