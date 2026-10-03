"""Full observed-marker coverage on all closed candidate species-tree views.

Pruning diagnoses observation geometry, never statistical estimability or
accepted rooting. Original branch units remain source specific.
"""
from collections import defaultdict
import json
import math
from pathlib import Path

import numpy as np

from audit_selected_taxon_identity_snapshot_v2 import rows, sha
from coalescent_tree_comparison_sources import load as load_tree_sources

SOURCES = {
    'AlphaFold': ('results/phylogeny/paired-inputs-afdb-recovered-20260925-v1',
                 'metadata/recovered_afdb_paired_inputs_readback_20260925.json'),
    'ESMFold': ('results/phylogeny/paired-inputs-esmfold-all-completed-20260922-v1',
                'metadata/esmfold_all_completed_paired_inputs_completed_readback.json'),
}
STATUSES = ['insufficient_family_taxa','no_observations_on_one_side',
            'terminal_projection','unique_internal_projection','shared_internal_projection']
DTYPES = {'status':'uint8','observed_side':'uint16','observed_complement':'uint16',
          'path_internal_branch_count':'uint16','path_length_sum':'float64'}


def fasta(path):
    result, label, sequence = {}, None, []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if label is not None:
                    assert label not in result
                    result[label] = ''.join(sequence)
                label, sequence = line[1:].strip(), []
                assert label and not any(c.isspace() for c in label)
            else:
                assert label is not None
                sequence.append(line.strip())
    if label is not None:
        assert label not in result
        result[label] = ''.join(sequence)
    return result


def normalize(mask, full):
    other = full ^ mask
    assert mask & ~full == 0
    return min([mask,other], key=lambda value:(value.bit_count(),value))


def mask_for(taxa, positions):
    return sum(1 << positions[t] for t in taxa)


def project(branches, retained):
    """All original internal branches; terminal/lost paths have no length claim."""
    n = retained.bit_count()
    arrays = {key:np.zeros(len(branches),dtype=dtype) for key,dtype in DTYPES.items()}
    arrays['path_length_sum'].fill(np.nan)
    groups = defaultdict(list)
    for i,branch in enumerate(branches):
        side = int(branch['split_mask_hex'],16) & retained
        a,b = side.bit_count(),n-side.bit_count()
        arrays['observed_side'][i], arrays['observed_complement'][i] = a,b
        if n < 4: continue
        if min(a,b) == 0: arrays['status'][i] = 1
        elif min(a,b) == 1: arrays['status'][i] = 2
        else: groups[normalize(side,retained)].append(i)
    for indices in groups.values():
        length = math.fsum(branches[i]['branch_length'] for i in indices)
        for i in indices:
            arrays['status'][i] = 3 if len(indices)==1 else 4
            arrays['path_internal_branch_count'][i] = len(indices)
            arrays['path_length_sum'][i] = length
    if n >= 4:
        assert len(groups) == n-3
    return arrays,len(groups)


def load():
    bindings = {}
    def bind(path, expected=None):
        path = str(path); digest = sha(path)
        assert expected is None or expected == digest,path
        assert path not in bindings or bindings[path] == digest,path
        bindings[path] = digest
    tree_plan_path = Path('metadata/coalescent_tree_comparison_plan_20261002_v2.json')
    tree_plan = json.loads(tree_plan_path.read_text())
    completion_path = Path('metadata/coalescent_tree_comparisons_completed_20261002_v2.json')
    completion = json.loads(completion_path.read_text())
    assert completion['status']=='complete_verified_full_coalescent_reference_tree_comparisons'
    assert (completion['cohorts'],completion['tree_views']) == (5,70)
    bind(completion_path); bind(completion['full_hash_archive'],completion['full_hash_archive_sha256'])
    archive=json.loads(Path(completion['full_hash_archive']).read_text())
    assert len(archive['services'])==completion['exact_process_journals_checked']==2
    for service in archive['services']:
        assert service['captured_process_messages']>0 and service['completion_resource_records']
        assert service['observed_state']['Result']=='success' and service['observed_state']['ExecMainStatus']=='0'
    for path,digest in archive['source_hashes'].items(): bind(path,digest)
    groups,universe,prior=load_tree_sources(tree_plan,tree_plan_path)
    for path,digest in prior.items(): bind(path,digest)
    manifest_path=Path('metadata/analysis_manifest.tsv');bind(manifest_path)
    manifest={r['taxon_id']:r for r in rows(manifest_path)}
    assert universe==set(manifest) and len(universe)==526
    positions={t:i for i,t in enumerate(sorted(universe))}
    overlay_path=Path('config/assembly_identity_reviews_20261003_v2.json');bind(overlay_path)
    overlay=json.loads(overlay_path.read_text())['records']
    for taxon,record in overlay.items():
        assert record['assembly_accession']==manifest[taxon]['assembly_accession']
        bind(record['evidence_path'],record['evidence_sha256'])
    identity_path=Path('results/taxonomy/selected-identity-snapshot-20261003-v2/taxon_identity.tsv')
    identity_receipt=Path('metadata/selected_taxon_identity_snapshot_completed_20261003_v2.json')
    bind(identity_receipt);identity_proof=json.loads(identity_receipt.read_text())
    bind(identity_path,identity_proof['artifacts'][str(identity_path)])
    identities={r['taxon_id']:r for r in rows(identity_path)}
    assert set(identities)==universe
    coverage_path=Path('metadata/taxon_matrix_coverage.tsv')
    coverage_receipt=Path('metadata/taxon_matrix_coverage_sensitivity_receipt.json')
    bind(coverage_receipt)
    bind(coverage_path,json.loads(coverage_receipt.read_text())['artifacts']['taxon_matrix_coverage.tsv'])
    matrix_coverage={r['taxon_id']:r for r in rows(coverage_path)}
    assert set(matrix_coverage)==universe
    datasets,marker_ids={},None
    for source,(directory,readback_path) in SOURCES.items():
        root=Path(directory);receipt=root/'receipt.json';bind(receipt)
        record=json.loads(receipt.read_text());readback=json.loads(Path(readback_path).read_text());bind(readback_path)
        assert readback['status']=='passed_complete_paired_inputs_from_qualified_arrays_readback'
        assert readback['source_receipt_sha256']==sha(receipt)
        assert (record['taxa_audited'],record['markers_audited'],record['taxon_marker_rows'])==(526,125,65750)
        for name,digest in record['artifacts'].items():bind(root/name,digest)
        for path,digest in readback['input_receipts'].items():bind(path,digest)
        summaries={r['marker']:r for r in rows(root/'marker_summary.tsv')}
        if marker_ids is None: marker_ids=sorted(summaries)
        assert sorted(summaries)==marker_ids and len(marker_ids)==125
        cell_rows=rows(root/'taxon_coverage.tsv')
        cells={(r['marker'],r['taxon_id']):r for r in cell_rows}
        assert len(cells)==len(cell_rows)==65750
        assert set(cells)=={(m,t) for m in marker_ids for t in universe}
        eligible={m:set() for m in marker_ids}
        categories=['observed','noncanonical_or_missing_sequence','no_structural_mapping',
                    'invalid_native_feature','low_feature_plddt','high_feature_pae']
        for (m,t),row in cells.items():
            assert sum(int(row[k]) for k in categories)==int(row['original_marker_columns'])
            required=max(50,math.ceil(.3*int(row['original_marker_columns'])))
            assert required==int(row['required_observed_columns'])
            assert row['taxon_eligible'] in ['True','False']
            assert (int(row['observed'])>=required)==(row['taxon_eligible']=='True')
            if row['taxon_eligible']=='True':eligible[m].add(t)
        usable={}
        for marker in marker_ids:
            summary=summaries[marker]
            ready=summary['status']=='ready_for_inference'
            assert int(summary['eligible_taxa'])==len(eligible[marker])
            assert ready==(len(eligible[marker])>=4)
            if ready:
                aa=fasta(root/marker/'aa.faa');states=fasta(root/marker/'3di.faa')
                assert set(aa)==set(states)==eligible[marker]
                assert all(len(v)==int(summary['retained_columns']) for v in [*aa.values(),*states.values()])
                for taxon in aa:
                    assert all((a=='?')==(b=='?') for a,b in zip(aa[taxon],states[taxon]))
                    assert sum(a!='?' for a in aa[taxon])==int(cells[marker,taxon]['observed'])
            usable[marker]=eligible[marker] if ready else set()
        assert len(set().union(*usable.values()))==readback['taxa']
        assert sum(map(len,usable.values()))==readback['marker_taxon_cells']
        datasets[source]=dict(eligible=eligible,usable=usable,receipt=str(receipt),
            mask=record['mask'],eligibility=record['eligibility'])
    assert len({v['mask'] for v in datasets.values()})==len({v['eligibility'] for v in datasets.values()})==1
    views=[]
    for cohort,group in sorted(groups.items()):
        assert set(group['taxa'])<=universe and len(group['views'])==14
        full=mask_for(group['taxa'],positions)
        for name,view in sorted(group['views'].items()):
            branches=[]
            for split,metric in sorted(view['splits'].items()):
                branch_mask=normalize(mask_for(split,positions),full)
                length=metric['branch_length']
                assert math.isfinite(length) and length>=0
                branches.append(dict(split_mask_hex=hex(branch_mask),branch_length=length,
                    branch_length_unit=metric['branch_length_unit']))
            assert len(branches)==len(group['taxa'])-3
            assert len({r['split_mask_hex'] for r in branches})==len(branches)
            views.append(dict(cohort=cohort,view=name,family=view['family'],
                source_tree=view['source_tree'],prune=view['prune'],taxa=sorted(group['taxa']),
                full_mask_hex=hex(full),roles=group['roles'],branches=branches))
    assert len(views)==70 and sum(len(v['branches']) for v in views)==36190
    taxa=[]
    for taxon,row in sorted(manifest.items()):
        contributions={source:sum(taxon in dataset['usable'][m] for m in marker_ids)
                       for source,dataset in datasets.items()}
        flags=identities[taxon]['identity_review_flags'].split(';')
        if taxon in overlay:flags.append(overlay[taxon]['review_class'])
        taxa.append(dict(taxon_id=taxon,frozen_species_name=row['species_name'],
            display_name=overlay.get(taxon,{}).get('display_name',row['species_name']),
            study_role=row['study_role'],manifest_lineage_group=row['lineage'].split(';')[0],
            identity_review_flags=';'.join(flag for flag in flags if flag),
            profile_canonical_fraction=matrix_coverage[taxon]['profile_fraction'],
            mafft_canonical_fraction=matrix_coverage[taxon]['mafft_fraction'],
            AlphaFold_usable_markers=contributions['AlphaFold'],ESMFold_usable_markers=contributions['ESMFold'],
            same_marker_predictor_overlap=sum(taxon in datasets['AlphaFold']['usable'][m] and
                taxon in datasets['ESMFold']['usable'][m] for m in marker_ids),
            either_predictor_markers=sum(any(taxon in datasets[s]['usable'][m] for s in SOURCES) for m in marker_ids)))
    return views,positions,marker_ids,datasets,taxa,bindings
