"""Full-grid joint-logger job construction and native output qualification."""
from collections import Counter, defaultdict
import copy
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_node_logger_v5 import restore
from baliphy_reference_sampler_qualification import inspect as legacy_inspect, summarize as legacy_summary
from baliphy_reference_sampler_qualification import SUMMARY_FIELDS as LEGACY_FIELDS
from independent_joint_ancestral_frames import decode, verify_arrays, write_arrays
from independent_native_ancestral_alignment import fasta_records
from independent_native_ancestral_topology import match_trees
from independent_short_sampler_outputs_v2 import mapping_rows, strict_json


SUCCESS = 'full_joint_short_sampler_output_integrity_checked_not_posterior'
INVALID = 'invalid_joint_sampler_qualification_output_retained'
SUMMARY_FIELDS = LEGACY_FIELDS + ['joint_saved_frames', 'joint_ancestral_records',
    'joint_ancestral_residue_category_pairs', 'joint_tip_residue_category_pairs',
    'joint_anchored_candidate_cells', 'joint_unanchored_candidate_residues',
    'joint_category_data_available_for_checked_roles']


def build_jobs(roles, originals):
    assert len(roles) == len(originals) == 1620
    by_id = {j['chain']['chain_id']:j for j in originals}; assert len(by_id) == 1620
    assert {r['source_chain_id'] for r in roles} == set(by_id)
    assert len({r['chain']['chain_id'] for r in roles}) == 1620
    seeds = {r['chain']['seed'] for r in roles}
    forbidden = {j['chain']['seed'] for j in originals} | {j['source_seed'] for j in originals} | {20267001,20267002,20267003}
    assert len(seeds) == 1620 and seeds.isdisjoint(forbidden)
    assert all(type(seed) is int and 1 <= seed < 2**31 for seed in seeds)
    jobs = []; groups = defaultdict(list)
    for role in roles:
        chain = role['chain']; source = by_id[role['source_chain_id']]; old = source['chain']
        assert role['native_execution_launched'] is role['scientific_eligibility'] is role['posterior_qualified'] is False
        assert role['current_short_sampler_seed'] == old['seed']
        for key in ['effective_input_group','prior_label','chain','original_configuration_ids','family',
                    'proteins','alignment','alignment_sha256','tree','tree_sha256']:
            assert chain[key] == old[key], key
        for c in [old,chain]: assert sha(c['program']) == c['program_sha256']
        assert restore(Path(chain['program']).read_text()) == Path(old['program']).read_text()
        config = copy.deepcopy(source['config']); command = config['command']
        assert command[command.index('--seed')+1] == str(old['seed'])
        assert command[command.index('run')+1] == old['program']
        assert command[command.index('--iterations')+1] == '20' and '--test' not in command
        assert source['memory_reservation_bytes'] in [12*2**30,48*2**30]
        assert ('--as='+str(source['memory_reservation_bytes'])) in command
        assert 900 <= config['timeout_seconds'] <= 7200
        command[command.index('--seed')+1] = str(chain['seed'])
        command[command.index('run')+1] = chain['program']
        assert config['pins'].pop(old['program']) == old['program_sha256']
        config['pins'][chain['program']] = chain['program_sha256']
        api = Path(command[5]).parent.parent/'lib/bali-phy/haskell'
        for library in ['Data/JSON.hs','Data/JSON/Encoding.hs','Data/JSON/Types/Foreign.hs',
                        'Data/JSON/Types/ToJSON.hs','SModel/Property.hs']:
            path = api/library
            config['pins'][str(path)] = sha(path)
        for path,h in config['pins'].items(): assert sha(path) == h,path
        job = dict(source,chain=dict(chain),config=config,source_seed=old['seed'],
                   source_chain_id=old['chain_id'],initialization_only_seed_reused_for_first_mcmc=True)
        jobs.append(job); groups[chain['effective_input_group']+'-'+chain['prior_label']].append(chain)
    assert len(groups) == 405 and all(len(cs)==4 and {c['chain'] for c in cs}=={1,2,3,4} for cs in groups.values())
    assert Counter(j['chain']['prior_label'] for j in jobs) == {'broad':540,'centered':540,'package':540}
    assert len({j['chain']['effective_input_group'] for j in jobs}) == 135
    assert len({a for j in jobs for a in j['chain']['original_configuration_ids']}) == 324
    assert Counter(j['memory_reservation_bytes']//2**30 for j in jobs) == {12:1512,48:108}
    return jobs


def joint_frames(chain, directory, mapping):
    matched = match_trees(Path(chain['tree']).read_text(),(directory/'runtime-tree.nwk').read_text())
    observed = {label:seq.replace('-','') for label,seq in fasta_records(Path(chain['alignment']).read_text().splitlines()).items()}
    assert sorted(observed) == matched['tips'] and len(observed) == chain['proteins']
    bits = {tip:1<<i for i,tip in enumerate(matched['tips'])}; candidates = {}
    for row in mapping_rows(mapping,chain):
        descendants = strict_json(row['retained_set_json']); node = row['source_node']
        assert descendants and len(descendants) == len(set(descendants)) and set(descendants) <= set(bits)
        mask = sum(bits[tip] for tip in descendants)
        assert node not in candidates and matched['source_labels'][node]['mask'] == mask
        assert (int(row['level'])==3) == (mask==(1<<len(bits))-1)
        if 'retained_descendants' in row: assert int(row['retained_descendants']) == len(descendants)
        candidates[node] = matched['runtime_index'][mask]['label']
    assert len(candidates) == len(set(candidates.values())) == 4
    decoded = []
    with (directory/'C1.P1.site-property-samples.jsonl').open() as handle:
        for index,line in enumerate(handle):
            assert index < 3,'Extra joint sample'
            frame = strict_json(line)
            decoded.append(decode(frame,[0,10,20][index],observed,matched['runtime_labels'],candidates))
    assert len(decoded) == 3,'Missing joint sample'
    return decoded,matched,candidates


def inspect(job, receipt_path, plan_digest, mapping, export_root, allow_export_creation=True):
    """Retain native/decoder failures; export failures abort the controller."""
    export_root = Path(export_root)
    row = legacy_inspect(job,receipt_path,plan_digest,mapping)
    row.update(legacy_output_status=row['status'],legacy_saved_alignments=row['saved_alignments'],
        legacy_candidate_frames=row['candidate_frames'],joint_frames=[],
        ancestral_categories_available=False,same_record_sequence_category_correspondence=False)
    if row['status'] != 'full_short_sampler_output_integrity_checked_not_posterior':
        assert not export_root.exists(),'Unresolved role cannot acquire arrays'
        return row
    try:
        directories = list(Path(receipt_path).parent.glob('independent-chain-*')); assert len(directories) == 1
        decoded,matched,candidates = joint_frames(job['chain'],directories[0],mapping)
    except (AssertionError,KeyError,ValueError,ArithmeticError,OSError,TypeError) as error:
        assert not export_root.exists(),'Invalid role cannot acquire arrays'
        row.update(status=INVALID,saved_alignments=0,candidate_frames=0,
                   joint_error_type=type(error).__name__,joint_error=str(error))
        return row
    # Only after all native frames pass can a role obtain successful exports.
    # Existing exports are verified byte/value/dtype/shape-wise, never replaced.
    assert type(allow_export_creation) is bool
    exists = export_root.exists()
    assert exists or allow_export_creation, "Readback must never recreate missing exports"
    if not exists: export_root.mkdir(parents=True,exist_ok=False)
    expected = set()
    for summary,arrays in decoded:
        path = export_root/('frame-'+str(summary['iteration'])+'.npz'); expected.add(path)
        if not exists: write_arrays(path,arrays)
        verify_arrays(path,arrays)
        row['joint_frames'].append(dict(summary,projection_array=str(path),projection_array_sha256=sha(path)))
    assert set(export_root.iterdir()) == expected,'Foreign/missing joint-frame artifact'
    row.update(status=SUCCESS,ancestral_categories_available=True,
        same_record_sequence_category_correspondence=True,joint_candidate_mapping=candidates,
        joint_node_rows=matched['rows'],negative_branch_pairs_require_review=matched['negative_branch_pairs_require_review'])
    return row


def summarize(rows):
    assert len(rows) == len({r['chain_id'] for r in rows}) == 1620
    compatible = []
    for row in rows:
        assert row['scientific_eligibility'] is row['posterior_qualified'] is False
        if row['status'] == SUCCESS:
            assert row['ancestral_categories_available'] is row['same_record_sequence_category_correspondence'] is True
            assert [f['iteration'] for f in row['joint_frames']] == [0,10,20]
            assert all(f['posterior_qualified'] is f['scientific_eligibility'] is False for f in row['joint_frames'])
            compatible.append(dict(row,status='full_short_sampler_output_integrity_checked_not_posterior'))
        else:
            assert row['status'] in ['unsuccessful_sampler_qualification_attempt_retained','invalid_sampler_qualification_output_retained',INVALID]
            assert row['joint_frames'] == [] and row['ancestral_categories_available'] is False
            compatible.append(dict(row,status='invalid_sampler_qualification_output_retained'))
    base = legacy_summary(compatible); base['status_counts'] = dict(Counter(r['status'] for r in rows))
    frames = [frame for row in rows for frame in row['joint_frames']]
    return dict(base,joint_saved_frames=len(frames),
        joint_ancestral_records=sum(f['native_ancestors'] for f in frames),
        joint_ancestral_residue_category_pairs=sum(f['ancestral_pairs'] for f in frames),
        joint_tip_residue_category_pairs=sum(f['tip_pairs'] for f in frames),
        joint_anchored_candidate_cells=sum(f['anchored_observations'] for f in frames),
        joint_unanchored_candidate_residues=sum(f['unanchored_candidate_residues'] for f in frames),
        joint_category_data_available_for_checked_roles=bool(frames))
