"""V6 scalar admission composed with preserved joint short-sampler checks.

This adapter does not launch samplers or adopt a posterior. Full-grid native
startup, resource, serialization and closure gates still apply separately.
"""
from collections import Counter
import copy
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_qualification_v3 import (
    inspect as joint_inspect, summarize as joint_summary, SUCCESS, INVALID,
)
from baliphy_reference_sampler_qualification import inspect as legacy_inspect
from baliphy_scalar_json_logger_v6c import restore, validate_record, SCHEMA
from read_baliphy_scalar_json_v6b import load, read_record, compare_tsv

REVIEW = 'explicit_nonfinite_or_literal_null_scalar_output_retained_for_review'
FIELDS = ['effective_input_group','prior_label','chain','original_configuration_ids',
          'family','proteins','alignment','alignment_sha256','tree','tree_sha256']


def build_jobs(roles, previous_jobs, forbidden_seeds):
    assert len(roles) == len(previous_jobs) == 1620
    originals = {j['chain']['chain_id']:j for j in previous_jobs}
    assert len(originals) == 1620
    assert {r['source_v5_chain_id'] for r in roles} == set(originals)
    ids = {r['chain']['chain_id'] for r in roles}
    seeds = {r['chain']['seed'] for r in roles}
    assert len(ids) == len(seeds) == 1620 and seeds.isdisjoint(forbidden_seeds)
    assert all(type(seed) is int and 1 <= seed < 2**31 for seed in seeds)
    jobs = [];groups = {}
    for role in roles:
        chain = role['chain']; previous = originals[role['source_v5_chain_id']]; old = previous['chain']
        assert role['native_execution_launched'] is role['scientific_eligibility'] is role['posterior_qualified'] is False
        assert role['seed_namespace'] == 'fungal-scalar-cjson-explicit-nonfinite-future-20261004-v6'
        assert role['source_v5_seed'] == old['seed']
        assert chain['chain_id'] == chain['effective_input_group']+'-'+chain['prior_label']+'-scalar-cjson-v6-chain'+str(chain['chain'])
        assert all(chain[k] == old[k] for k in FIELDS)
        assert sha(chain['program']) == chain['program_sha256'] and sha(old['program']) == old['program_sha256']
        assert restore(Path(chain['program']).read_text()) == Path(old['program']).read_text()
        config = copy.deepcopy(previous['config']); command = config['command']
        assert command[command.index('--seed')+1] == str(old['seed'])
        assert command[command.index('run')+1] == old['program']
        assert command[command.index('--iterations')+1] == '20' and '--test' not in command
        assert previous['memory_reservation_bytes'] in [12*2**30,48*2**30]
        assert '--as='+str(previous['memory_reservation_bytes']) in command
        assert 900 <= config['timeout_seconds'] <= 7200
        command[command.index('--seed')+1] = str(chain['seed'])
        command[command.index('run')+1] = chain['program']
        assert config['pins'].pop(old['program']) == old['program_sha256']
        config['pins'][chain['program']] = chain['program_sha256']
        # The complete installed source tree is already globally bound by the
        # separately required V6 manifest; bind the newly imported modules here
        # without duplicating that entire archive in every native configuration.
        binary = Path(command[command.index('--')+1])
        api = binary.parent.parent/'lib/bali-phy/haskell'
        for name in ['Data/JSON/Types/Internal.hs','MCMC/Types.hs','Compiler/RealFloat.hs']:
            library = api/name
            config['pins'][str(library)] = sha(library)
        for path,h in config['pins'].items(): assert sha(path) == h,path
        job = dict(previous,chain=dict(chain),config=config,source_seed=old['seed'],
                   source_chain_id=old['chain_id'],project_scalar_schema=SCHEMA)
        jobs.append(job)
        group = chain['effective_input_group']+'-'+chain['prior_label']
        groups.setdefault(group,[]).append(chain['chain'])
    assert len(groups) == 405 and all(sorted(v) == [1,2,3,4] for v in groups.values())
    assert Counter(j['chain']['prior_label'] for j in jobs) == {'broad':540,'centered':540,'package':540}
    assert Counter(j['memory_reservation_bytes']//2**30 for j in jobs) == {12:1512,48:108}
    assert len({j['chain']['effective_input_group'] for j in jobs}) == 135
    assert len({x for j in jobs for x in j['chain']['original_configuration_ids']}) == 324
    return jobs


def scalar_audit(directory):
    directory = Path(directory); path = directory/'C1.log.json'
    values = [load(line) for line in path.read_text().splitlines()]
    assert values[0] == {'fields':['iter','prior','likelihood','posterior'],'nested':True,
                         'format':'MCON','version':'0.2','projectScalarSchema':SCHEMA}
    rows = values[1:];assert len(rows) == 21
    reviews = [];literal_rows = []
    for i,row in enumerate(rows):
        primary = validate_record(row); independent = read_record(row)
        assert primary['iteration'] == independent['iteration'] == i
        assert primary['finite_record'] == independent['finite_record']
        assert primary['context']['numeric_leaves'] + primary['parameters']['numeric_leaves'] == independent['numeric_leaves']
        reviews += [dict(iteration=i,**r) for r in independent['nonfinite_reviews']]
        if primary['context']['literal_null_paths'] or primary['parameters']['literal_null_paths']:literal_rows.append(i)
    mapped = None
    if not reviews and not literal_rows: mapped = compare_tsv(directory)
    return dict(schema=SCHEMA,rows=len(rows),nonfinite_reviews=reviews,literal_null_iterations=literal_rows,
                mapped_values_compared=0 if mapped is None else mapped['mapped_values_compared'],
                scalar_json_sha256=sha(path),scalar_tsv_sha256=sha(directory/'C1.log'),
                fully_finite_mapped_integrity_checked=mapped is not None,scientific_eligibility=False)


def inspect(job, receipt, plan_digest, mapping, export_root, allow_export_creation=True):
    # Native configuration/process/artifact custody is invariant-fatal. It is
    # checked before malformed scalar data can become a retained disposition.
    base = legacy_inspect(job,receipt,plan_digest,mapping)
    if base['status'] != 'full_short_sampler_output_integrity_checked_not_posterior':
        row = joint_inspect(job,receipt,plan_digest,mapping,export_root,allow_export_creation)
        return dict(row,scalar_v6_audit=None,scalar_v6_error=None,scalar_integrity_accepted=False)
    try:
        directories = list(Path(receipt).parent.glob('independent-chain-*'));assert len(directories) == 1
        audit = scalar_audit(directories[0])
    except (AssertionError,KeyError,ValueError,ArithmeticError,OSError,TypeError) as error:
        assert not Path(export_root).exists(),'Invalid scalar role cannot acquire or retain qualified arrays'
        return dict(base,status=INVALID,saved_alignments=0,candidate_frames=0,
            legacy_output_status=base['status'],legacy_saved_alignments=base['saved_alignments'],
            legacy_candidate_frames=base['candidate_frames'],joint_frames=[],ancestral_categories_available=False,
            same_record_sequence_category_correspondence=False,scalar_v6_audit=None,
            scalar_v6_error=dict(error_type=type(error).__name__,error=str(error)),scalar_integrity_accepted=False)
    if not audit['fully_finite_mapped_integrity_checked']:
        assert not Path(export_root).exists(),'Scalar review must not silently retain qualified arrays'
        return dict(base,status=REVIEW,saved_alignments=0,candidate_frames=0,
            legacy_output_status=base['status'],legacy_saved_alignments=base['saved_alignments'],
            legacy_candidate_frames=base['candidate_frames'],joint_frames=[],ancestral_categories_available=False,
            same_record_sequence_category_correspondence=False,scalar_v6_audit=audit,
            scalar_v6_error=None,scalar_integrity_accepted=False)
    row = joint_inspect(job,receipt,plan_digest,mapping,export_root,allow_export_creation)
    return dict(row,scalar_v6_audit=audit,scalar_v6_error=None,scalar_integrity_accepted=row['status']==SUCCESS)


def summarize(rows):
    assert len(rows) == len({r['chain_id'] for r in rows}) == 1620
    compatible = []
    for row in rows:
        if row['status'] == SUCCESS:
            assert row['scalar_integrity_accepted'] and row['scalar_v6_audit']['fully_finite_mapped_integrity_checked']
        if row['status'] == REVIEW:
            assert not row['scalar_integrity_accepted'] and not row['joint_frames']
            compatible.append(dict(row,status=INVALID))
        else: compatible.append(row)
    result = joint_summary(compatible)
    result['status_counts'] = dict(Counter(r['status'] for r in rows))
    result['scalar_review_roles'] = sum(r['status']==REVIEW for r in rows)
    result['scalar_v6_finite_checked_roles'] = sum(r['scalar_integrity_accepted'] for r in rows)
    result['scalar_v6_mapped_values_checked'] = sum((r['scalar_v6_audit'] or {}).get('mapped_values_compared',0) for r in rows)
    return result
