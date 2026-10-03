#!/usr/bin/env python3
"""Paired same-seed native logger checks across all three priors and full source grid."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import run_attempt, sha, write_json
from baliphy_joint_node_logger_v3 import BEFORE, AFTER, fresh_seeds, transform, restore, validate_frame
from independent_native_ancestral_alignment import fasta_records, native_alignments, tree_labels
from independent_short_sampler_outputs import strict_json
from run_baliphy_reference_preflight import verify


def rejected(action):
    try: action()
    except (AssertionError, KeyError, ValueError): return
    raise AssertionError('Invalid logger accepted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True); parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); root = args.output.resolve(); root.mkdir(exist_ok=False)
    native_path = Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    native = json.loads(native_path.read_text()); verify(native)
    jobs_path = Path(native['jobs']); jobs = json.loads(jobs_path.read_text()); assert len(jobs) == 1620
    programs = {j['chain']['program']: j['chain']['program_sha256'] for j in jobs}; assert len(programs) == 405
    bindings = {str(native_path): sha(native_path), str(jobs_path): sha(jobs_path), str(Path(__file__)): sha(__file__)}
    for source_path, expected in programs.items():
        path = Path(source_path); assert sha(path) == expected
        original = path.read_text(); changed = transform(original)
        assert restore(changed) == original
        bindings[source_path] = expected
    ids = [j['chain']['chain_id'] for j in jobs]
    used = {j['chain']['seed'] for j in jobs} | {j['source_seed'] for j in jobs}
    seeds = fresh_seeds(ids, used)
    assert seeds == fresh_seeds(list(reversed(ids)), used)
    assert len(seeds) == 1620 and set(seeds.values()).isdisjoint(used)
    rejected(lambda: fresh_seeds(ids[:-1]+[ids[0]], used))
    original = Path(jobs[0]['chain']['program']).read_text()
    rejected(lambda: transform(transform(original)))
    rejected(lambda: transform(original.replace(BEFORE, '')))
    rejected(lambda: transform(original+BEFORE))
    input_labels = set()
    for alignment_path in {j['chain']['alignment'] for j in jobs}:
        with Path(alignment_path).open() as handle:
            sequences = fasta_records(handle)
        assert all(not any(ord(c) < 32 or c in '\"\\' for c in label) for label in sequences)
        input_labels.update(sequences)
    alignment = root/'alignment.faa'; tree = root/'tree.nwk'
    alignment.write_text('>a\nA-CDE-FG\n>b\nA-C-EYFG\n>c\nATCDE-F-\n>d\n--CDEF-G\n>e\nAX-DE-FG\n')
    tree.write_text('((((a:0.1,b:0.1)n0:0.1,c:0.1)n1:0.1,d:0.1)n2:0.1,e:0.1)n3;\n')
    bindings[str(alignment)] = sha(alignment); bindings[str(tree)] = sha(tree)
    binary = Path(jobs[0]['config']['command'][5]); prlimit = Path('/usr/bin/prlimit')
    api = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/lib/bali-phy/haskell').resolve()
    apis = [api / x for x in ['Graph.hs','Tree.hs','Bio/Alignment.hs','Bio/Alphabet.hs',
             'Probability/Distribution/PhyloCTMC/Properties.hs','SModel/ASRV.hs','Probability/Logger.hs',
             'Data/Text.hs','Data/OldList.hs','Data/JSON/Encoding.hs','Data/JSON/Types/ToJSON.hs']]
    for path in [binary,prlimit,*apis]: bindings[str(path)] = sha(path)
    pairs = []; samples = []; negatives = []; legacy_differences = []
    for index, prior in enumerate(['broad','centered','package']):
        example = next(j['chain'] for j in jobs if j['chain']['prior_label'] == prior)
        text = Path(example['program']).read_text()
        assert text.count(str(Path(example['alignment']).resolve())) == text.count(str(Path(example['tree']).resolve())) == 1
        base = text.replace(str(Path(example['alignment']).resolve()),str(alignment)).replace(str(Path(example['tree']).resolve()),str(tree))
        directories = {}; receipts = {}
        for version, code in [('original_tip_logger',base),('joint_node_logger',transform(base))]:
            program = root/(prior+'-'+version+'.hs'); program.write_text(code)
            config = dict(command=[str(prlimit),'--as='+str(12*2**30),'--cpu=300','--fsize='+str(64*2**20),
                '--',str(binary),'--seed',str(20267001+index),'run',str(program), '--iterations','20',
                '--log-format','json,tsv','--name','independent-chain'],timeout_seconds=300,
                pins={str(p):sha(p) for p in [prlimit,binary,program,alignment,tree,*apis]})
            receipt_path = run_attempt(root/'native'/prior/version,config)
            receipt = json.loads(receipt_path.read_text()); assert receipt['exit_code'] == 0
            dirs = list(receipt_path.parent.glob('independent-chain-*')); assert len(dirs) == 1
            directories[version] = dirs[0]; receipts[version] = str(receipt_path)
            bindings[str(receipt_path)] = sha(receipt_path); bindings[str(program)] = sha(program)
            bindings[str(receipt_path.parent.parent/'configuration.json')] = sha(receipt_path.parent.parent/'configuration.json')
            for name,h in receipt['artifacts'].items(): bindings[str(receipt_path.parent/name)] = h
            print('native_logger_fixture',prior,version,'exited_zero',flush=True)
        left,right = [directories[v] for v in ['original_tip_logger','joint_node_logger']]
        # Exact native comparisons include every scalar row, all recorded
        # ordinary and nonfinite parameter tokens, and every alignment draw.
        for name in ['C1.log','C1.log.column-map.json','C1.log.json','runtime-tree.nwk','C1.P1.fastas']:
            assert (left/name).read_bytes() == (right/name).read_bytes(), ('Logger changed same-seed output',prior,name)
        labels,tips = tree_labels((right/'runtime-tree.nwk').read_text())
        assert len(labels) == 9 and len(tips) == 5
        with (left/'C1.P1.site-property-samples.jsonl').open() as handle: old = [strict_json(line) for line in handle]
        with (right/'C1.P1.site-property-samples.jsonl').open() as handle: new = [strict_json(line) for line in handle]
        assert [f['iter'] for f in old] == [f['iter'] for f in new] == [0,10,20]
        with (right/'C1.P1.fastas').open() as handle: alignments = list(native_alignments(handle))
        for previous, frame, (iteration, sequences) in zip(old,new,alignments):
            assert set(sequences) == labels
            assert set(previous['catStates']) == tips
            assert {tip:frame['catStates'][tip] for tip in tips} == previous['catStates']
            assert frame['properties'] == previous['properties'] and frame['conditions'] == previous['conditions']
            joint = fasta_records(frame['alignmentLines'])
            assert set(joint) == labels
            assert {label:[c == '-' for c in seq] for label,seq in joint.items()} == {label:[c == '-' for c in seq] for label,seq in sequences.items()}
            legacy_differences.extend(dict(prior=prior,iteration=iteration,label=label,legacy=sequences[label],joint=joint[label]) for label in labels if sequences[label] != joint[label])
            samples.append(dict(prior=prior, iteration=iteration, **validate_frame(frame,iteration,joint,tips)))
        pairs.append(dict(prior=prior, seed=20267001+index, original_receipt=receipts['original_tip_logger'],
            full_node_receipt=receipts['joint_node_logger'], exact_complete_scalar_and_alignment_bytes_equal=True,
            original_tip_category_states_unchanged=True, native_nodes=9,native_ancestors=4))
        ancestor = sorted(labels-tips)[0]; frame=new[0]; iteration=alignments[0][0]; sequences=fasta_records(frame['alignmentLines'])
        for case in ['missing_ancestor','invented_node','ancestor_state_letter','ancestor_category_support','ancestor_length','false_iteration','missing_alignment_lines','empty_alignment_lines','unsafe_alignment_line','alignment_state_mismatch','extra_alignment_label']:
            bad=copy.deepcopy(frame)
            if case == 'missing_ancestor': bad['catStates'].pop(ancestor)
            elif case == 'invented_node': bad['catStates']['invented']=bad['catStates'][ancestor]
            elif case == 'ancestor_state_letter': bad['catStates'][ancestor]['states'][0]=(bad['catStates'][ancestor]['states'][0]+1)%20
            elif case == 'ancestor_category_support': bad['catStates'][ancestor]['categories'][0]=4
            elif case == 'ancestor_length': bad['catStates'][ancestor]['states'].pop()
            elif case == 'false_iteration': bad['iter']=1
            elif case == 'missing_alignment_lines': bad.pop('alignmentLines')
            elif case == 'empty_alignment_lines': bad['alignmentLines']=[]
            elif case == 'unsafe_alignment_line': bad['alignmentLines'][0]+='\n'
            elif case == 'alignment_state_mismatch':
                position=bad['alignmentLines'].index('>'+ancestor)+1
                value=bad['alignmentLines'][position]
                offset=next(i for i,c in enumerate(value) if c!='-')
                bad['alignmentLines'][position]=value[:offset]+('C' if value[offset]!='C' else 'A')+value[offset+1:]
            elif case == 'extra_alignment_label': bad['alignmentLines'].extend(['>invented','A'])
            else: raise AssertionError(case)
            rejected(lambda:validate_frame(bad,iteration,sequences,tips)); negatives.append(prior+':'+case)
    write_json(root/'paired_native_results.json',pairs); write_json(root/'native_frame_results.json',samples)
    for path in [root/'paired_native_results.json',root/'native_frame_results.json',
                 Path('scripts/baliphy_joint_node_logger_v3.py'),Path('scripts/baliphy_full_node_logger.py'),Path('scripts/independent_short_sampler_outputs.py'),
                 Path('scripts/independent_native_ancestral_alignment.py')]: bindings[str(path)] = sha(path)
    receipt = dict(status='passed_joint_node_logger_v3_static_and_native_software_qualification',
        checked_utc=datetime.now(timezone.utc).isoformat(),full_programs_checked=405,full_roles_checked=1620,
        deterministic_disjoint_future_seeds=1620,original_seed_namespace_size=len(used),
        native_joint_logger_fixture_priors=['broad','centered','package'],native_joint_logger_fixture_runs=6,
        exact_same_seed_scalar_and_alignment_pairs=3,full_node_saved_frames=9,full_node_ancestor_frames=36,same_record_alignment_state_correspondence_checked=True,legacy_separate_logger_sequence_differences=legacy_differences,
        full_grid_input_labels_safe_for_native_json=len(input_labels),native_ancestor_category_state_pairs=sum(x['ancestral_pairs'] for x in samples),
        original_tip_records_preserved=True,malformed_ancestor_records_rejected=negatives,
        source_hashes=bindings,scientific_eligibility=False,posterior_qualified=False,
        scope='All405existingprograms and1620realrolemetadata statically checked; three exact reversible logger-only edits. Six capped synthetic five-tip20iteration native runs acrossallthreepriors compare identical seeds: complete scalar/JSON/column-map/runtime-tree/FASTA bytes and tip category/state/property records identical; all9same-recordfull-nodeframes and36ancestorframes checked against the alignment encoded in the same record. Legacy separate-logger alignments remain retained and differences explicitly reported; their joint identity is not asserted. Complete latent likelihood, model adequacy, root/tree uncertainty, allocation repair, longer sampling and production horizon remain unqualified. No fungal pilot, existing source/job mutation, GPU or charges.')
    with args.receipt.open('x') as handle: handle.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'}),flush=True)


if __name__ == '__main__': main()
