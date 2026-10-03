#!/usr/bin/env python3
"""Full-grid serialization plus adversarial matched-mask/topology contracts."""
import argparse
import copy
from datetime import datetime, timezone
import gzip
from io import StringIO
import json
from pathlib import Path
import shutil

from Bio import Phylo

from audit_selected_taxon_identity_snapshot_v2 import sha
from matched_predictor_branch_inputs import load, joint_rows, verify
from prepare_matched_predictor_branch_inputs import produce
from readback_matched_predictor_branch_inputs import (check_grid, raw_splits, independent_marker,
                                                    expected_case, check_ready_input)


def rejected(action):
    try: action()
    except (AssertionError,ValueError,KeyError,FileNotFoundError): return
    raise AssertionError('Malformed input accepted')


def fixture(root):
    names = list('abcdefghi'); positions = {t:i for i,t in enumerate(names)}
    raw = '((a:1,b:1):2,(c:1,d:1):3,((e:1,f:1):4,(g:1,(h:1,i:1):5):6):7);'
    source_tree = root/'synthetic_raw_tree.nwk'; source_tree.write_text(raw)
    memberships = [names,names[:-1],names[2:],list('abhi'),list('acfh')]
    views = []
    for vi in range(70):
        taxa = memberships[vi%len(memberships)]
        tree = Phylo.read(StringIO(raw),'newick')
        for tip in list(tree.get_terminals()):
            if tip.name not in taxa: tree.prune(tip)
        branches = [dict(split_mask_hex=hex(k),branch_length=1.,branch_length_unit='synthetic')
                    for k in sorted(raw_splits(tree,taxa,positions))]
        views.append(dict(cohort='synthetic_'+str(vi%5),view='synthetic_'+str(vi),family='synthetic',
                          source_tree=str(source_tree),prune=taxa!=names,taxa=taxa,branches=branches))
    markers = [f'synthetic_{i:03d}' for i in range(125)]; contexts = {}; alignments = {}; required = {}
    for mi, marker in enumerate(markers):
        if mi%8 == 7: alignments[marker] = None; required[marker] = None; continue
        count = [9,8,6,4,3,1,0][mi%8]; eligible = set(names[:count])
        left = dict(aa={}, **{'3di':{}}); right = dict(aa={}, **{'3di':{}})
        pa,pb = list(range(1,65)),list(range(2,66))
        for t in names:
            left['aa'][t] = 'A'*64; left['3di'][t] = 'G'*64
            observed = 63 if t in eligible else 20
            right['aa'][t] = 'A'*observed+'?'*(64-observed)
            right['3di'][t] = 'K'*observed+'?'*(64-observed)
            contexts[marker,t] = dict(jointly_observed_positions=str(observed),same_aa_positions=str(observed),state_mismatches=str(observed))
        alignments[marker] = [(pa,left),(pb,right)]; required[marker] = 50
    return dict(views=views,positions=positions,markers=markers,taxa=[dict(taxon_id=t,study_role='synthetic') for t in names],
                contexts=contexts,alignments=alignments,required=required)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists() and not a.receipt.exists(); a.output.mkdir(parents=True)
    source = fixture(a.output); root=a.output/'full_grid';root.mkdir()
    producer, _ = produce(source,root,verbose=False); reader = check_grid(source,root,verbose=False)
    assert producer == reader and producer['comparison_cases'] == 8750
    assert all(n > 0 for n in producer['branch_status_counts'].values())
    first = source['markers'][0]
    columns,data,_ = independent_marker(source,first)
    record, alignments, splits = expected_case(source,0,first,data)
    key = record['input_id']; negatives=[]
    for label in ['different_aa','different_mask','wrong_columns']:
        pair=copy.deepcopy(source['alignments'][first])
        if label=='different_aa':pair[1][1]['aa']['a']='C'+pair[1][1]['aa']['a'][1:]
        elif label=='different_mask':pair[1][1]['3di']['a']='?'+pair[1][1]['3di']['a'][1:]
        else:pair[1][0][0]=99
        ctx={t:r for (m,t),r in source['contexts'].items() if m==first}
        rejected(lambda:joint_rows(*pair,ctx,50));negatives.append(label)
    with gzip.open(root/'comparison_cases.jsonl.gz','rt') as f: cases=[json.loads(line) for line in f]
    for label in ['missing_case','duplicate_case','wrong_status','wrong_identity','changed_column','changed_taxa','changed_original_branch_map','false_branch_count']:
        bad=copy.deepcopy(cases)
        if label=='missing_case':bad.pop()
        elif label=='duplicate_case':bad[-1]=bad[0]
        elif label=='wrong_status':bad[0]['status']='accepted'
        elif label=='wrong_identity':bad[0]['input_id']='0'*64
        elif label=='changed_column':bad[0]['columns'][0]+=1
        elif label=='changed_taxa':bad[0]['taxa'][0]='invented'
        elif label=='changed_original_branch_map':bad[0]['internal_branch_mapping'][0]['original_branch_indices'].append(999)
        else:bad[0]['branch_status_counts']['unique_internal_projection']+=1
        folder=a.output/label;folder.mkdir()
        for path in root.iterdir():
            if path.name!='comparison_cases.jsonl.gz': (folder/path.name).symlink_to(path.resolve(),target_is_directory=path.is_dir())
        with gzip.open(folder/'comparison_cases.jsonl.gz','wt') as f:
            for row in bad:f.write(json.dumps(row)+'\n')
        rejected(lambda:check_grid(source,folder,verbose=False));negatives.append(label)
    for label in ['missing_alignment','changed_aa_export','changed_state_mask','changed_topology','wrong_config_columns','extra_input_file']:
        folder=a.output/label; (folder/'inputs').mkdir(parents=True)
        target=folder/'inputs'/key;shutil.copytree(root/'inputs'/key,target)
        if label=='missing_alignment':(target/'ESMFold.faa').unlink()
        elif label=='changed_aa_export':
            path=target/'aa.faa';path.write_text(path.read_text().replace('A','C',1))
        elif label=='changed_state_mask':
            path=target/'AlphaFold.faa';path.write_text(path.read_text().replace('G','?',1))
        elif label=='changed_topology':(target/'topology.nwk').write_text('(a:0.1,b:0.1,c:0.1,d:0.1);\n')
        elif label=='wrong_config_columns':
            path=target/'config.json';config=json.loads(path.read_text());config['columns'][0]+=1;path.write_text(json.dumps(config))
        else:(target/'invented.txt').write_text('invalid')
        rejected(lambda:check_ready_input(folder,key,record,alignments,splits,source['positions']));negatives.append(label)
    real,bindings=load()
    assert (len(real['views']),len(real['markers']),len(real['positions']),len(real['contexts']))==(70,125,526,673)
    own=['matched_predictor_branch_inputs','prepare_matched_predictor_branch_inputs','readback_matched_predictor_branch_inputs',
         'check_matched_predictor_branch_inputs']
    bindings.update({f'scripts/{name}.py':sha(f'scripts/{name}.py') for name in own});verify(bindings)
    result=dict(status='passed_full_matched_predictor_branch_input_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        synthetic_grid=producer,synthetic_full_serialized_comparison_cases=8750,malformed_cases_rejected=negatives,
        independent_raw_tree_pruning_checked=True,full_real_tree_views=70,full_real_marker_slots=125,full_real_taxon_entries=526,
        exact_complete_protein_context_cells=673,source_hashes=bindings,
        artifacts={str(p):sha(p) for p in a.output.rglob('*') if p.is_file() and not p.is_symlink()},
        scientific_eligibility=False,
        scope='Complete 70-view/125-marker synthetic grid serialized and independently rebuilt; all five branch dispositions and 17 altered mask/identity/mapping/file/topology cases checked. Full actual source scope and 673 complete protein identities bound. Synthetic taxa are nine fixture tips, not the 526-entry fungal study. No native fitting, accepted species framework, predictor experimental error, calibrated branch effect or aim completion.')
    with a.receipt.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2),flush=True)


if __name__=='__main__':main()
