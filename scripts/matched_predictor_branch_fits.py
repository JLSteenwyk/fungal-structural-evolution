"""Bounded matched-predictor native point fits; no topology search or inference."""
import hashlib
from io import StringIO
import json
import math
from pathlib import Path
import re

from Bio import Phylo

from ancestral_chain_attempt import run_attempt, sha
from background_measurement_union_sources import closed_source
from matched_predictor_branch_inputs import verify
from readback_matched_predictor_branch_inputs import raw_fasta
from run_paired_marker_fits import tree_edges

ROLES = [('aa','aa','LG+F+G4'),
         ('AlphaFold_af','AlphaFold','af'),('ESMFold_af','ESMFold','af'),
         ('AlphaFold_af_empirical','AlphaFold','af_empirical'),('ESMFold_af_empirical','ESMFold','af_empirical'),
         ('AlphaFold_llm','AlphaFold','llm'),('ESMFold_llm','ESMFold','llm')]
SUMMARY_FIELDS = ['input_comparison_cases','unique_ready_inputs','native_roles','native_status_counts',
                  'intact_inputs','unresolved_inputs','point_branch_values','native_seconds_sum']


def load(plan,path):
    bindings=dict(plan['pins']);bindings[str(path)]=sha(path)
    completion=closed_source(plan['input_completion'],'complete_verified_full_matched_predictor_branch_inputs',
        'complete_verified_full_matched_predictor_branch_inputs_archive',2,bindings)
    assert (completion['comparison_cases'],completion['unique_ready_inputs'],completion['future_native_roles'])==(8750,133,931)
    assert completion['scientific_eligibility'] is False
    input_plan=json.loads(Path(plan['input_plan']).read_text());root=Path(input_plan['output'])
    assert completion['producer_receipt']==str(root/'receipt.json')
    manifest=json.loads((root/'input_manifest.json').read_text());assert len(manifest)==133
    axes=json.loads((root/'input_axes.json').read_text());configs={}
    for item in manifest:
        key=item['input_id'];config=json.loads((root/'inputs'/key/'config.json').read_text())
        assert key==config['input_id'] and item['taxa']==len(config['taxa']) and item['columns']==len(config['columns'])
        configs[key]=config
    assert len(configs)==133
    verify(bindings)
    return dict(root=root,configs=configs,axes=axes,completion=completion),bindings


def config(plan,input_root,key,source_config,role):
    name,alphabet,model=role;folder=(Path(input_root)/'inputs'/key).resolve()
    alignment=folder/(alphabet+'.faa');topology=folder/'topology.nwk';binary=Path(plan['executable']).resolve()
    definitions={'af':str(Path(plan['models']['AF']).resolve())+'+G4',
                 'af_empirical':str(Path(plan['models']['AF']).resolve())+'+F+G4',
                 'llm':str(Path(plan['models']['LLM']).resolve())+'+G4'}
    native_model=definitions.get(model,model)
    seed=int(hashlib.sha256(('matched-predictor-fixed-topology-v1:'+key+':'+model).encode()).hexdigest()[:12],16)%2147483646+1
    limits=plan['resources']
    native=[str(binary),'-s',str(alignment),'-st','AA','-m',native_model,'-te',str(topology),'-T','1','--mem','2G',
            '--seed',str(seed),'-keep-ident','--prefix','{attempt}/fit']
    command=['/usr/bin/prlimit','--as='+str(limits['native_address_space_gib']*2**30),
             '--cpu='+str(limits['native_cpu_seconds']),'--fsize='+str(limits['native_per_file_limit_mib']*2**20),'--',*native]
    pins={str(p):sha(p) for p in [Path('/usr/bin/prlimit'),binary,alignment,topology,folder/'config.json']}
    if model!='LG+F+G4':
        p=Path(plan['models']['LLM' if model=='llm' else 'AF']).resolve();pins[str(p)]=sha(p)
    assert sha(alignment)==source_config['alignment_sha256'][alphabet]
    assert sha(topology)==source_config['topology_sha256']
    return dict(command=command,pins=pins,timeout_seconds=limits['native_wall_seconds'],
                input_id=key,role=name,alphabet=alphabet,native_model=native_model,seed=seed,
                source_input_config_sha256=sha(folder/'config.json'),
                scientific_eligibility=False,interpretation='Fixed-topology point fit; identical-protein matched-source control only.')


def decode(folder,source_config,positions):
    tree=folder/'fit.treefile';report=(folder/'fit.iqtree').read_text();log=(folder/'fit.log').read_text()
    taxa=set(source_config['taxa']);edges=tree_edges(tree,taxa)
    mask=sum(1<<positions[t] for t in taxa);saved=[]
    for names,length in edges.items():
        side=sum(1<<positions[t] for t in names);other=mask^side
        split=min([side,other],key=lambda x:(x.bit_count(),x))
        saved.append(dict(split_mask_hex=hex(split),length=length,internal=min(side.bit_count(),other.bit_count())>=2))
    assert {r['split_mask_hex'] for r in saved if r['internal']}==set(source_config['internal_splits'])
    data=re.search(r'Input data:\s*(\d+) sequences with (\d+) amino-acid sites',report)
    assert data and tuple(map(int,data.groups()))==(len(taxa),len(source_config['columns']))
    match=re.search(r'Log-likelihood of the tree:\s*([-+\d.eE]+)',report);assert match
    likelihood=float(match[1]);assert math.isfinite(likelihood)
    reported=[line.strip() for line in report.splitlines() if line.startswith('(') and line.rstrip().endswith(';')]
    assert len(reported)==1
    # The serialized report and treefile are separate output paths.
    tmp=Phylo.read(StringIO(reported[0]),'newick')
    tips=[t.name for t in tmp.get_terminals()];assert len(tips)==len(taxa) and set(tips)==taxa
    notes=[line.strip() for line in (report+'\n'+log).splitlines() if 'WARNING' in line.upper() or 'NOTE:' in line.upper()]
    return dict(log_likelihood_reported=likelihood,branches=sorted(saved,key=lambda r:int(r['split_mask_hex'],16)),
                diagnostic_messages=notes,branch_unit='expected_state_substitutions_per_site',
                zero_branch_values=sum(r['length']==0 for r in saved),scientific_eligibility=False)


def execute(plan,input_root,output_root,key,source_config,role,positions):
    name=role[0];root=Path(output_root)/'native'/key/name
    assert not root.exists(),'No native attempt retry or old-output adoption'
    cfg=config(plan,input_root,key,source_config,role)
    rp=run_attempt(root,cfg);native=json.loads(rp.read_text());folder=rp.parent
    value=dict(input_id=key,role=name,seed=cfg['seed'],native_receipt=str(rp),native_receipt_sha256=sha(rp),
               elapsed_seconds=native['elapsed_seconds'],exit_code=native['exit_code'],native_status=native['status'],
               source_config_sha256=sha(Path(input_root)/'inputs'/key/'config.json'),point_estimate=None,scientific_eligibility=False)
    if native['exit_code']!=0:value['status']='native_unsuccessful_retained'
    else:
        try:
            value['point_estimate']=decode(folder,source_config,positions)
            value['status']='native_point_output_integrity_checked_not_model_qualified'
        except (AssertionError,ValueError,KeyError,FileNotFoundError) as error:
            value.update(status='native_output_requires_review',error_type=type(error).__name__,error_message=str(error))
    record=Path(output_root)/'roles'/(key+'-'+name+'.json')
    with record.open('x') as f:f.write(json.dumps(value,indent=2)+'\n')
    return value
