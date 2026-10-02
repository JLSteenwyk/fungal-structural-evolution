#!/usr/bin/env python3
"""Full operator-bank software contracts with synthetic closed sources/journals."""
import argparse
import copy
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from collections import Counter
import numpy as np
from scipy import sparse
from full_entity_operator_sources import KINDS,MODES,digest
from prepare_full_entity_operators import run
from run_ortholog_pair_guide_comparison import sha


def write(path,j):Path(path).write_text(json.dumps(j,indent=2)+'\n')


def table(path,rows):
    with gzip.open(path,'wt') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def eid(kind,token):return hashlib.sha256(json.dumps(['expanded-covariance-'+kind+'-v1',token],separators=(',',':')).encode()).hexdigest()


def setup(root):
    source=root/'covariance';source.mkdir();cases=[];occurrences=[];entities={k:set() for k in KINDS}
    for i in range(24):
        fam='F'+str(i//12);cid=digest(['fixture-case',i]);cases.append(dict(case_id=cid,target_family=fam,background_family=fam,family_component=digest(['component',fam]),species_pattern_row=i%7))
        for side,sign in [('target',1),('background',-1)]:
            tokens={side+'_node':[digest([side,fam,i if side=='target' else i%3])],
                'model_pair':[digest(['pair',side,fam,i%3])],'family':[fam]}
            for kind,token in tokens.items():
                identifier=eid(kind,token);entities[kind].add(identifier)
                occurrences.append(dict(case_id=cid,entity_kind=kind,entity_id=identifier,entity_value=json.dumps(token,separators=(',',':')),side=side,endpoint='',signed_loading=sign,unsigned_loading=1))
            for end in ['a','b']:
                gene=['gene-'+fam+'-'+side+'-'+end+'-'+str(i if side=='target' else i%3)]
                model=['model-'+fam+'-'+side+'-'+end+'-'+str(i%3),1,digest(['coordinate',fam,side,end,i%3])]
                for kind,token in [('gene',gene),('model',model)]:
                    identifier=eid(kind,token);entities[kind].add(identifier)
                    occurrences.append(dict(case_id=cid,entity_kind=kind,entity_id=identifier,entity_value=json.dumps(token,separators=(',',':')),side=side,endpoint=end,signed_loading=sign*.5,unsigned_loading=.5))
    table(source/'case_covariance_index.tsv.gz',cases);table(source/'entity_incidence.tsv.gz',occurrences)
    rng=np.random.default_rng(920);factor=rng.normal(size=(7,3));trees=['t1','t2','t3','t4','t5']
    for i,t in enumerate(trees):np.savez_compressed(source/(t+'.npz'),factor=factor*(.8+.05*i))
    config=root/'covariance-plan.json';write(config,dict(output=str(source),trees=trees))
    summary=dict(logical_cases=24,entity_occurrences=336,family_components=2,trees=trees,unique_entities={k:len(v) for k,v in entities.items()})
    rp=source/'receipt.json';write(rp,dict(status='complete_full_expanded_covariance_pending_independent_readback',plan_sha256=sha(config),scientific_eligibility=False,
        artifacts={p.name:sha(p) for p in source.iterdir() if p.is_file()},**summary))
    bindings={str(p):sha(p) for p in root.rglob('*') if p.is_file()}
    ap=source/'completion_archive.json';write(ap,dict(status='complete_verified_full_expanded_covariance_archive',services=[{},{}],summary=summary,source_hashes=bindings))
    cp=root/'covariance-completion.json';write(cp,dict(status='complete_verified_full_expanded_covariance',exact_process_journals_checked=2,bound_source_hashes=len(bindings),full_hash_archive=str(ap),full_hash_archive_sha256=sha(ap),producer_receipt=str(rp),producer_receipt_sha256=sha(rp),scientific_eligibility=False,**summary))
    plan=root/'plan.json';write(plan,dict(covariance_completion=str(cp),covariance_plan=str(config),expected={k:summary[k] for k in ['logical_cases','entity_occurrences','family_components']},trees=trees,output=str(root/'output'),pins={},resources=dict(minimum_free_disk_gib=0),scope='Synthetic closure/journal source fixture; full software contracts only, no biological pilot.'))
    return plan


def invoke(plan,output,good=True):
    r=subprocess.run([sys.executable,'scripts/readback_full_entity_operators.py','--plan',str(plan),'--output',str(output)],capture_output=True,text=True)
    if good and r.returncode:raise RuntimeError(r.stdout+r.stderr)
    if not good:assert r.returncode!=0,'Altered operator bank accepted'


def mutate(root,name):
    out=root/'output';rp=out/'receipt.json';receipt=json.loads(rp.read_text())
    def changed(fp):receipt['artifacts'][str(fp.relative_to(out))]=sha(fp)
    if name in ['case_order','block_label','entity_labels']:
        fp=out/('case_ids.json' if name=='case_order' else 'block_labels.npy' if name=='block_label' else 'entity_labels.json')
        if name=='case_order':rows=json.load(open(fp));rows.reverse();write(fp,rows)
        elif name=='block_label':a=np.load(fp);a[0]=a[-1];np.save(fp,a)
        else:labels=json.load(open(fp));labels['gene'][0]='wrong';write(fp,labels)
        changed(fp)
    elif name in ['loading_sign','signed_family_promoted','family_intercept_doubled','gene_columns_collapsed','missing_operator']:
        mp=out/'operator_manifest.json';entries=json.load(open(mp))
        if name=='missing_operator':entries.pop()
        else:
            mode,kind=('signed','gene') if name=='loading_sign' else ('signed','family') if name=='signed_family_promoted' else ('contrast','family_intercept') if name=='family_intercept_doubled' else ('unsigned','gene')
            entry=next(e for e in entries if (e['mode'],e['kind'])==(mode,kind));fp=out/entry['path'];z=sparse.load_npz(fp).tolil()
            if name=='loading_sign':z[0,0]=.25
            elif name=='signed_family_promoted':z[0,0]=1.
            elif name=='family_intercept_doubled':z=z*2
            else:z[:,0]=z[:,0]+z[:,1];z[:,1]=0
            z=z.tocsr();z.eliminate_zeros();sparse.save_npz(fp,z);entry.update(sha256=sha(fp),nnz=z.nnz);changed(fp)
        write(mp,entries);changed(mp)
    elif name in ['gram_value','gram_rank','zero_kernel_removed']:
        mp=out/'kernel_gram_manifest.json';rows=json.load(open(mp));entry=rows[0]
        if name=='gram_value':
            fp=out/entry['path']
            with np.load(fp) as a:g=a['gram'].copy()
            g[0,0]+=1;np.savez_compressed(fp,gram=g);entry['sha256']=sha(fp);changed(fp)
        elif name=='gram_rank':entry['normalized_kernel_gram_rank']+=1
        else:entry['zero_kernels']=[]
        write(mp,rows);changed(mp)
    else:
        mp=out/'benchmark_manifest.json';rows=json.load(open(mp));entry=rows[0]
        if name=='benchmark_solution':
            fp=out/entry['path']
            with np.load(fp) as a:v=a['solution'].copy()
            v[0,0]+=.1;np.savez_compressed(fp,solution=v);entry['sha256']=sha(fp);changed(fp)
        elif name=='benchmark_variance':entry['entity_variances']['gene']+=.1
        elif name=='benchmark_logdet':entry['logdet']+=1
        elif name=='omit_final_tree':rows=[r for r in rows if r['tree']!='t5']
        else:raise AssertionError(name)
        write(mp,rows);changed(mp)
    write(rp,receipt)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    changes=['case_order','block_label','entity_labels','loading_sign','signed_family_promoted','family_intercept_doubled','gene_columns_collapsed','missing_operator','gram_value','gram_rank','zero_kernel_removed','benchmark_solution','benchmark_variance','benchmark_logdet','omit_final_tree']
    with tempfile.TemporaryDirectory(prefix='full-entity-bank-contract-') as directory:
        root=Path(directory);plan=setup(root);result=run(plan);invoke(plan,root/'readback.json')
        assert result['matrices']==13 and result['numerical_benchmarks']==10
        assert all('family' in v for k,v in result['zero_kernels'].items() if k.startswith('signed|'))
        try:run(plan)
        except AssertionError:pass
        else:raise AssertionError('Completed bank restart accepted')
        snapshot={str(p.relative_to(root/'output')):p.read_bytes() for p in (root/'output').rglob('*') if p.is_file()}
        for change in changes:
            mutate(root,change);invoke(plan,root/'bad-readback.json',False)
            for p in (root/'output').rglob('*'):
                if p.is_file() and str(p.relative_to(root/'output')) not in snapshot:p.unlink()
            for name,data in snapshot.items():(root/'output'/name).write_bytes(data)
        root2=root/'interrupt';root2.mkdir();pp=setup(root2)
        try:run(pp,stop_after_exports=True)
        except InterruptedError:pass
        else:raise AssertionError('Interruption not exercised')
        run(pp);invoke(pp,root2/'readback.json')
    sources=[Path(__file__),*[Path('scripts')/n for n in ['full_entity_operator_sources.py','prepare_full_entity_operators.py','readback_full_entity_operators.py','shared_entity_covariance.py']]]
    result=dict(status='passed_full_entity_operator_bank_software_contracts',rejected_rehashed_exports=changes,full_interrupt_replay_passed=True,completed_restart_refused=True,
        source_and_journal_fixtures_synthetic=True,all_five_trees_both_loadings_preserved=True,distinct_genes_shared_models_and_cancelling_family_preserved=True,
        source_hashes={str(p):sha(p) for p in sources},scope='24 distinct logical gene contexts, shared model coordinates, reused backgrounds, both loadings, all13operators and10 full-case fixed benchmarks. Independent raw Counter loadings, latent-space Grams and explicit outer-product/LU benchmark reconstruction. Synthetic source/journal software contracts, not a pilot, production fit or inferential calibration.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
