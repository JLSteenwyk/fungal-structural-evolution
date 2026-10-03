"""Paired AA/structural-state resampling and lossless native archive contracts.

Fixed trees and already predicted states remain conditioning inputs. Circular
blocks are adjacent retained alignment columns, not physical residue distances.
"""
from collections import Counter
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tarfile

import numpy as np

from ancestral_chain_attempt import sha
from background_measurement_union_sources import closed_source
from matched_predictor_branch_fits import load as load_inputs, ROLES, execute
from matched_predictor_branch_inputs import verify
from readback_matched_predictor_branch_inputs import raw_fasta

SCHEMA='matched-predictor-paired-resampling-20261003-v1'
MODES=['iid_sites','circular_blocks10']
SUMMARY_FIELDS=['original_comparison_cases','unique_inputs','replicates_per_mode','modes',
                'resampling_cases','native_roles','native_status_counts','complete_cases','unresolved_cases',
                'serialized_branch_value_slots','finite_branch_values','archived_native_files']


def seed(key,mode,replicate):
    return int.from_bytes(hashlib.sha256((SCHEMA+':'+key+':'+mode+':'+str(replicate)).encode()).digest()[:16],'big')


def indices(length,key,mode,replicate):
    assert isinstance(length,int) and 0<length<65536 and mode in MODES
    assert isinstance(replicate,int) and replicate>=0
    generator=np.random.Generator(np.random.PCG64(seed(key,mode,replicate)))
    if mode=='iid_sites':values=generator.integers(length,size=length)
    else:
        starts=generator.integers(length,size=(length+9)//10)
        values=((starts[:,None]+np.arange(10))%length).reshape(-1)[:length]
    return values.astype('<u2')


def case_id(key,mode,replicate,draw):
    value=json.dumps([SCHEMA,key,mode,replicate,sha_array(draw)],separators=(',',':'))
    return hashlib.sha256(value.encode()).hexdigest()


def sha_array(draw):
    return hashlib.sha256(np.ascontiguousarray(draw,dtype='<u2').tobytes()).hexdigest()


def load(plan,path):
    source,bindings=load_inputs(plan,path)
    point=closed_source(plan['point_completion'],'complete_verified_full_matched_predictor_native_point_fits',
        'complete_verified_full_matched_predictor_native_point_fits_archive',2,bindings)
    assert point['native_roles']==931 and point['unresolved_inputs']==0
    assert plan['replicates_per_mode']==200 and plan['modes']==MODES
    source['data']={};source['resampling_groups']={}
    for key,cfg in source['configs'].items():
        folder=source['root']/'inputs'/key
        data={label:raw_fasta(folder/(label+'.faa')) for label in ['aa','AlphaFold','ESMFold']}
        assert all(set(v)==set(cfg['taxa']) for v in data.values())
        assert all(len(s)==len(cfg['columns']) for v in data.values() for s in v.values())
        for t in cfg['taxa']:
            assert [x=='?' for x in data['aa'][t]]==[x=='?' for x in data['AlphaFold'][t]]==[x=='?' for x in data['ESMFold'][t]]
        source['data'][key]=data
        source['resampling_groups'][key]=hashlib.sha256(json.dumps([SCHEMA,'alignment_draw_group',
            cfg['marker'],cfg['taxa'],cfg['columns'],cfg['alignment_sha256']],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    source['point_completion']=point;verify(bindings);return source,bindings


def build_input(source,key,mode,replicate,root):
    original=source['configs'][key];group=source['resampling_groups'][key]
    draw=indices(len(original['columns']),group,mode,replicate)
    identity=case_id(key,mode,replicate,draw);folder=root/'inputs'/identity;folder.mkdir(parents=True,exist_ok=False)
    for label,rows in source['data'][key].items():
        text=''.join('>'+t+'\n'+''.join(rows[t][int(i)] for i in draw)+'\n' for t in original['taxa'])
        (folder/(label+'.faa')).write_text(text)
    shutil.copyfile(source['root']/'inputs'/key/'topology.nwk',folder/'topology.nwk')
    cfg=dict(original,input_id=identity,columns=[original['columns'][int(i)] for i in draw],
        alignment_sha256={label:sha(folder/(label+'.faa')) for label in source['data'][key]},
        original_input_id=key,resampling_mode=mode,replicate=replicate,resampling_group=group,resampling_seed=seed(group,mode,replicate),
        drawn_source_indices_sha256=sha_array(draw),schema=SCHEMA)
    cfg['topology_sha256']=sha(folder/'topology.nwk')
    with (folder/'config.json').open('x') as f:f.write(json.dumps(cfg,indent=2)+'\n')
    return identity,cfg,draw


def values(rows,splits):
    assert len(rows)==7 and {r['role'] for r in rows}=={r[0] for r in ROLES}
    by_role={r['role']:r for r in rows};status=np.zeros(7,dtype='uint8')
    result=np.full((7,len(splits)),np.nan,dtype='<f8')
    for ri,role in enumerate(ROLES):
        row=by_role[role[0]]
        if row['status']=='native_unsuccessful_retained':status[ri]=1;assert row['point_estimate'] is None
        elif row['status']=='native_output_requires_review':status[ri]=2;assert row['point_estimate'] is None
        else:
            assert row['status']=='native_point_output_integrity_checked_not_model_qualified'
            points={r['split_mask_hex']:r['length'] for r in row['point_estimate']['branches']}
            assert set(points)==set(splits)
            result[ri]=[points[k] for k in splits]
            assert np.isfinite(result[ri]).all() and np.all(result[ri]>=0)
    return status,result


def read_arrays(path,draw,splits,expected_status=None,expected_lengths=None):
    with np.load(path,allow_pickle=False) as a:
        assert len(a.files)==len(set(a.files))==3 and set(a.files)=={'source_indices','role_status','branch_lengths'}
        out={k:a[k] for k in a.files}
    assert out['source_indices'].dtype==np.dtype('<u2') and out['source_indices'].shape==draw.shape
    assert np.array_equal(out['source_indices'],draw)
    assert out['role_status'].dtype==np.dtype('uint8') and out['role_status'].shape==(7,)
    assert np.all(out['role_status']<=2)
    assert out['branch_lengths'].dtype==np.dtype('<f8') and out['branch_lengths'].shape==(7,len(splits))
    for code,row in zip(out['role_status'],out['branch_lengths']):
        assert np.isfinite(row).all() and np.all(row>=0) if code==0 else np.isnan(row).all()
    if expected_status is not None:assert np.array_equal(out['role_status'],expected_status)
    if expected_lengths is not None:assert np.array_equal(out['branch_lengths'],expected_lengths,equal_nan=True)
    return out


def pack(work,path):
    members={str(p.relative_to(work)):sha(p) for p in work.rglob('*') if p.is_file()}
    assert all(not p.is_symlink() for p in work.rglob('*'))
    with tarfile.open(path,'w:gz',compresslevel=3) as archive:
        for name in sorted(members):archive.add(work/name,arcname=name,recursive=False)
    with tarfile.open(path,'r:gz') as archive:
        seen=set()
        for member in archive:
            assert member.isfile() and member.name in members and member.name not in seen
            value=archive.extractfile(member)
            assert hashlib.sha256(value.read()).hexdigest()==members[member.name]
            seen.add(member.name)
        assert seen==set(members)
    return members


def unpack(path,destination,expected):
    assert not destination.exists();destination.mkdir(parents=True)
    seen=set()
    with tarfile.open(path,'r:gz') as archive:
        for member in archive:
            p=PurePosixPath(member.name)
            assert member.isfile() and not p.is_absolute() and '..' not in p.parts
            assert member.name in expected and member.name not in seen
            raw=archive.extractfile(member).read();assert hashlib.sha256(raw).hexdigest()==expected[member.name]
            target=destination/member.name;target.parent.mkdir(parents=True,exist_ok=True)
            with target.open('xb') as f:f.write(raw)
            seen.add(member.name)
    assert seen==set(expected)


def perform(plan,source,key,mode,replicate,output):
    group=source['resampling_groups'][key]
    draw=indices(len(source['configs'][key]['columns']),group,mode,replicate)
    identity=case_id(key,mode,replicate,draw);work=output/'work'/identity
    assert not work.exists();work.mkdir(parents=True)
    newkey,cfg,draw=build_input(source,key,mode,replicate,work);assert newkey==identity
    (work/'roles').mkdir();(work/'native').mkdir();rows=[]
    for role in ROLES:
        rows.append(execute(plan,work,work,identity,cfg,role,source['axes']['positions']))
    splits=source['splits'][key];status,lengths=values(rows,splits)
    destination=output/'cases'/key/mode;destination.mkdir(parents=True,exist_ok=True)
    stem=f'{replicate:04d}';npz=destination/(stem+'.npz');archive=destination/(stem+'.tar.gz')
    with npz.open('xb') as f:np.savez_compressed(f,source_indices=draw,role_status=status,branch_lengths=lengths)
    read_arrays(npz,draw,splits,status,lengths)
    assert not archive.exists();members=pack(work,archive)
    # Only this new declared scratch directory is removed, after every byte is
    # recovered and verified in the retained archive; failed evidence is included.
    shutil.rmtree(work)
    result=dict(schema=SCHEMA,case_id=identity,original_input_id=key,mode=mode,replicate=replicate,
        original_config_sha256=sha(source['root']/'inputs'/key/'config.json'),
        original_native_work_root=str(work.resolve()),resampling_group=group,seed=seed(group,mode,replicate),source_indices_sha256=sha_array(draw),
        native_role_count=7,native_status_counts=dict(Counter(r['status'] for r in rows)),
        complete=bool(np.all(status==0)),serialized_branch_value_slots=lengths.size,
        finite_branch_values=int(np.isfinite(lengths).sum()),native_seconds_sum=sum(r['elapsed_seconds'] for r in rows),
        archive=str(archive),archive_sha256=sha(archive),archive_members=members,array=str(npz),array_sha256=sha(npz),
        scientific_eligibility=False)
    with (destination/(stem+'.json')).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    return result


def add_splits(source,point_root):
    assert Path(source['point_completion']['producer_receipt']).parent.resolve()==point_root.resolve()
    source['splits']={}
    for key in source['configs']:
        row=json.loads((point_root/'roles'/(key+'-aa.json')).read_text())
        source['splits'][key]=[r['split_mask_hex'] for r in row['point_estimate']['branches']]
        assert len(source['splits'][key])==2*len(source['configs'][key]['taxa'])-3
    return source


def runtime_caps(plan,reader=False):
    cpus=2 if reader else plan['resources']['cpus']
    memory=16 if reader else plan['resources']['memory_gib']
    group=next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    folder=Path('/sys/fs/cgroup')/group.lstrip('/')
    observed={key:(folder/key).read_text().strip() for key in ['cpu.max','memory.max','memory.swap.max']}
    assert observed=={'cpu.max':str(cpus*100000)+' 100000','memory.max':str(memory*2**30),'memory.swap.max':'0'}
    assert all(os.environ.get(key)=='1' for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    return observed
