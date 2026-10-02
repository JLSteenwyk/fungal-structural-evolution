"""Native sequence-only MSA execution and exact original-residue correspondence."""
import base64
import hashlib
import itertools
import json
import os
import signal
import subprocess
import time
from pathlib import Path

import psutil

METHODS = ['mafft_auto', 'famsa_default']
ORDERS = list(itertools.permutations(range(3)))


def digest(data):
    return hashlib.sha256(data).hexdigest()


def fasta(sequences, order):
    return ''.join('>m'+str(i)+'\n'+sequences[i]+'\n' for i in order)


def parse_alignment(data, sequences):
    records={};name=None
    for line in data.decode('ascii').splitlines():
        if line.startswith('>'):
            name=line[1:].split()[0]
            assert name in ['m0','m1','m2'] and name not in records
            records[name]=''
        elif line.strip():
            assert name is not None
            records[name]+=line.strip().upper()
    assert set(records)=={'m0','m1','m2'}
    aligned=[records['m'+str(i)] for i in range(3)]
    assert len({len(s) for s in aligned})==1 and len(aligned[0])>0
    for text,original in zip(aligned,sequences):
        assert set(text)<=set(original)|{'-'} and text.replace('-','')==original
    positions=[0,0,0];triples=[]
    for column in zip(*aligned):
        assert any(c!='-' for c in column)
        for i,c in enumerate(column):
            if c!='-':positions[i]+=1
        if '-' not in column:triples.append(positions.copy())
    assert positions==list(map(len,sequences))
    return dict(alignment_columns=len(aligned[0]),common_residues=len(triples),common_full_triples=triples)


def command(plan,method,input_path):
    if method=='mafft_auto':
        return [plan['tools']['mafft'],'--amino','--anysymbol','--thread','1','--auto',str(input_path)]
    assert method=='famsa_default'
    return [plan['tools']['famsa'],'-t','1',str(input_path),'STDOUT']


def execute(plan,record,method,order,root):
    sid=record['sequence_set_id'];permutation=''.join(map(str,order))
    folder=Path(root)/'native_work'/sid/(method+'-'+permutation)
    folder.mkdir(parents=True,exist_ok=True)
    src=folder/'input.faa';text=fasta(record['sequences'],order)
    assert set(folder.iterdir())<=({src} if src.exists() else set())
    if src.exists():assert src.read_text()==text
    else:src.write_text(text)
    cmd=command(plan,method,src);start=time.monotonic()
    env=dict(os.environ,MAFFT_BINARIES=plan['tools']['mafft_binaries'],
             OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env,start_new_session=True)
    created=psutil.Process(proc.pid).create_time()
    try:observed_cmd=psutil.Process(proc.pid).cmdline()
    except (psutil.NoSuchProcess,psutil.ZombieProcess):observed_cmd=[]
    timeout=False
    try:stdout,stderr=proc.communicate(timeout=plan['native_timeout_seconds'])
    except subprocess.TimeoutExpired:
        timeout=True
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        stdout,stderr=proc.communicate()
    result=dict(sequence_set_id=sid,method=method,permutation=permutation,
        models=record['models'],original_lengths=record['original_lengths'],
        input_fasta=text,input_sha256=digest(text.encode()),native_command=cmd,
        native_pid=proc.pid,native_created=created,native_observed_cmdline=observed_cmd,
        returncode=proc.returncode,timed_out=timeout,elapsed_seconds=time.monotonic()-start,
        stdout_base64=base64.b64encode(stdout).decode(),stderr_base64=base64.b64encode(stderr).decode(),
        stdout_sha256=digest(stdout),stderr_sha256=digest(stderr),
        alignment_columns=None,common_residues=None,common_full_triples=None)
    if timeout:result['status']='native_timeout'
    elif proc.returncode:result['status']='native_nonzero_exit'
    else:
        try:result.update(parse_alignment(stdout,record['sequences']));result['status']='valid_sequence_alignment'
        except (AssertionError,UnicodeError,ValueError,IndexError):result['status']='invalid_native_alignment'
    # Raw exact input/output/error bytes and command are stored durably in the
    # SQLite checkpoint; temporary native files are not the provenance archive.
    src.unlink();folder.rmdir()
    return result


def validate_checkpoint(plan,record,result,root):
    order=tuple(map(int,result['permutation']));assert order in ORDERS and result['method'] in METHODS
    assert result['sequence_set_id']==record['sequence_set_id'] and result['models']==record['models']
    assert result['original_lengths']==record['original_lengths']
    expected=fasta(record['sequences'],order);assert result['input_fasta']==expected
    assert result['input_sha256']==digest(expected.encode())
    path=Path(root)/'native_work'/record['sequence_set_id']/(result['method']+'-'+result['permutation'])/'input.faa'
    assert result['native_command']==command(plan,result['method'],path)
    stdout=base64.b64decode(result['stdout_base64'],validate=True)
    stderr=base64.b64decode(result['stderr_base64'],validate=True)
    assert digest(stdout)==result['stdout_sha256'] and digest(stderr)==result['stderr_sha256']
    if result['timed_out']:status='native_timeout';metrics=None
    elif result['returncode']:status='native_nonzero_exit';metrics=None
    else:
        try:metrics=parse_alignment(stdout,record['sequences']);status='valid_sequence_alignment'
        except (AssertionError,UnicodeError,ValueError,IndexError):metrics=None;status='invalid_native_alignment'
    assert result['status']==status
    for field in ['alignment_columns','common_residues','common_full_triples']:
        assert result[field]==(metrics[field] if metrics is not None else None)
    assert isinstance(result['native_pid'],int) and result['native_pid']>0
    assert result['native_created']>0 and result['elapsed_seconds']>=0
