#!/usr/bin/env python3
"""Native-success, checkpoint, exclusions and failure fixtures for pair execution."""
import hashlib,json,math,subprocess,tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from duplication_alignment_inputs import render_ca
from run_duplication_alignments import run_job
from run_ortholog_pair_guide_comparison import sha

with tempfile.TemporaryDirectory() as td:
    root=Path(td);seq='ACDEFGHIKLMNPQRSTVWY';record=dict(status='validated',sequence=seq,ca_xyz=[[3*math.cos(i),3*math.sin(i),i] for i in range(len(seq))],ca_plddt=[90]*len(seq))
    blob,_,_=render_ca(record);a=root/'a.pdb';b=root/'b.pdb';a.write_bytes(blob);b.write_bytes(blob)
    inputs={('A',1,'full'):dict(status='ready',path=str(a),sha256=sha(a),sequence=seq),('B',1,'full'):dict(status='ready',path=str(b),sha256=sha(b),sequence=seq),('C',1,'full'):dict(status='too_few_retained_residues')}
    plan={'output':str(root/'out'),'usalign':'/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign','options':['-mol','prot','-mm','0','-outfmt','0','-ter','2'],'per_pair_timeout_seconds':10}
    def job(tag,right='B'):return(hashlib.sha256(tag.encode()).hexdigest(),('A',1),(right,1),'full',0)
    dest,status=run_job(job('native'),inputs,plan,'plan','manifest');assert status=='aligned';r=json.loads(dest.read_text());assert r['metrics']['tm_left']>.999
    with patch('run_duplication_alignments.subprocess.run',side_effect=AssertionError('Unexpected rerun')):
        assert run_job(job('native'),inputs,plan,'plan','manifest')==(dest,'aligned')
        assert run_job(job('excluded','C'),inputs,plan,'plan','manifest')[1]=='input_unavailable'
    for tag,native,expected in [('error',SimpleNamespace(returncode=7,stdout='',stderr='fixture'),'native_error'),('parse',SimpleNamespace(returncode=0,stdout='invalid',stderr=''),'parse_error')]:
        with patch('run_duplication_alignments.subprocess.run',return_value=native):assert run_job(job(tag),inputs,plan,'plan','manifest')[1]==expected
    with patch('run_duplication_alignments.subprocess.run',side_effect=subprocess.TimeoutExpired('fixture',10,output=b'partial')):assert run_job(job('timeout'),inputs,plan,'plan','manifest')[1]=='timeout'
    bad=json.loads(dest.read_text());bad['metrics']['tm_left']=0;dest.write_text(json.dumps(bad))
    try:run_job(job('native'),inputs,plan,'plan','manifest')
    except ValueError:pass
    else:raise AssertionError('Changed metrics accepted')
    a.write_bytes(blob+b'\n')
    try:run_job(job('changed'),inputs,plan,'plan','manifest')
    except ValueError:pass
    else:raise AssertionError('Changed source accepted')
print('Native alignment, cached reuse, short-input exclusion, native/parse/timeout dispositions passed; altered metrics and coordinates rejected.')
