#!/usr/bin/env python3
"""End-to-end PAE readback checks with an asymmetric hand-specified matrix."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np

SCRIPT=Path(__file__).with_name('readback_afdb_pae_qualification.py').resolve()


def main():
    cases=[]
    for case in ['valid','maximum','invalid_sentinel','coordinate','pae_provenance','json_checksum']:
        with tempfile.TemporaryDirectory(prefix='afdb-pae-readback-') as temporary:
            root=Path(temporary)
            def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
            def js(p,data):p.write_text(json.dumps(data))
            def table(p,row):
                with p.open('w') as f:
                    w=csv.DictWriter(f,list(row),delimiter='\t');w.writeheader();w.writerow(row)
            folders={name:root/name for name in ['snapshot','coordinates','qualified','pae']}
            for folder in folders.values():folder.mkdir()
            sequence='ACDEFG';seqsha=hashlib.sha256(sequence.encode()).hexdigest()
            model=dict(model_id='AF-SYNTHETIC-F1',version=1,length=6,sequence_sha256=seqsha,pae_url='https://alphafold.ebi.ac.uk/files/AF-SYNTHETIC-F1-predicted_aligned_error_v1.json')
            js(folders['snapshot']/'model_provenance.json',[model])
            js(folders['snapshot']/'receipt.json',dict(artifacts={'model_provenance.json':sha(folders['snapshot']/'model_provenance.json')}))
            arrays=dict(sequence=np.array(sequence),valid=np.array([False,True,True,True,True,False]),partner_residue_1based=np.array([0,4,4,2,2,0]),ca_plddt=np.full(6,80.),feature_min_plddt=np.full(6,80.))
            before=folders['coordinates']/'model.npz';np.savez_compressed(before,**arrays)
            new={k:v.copy() for k,v in arrays.items()};new['feature_max_pae']=np.array([np.nan,13,0,13,13,np.nan])
            if case=='maximum':new['feature_max_pae'][1]=2
            if case=='invalid_sentinel':new['feature_max_pae'][0]=0
            if case=='coordinate':new['ca_plddt'][2]=79
            after=folders['qualified']/'model.npz';np.savez_compressed(after,**new)
            counts=dict(length=6,valid_states=4,invalid_states=2,valid_focal_plddt70=4,valid_feature_plddt70=4)
            prior=dict(model_id=model['model_id'],version=1,sequence_sha256=seqsha,**counts,encoding_path=str(before),encoding_sha256=sha(before))
            table(folders['coordinates']/'model_summary.tsv',prior)
            qualified=dict(prior,encoding_path=str(after),encoding_sha256=sha(after),valid_feature_plddt70_pae10=1)
            table(folders['qualified']/'model_summary.tsv',qualified)
            mapping_sha=sha(folders['snapshot']/'receipt.json')
            js(folders['coordinates']/'receipt.json',dict(status='complete_native_3di_coordinate_audit',mapping_receipt_sha256=mapping_sha,artifacts={'model_summary.tsv':sha(folders['coordinates']/'model_summary.tsv')}))
            pae=np.zeros((6,6));pae[4,0]=13;pae[0,4]=2
            raw=json.dumps([dict(predicted_aligned_error=pae.tolist(),max_predicted_aligned_error=13)]).encode()
            raw_path=folders['pae']/'model.json.gz';raw_path.write_bytes(gzip.compress(raw,mtime=0))
            pae_row=dict(model,status='verified',url=model['pae_url'],path=str(raw_path),gzip_sha256=sha(raw_path),json_sha256=hashlib.sha256(raw).hexdigest())
            if case=='pae_provenance':pae_row['sequence_sha256']='wrong'
            if case=='json_checksum':pae_row['json_sha256']='wrong'
            js(folders['pae']/'pae_manifest.json',[pae_row])
            js(folders['pae']/'receipt.json',dict(mapping_receipt_sha256=mapping_sha,models_failed=0,artifacts={'pae_manifest.json':sha(folders['pae']/'pae_manifest.json')}))
            js(folders['qualified']/'receipt.json',dict(status='complete_native_3di_feature_audit',mapping_receipt_sha256=mapping_sha,coordinate_audit_receipt_sha256=sha(folders['coordinates']/'receipt.json'),pae_receipt_sha256=sha(folders['pae']/'receipt.json'),models=1,totals=dict(counts,valid_feature_plddt70_pae10=1),artifacts={'model_summary.tsv':sha(folders['qualified']/'model_summary.tsv')}))
            command=[sys.executable,str(SCRIPT)]
            for name,folder in folders.items():command+=['--'+name,str(folder)]
            output=root/'readback.json';command+=['--output',str(output)]
            proc=subprocess.run(command,capture_output=True,text=True,timeout=30)
            if case=='valid':
                assert proc.returncode==0,proc.stderr
                assert json.loads(output.read_text())['totals']['valid_feature_plddt70_pae10']==1
            else:assert proc.returncode!=0 and not output.exists(),case+' was accepted'
            cases.append(dict(case=case,passed=True))
    print(json.dumps(dict(status='passed_afdb_pae_readback_synthetic_checks',cases=cases),indent=2))


if __name__=='__main__':main()
