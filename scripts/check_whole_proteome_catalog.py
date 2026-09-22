#!/usr/bin/env python3
"""Exercise source policy, latest status, model selection and byte integrity."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    script=Path(__file__).with_name('catalog_whole_proteome_structures.py').resolve()
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary)
        fasta=root/'proteins.faa';fasta.write_text('>p1\nAAAA\n>p2\nCCCC\n>p3\nDDDD\n>p4\nAAAA\n')
        manifest=root/'manifest.tsv';manifest.write_text('taxon_id\tspecies_name\tstudy_role\nT1\tSynthetic fixture\tingroup\n')
        digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        reps=root/'representatives.json'
        reps.write_text(json.dumps(dict(status='complete',processed_taxa=1,taxa=[dict(taxon_id='T1',path=str(fasta),sha256=digest(fasta),selected_proteins=4)])))
        def model(name,sequence,confidence,provider='GDM'):
            path=root/(name+'.cif');path.write_text('Fixture bytes; not a coordinate-parser test: '+name)
            return dict(model_id=name,version=1,sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),length=len(sequence),mean_ca_plddt=confidence,provider=provider,tool='AlphaFold Monomer v2.0 pipeline',path=str(path),sha256=digest(path))
        low=model('low','AAAA',70);high=model('high','AAAA',90)
        withdrawn=model('withdrawn','CCCC',99);wrong=model('wrong_source','DDDD',99,'OTHER')
        inventory=root/'inventory.jsonl'
        records=[dict(uniprot_accession=name,status='verified',models=[m]) for name,m in [('one',low),('two',high),('three',withdrawn),('four',wrong)]]
        records.append(dict(uniprot_accession='three',status='failed',models=[]))
        inventory.write_text(''.join(json.dumps(r)+'\n' for r in records))
        plan=dict(output=str(root/'output'),inventory=str(inventory),manifest=str(manifest),representatives=str(reps),expected_taxa=1,expected_proteins=4,provider='GDM',tool='AlphaFold Monomer v2.0 pipeline',pins={str(p):digest(p) for p in (inventory,manifest,reps,script)})
        plan_path=root/'plan.json';plan_path.write_text(json.dumps(plan))
        subprocess.run([sys.executable,str(script),'--plan',str(plan_path)],check=True,capture_output=True)
        receipt=json.loads((root/'output/receipt.json').read_text())
        assert receipt['proteins_screened']==4 and receipt['proteins_linked']==2 and receipt['unique_models']==1
        rows=list(csv.DictReader((root/'output/protein_model_links.tsv').open(),delimiter='\t'))
        assert {r['protein_id'] for r in rows}=={'p1','p4'} and {r['model_id'] for r in rows}=={'high'}
        Path(high['path']).write_text('altered bytes')
        plan['output']=str(root/'tampered');plan_path.write_text(json.dumps(plan))
        failed=subprocess.run([sys.executable,str(script),'--plan',str(plan_path)],capture_output=True,text=True)
        assert failed.returncode!=0 and 'Changed coordinate file' in failed.stderr
    print(json.dumps(dict(status='passed_whole_proteome_catalog_fixture',checks=['best confidence model selected','latest failed accession excluded','other provider excluded','duplicate sequences retain separate protein links','tampered selected coordinates rejected'],script_sha256=hashlib.sha256(script.read_bytes()).hexdigest()),indent=2))


if __name__=='__main__':
    main()
