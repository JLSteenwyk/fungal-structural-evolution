#!/usr/bin/env python3
"""Exercise exact recovered ASA union and reject provenance failures using fixtures."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from catalog_whole_proteome_structures import sha


def save(p, x):
    p.write_text(json.dumps(x, indent=2)+'\n')


def main():
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory); snapshot=root/'snapshot'; snapshot.mkdir()
        models=[]; cohorts=[]; expected_entries={}; tables=[]
        for index, names in enumerate([['shared','excluded'],['new']]):
            mapping=root/f'mapping{index}'; asa=root/f'asa{index}'; audit=root/f'audit{index}'
            for p in [mapping,asa,audit]:p.mkdir()
            local=[]
            for sid in names:
                coord=root/(sid+'.cif');coord.write_text('synthetic coordinates '+sid)
                m=dict(model_id=sid,length=2,sequence_sha256='sequence'+sid,sha256=sha(coord),path=str(coord),uniprot_accession=sid)
                local.append(m)
                if sid!='excluded':models.append(dict(m,**({'source_record_id':sid} if sid=='shared' else {})))
            save(mapping/'model_provenance.json',local)
            save(mapping/'receipt.json',dict(artifacts={'model_provenance.json':sha(mapping/'model_provenance.json')}))
            save(asa/'config.json',dict(snapshot_receipt_sha256=sha(mapping/'receipt.json'),models=len(local),residues=2*len(local),workers=1,method='synthetic'))
            entries={}; rows=[]
            for m in local:
                sid=m['model_id'];table=asa/(sid+'.residues.tsv.gz');table.write_bytes(b'synthetic table');tables.append(table)
                rp=asa/(sid+'.receipt.json');save(rp,dict(config_sha256=sha(asa/'config.json'),model_sha256=m['sha256'],sequence_sha256=m['sequence_sha256'],table_sha256=sha(table)))
                entries[sid]=sha(rp);rows.append(sid+'\t2\t'+sha(rp))
                if sid!='excluded':expected_entries[sid]=sha(rp)
            save(asa/'receipt.json',dict(status='complete_snapshot_predicted_accessibility',config_sha256=sha(asa/'config.json'),entry_receipts=entries))
            at=audit/'audited_models.tsv';at.write_text('model_id\tresidues\tentry_receipt_sha256\n'+'\n'.join(rows)+'\n')
            save(audit/'receipt.json',dict(status='passed_full_accessibility_snapshot',assessment_receipt_sha256=sha(asa/'receipt.json'),snapshot_receipt_sha256=sha(mapping/'receipt.json'),models_audited=len(local),residues_audited=2*len(local),artifacts={at.name:sha(at)}))
            cohorts.append(dict(mapping=str(mapping),mapping_receipt_sha256=sha(mapping/'receipt.json'),accessibility=str(asa),audit=str(audit),excluded_model_ids=['excluded'] if index==0 else []))
        save(snapshot/'model_provenance.json',models)
        def snapshot_receipt():save(snapshot/'receipt.json',dict(distinct_models=2,artifacts={'model_provenance.json':sha(snapshot/'model_provenance.json')}))
        snapshot_receipt(); manifest=root/'manifest.json'
        def write_manifest():save(manifest,dict(snapshot_receipt_sha256=sha(snapshot/'receipt.json'),cohorts=cohorts))
        write_manifest()
        arguments=['='.join(c[k] for k in ['mapping','accessibility','audit']) for c in cohorts]
        def run(name,selected=None):
            cmd=[sys.executable,'scripts/merge_recovered_accessibility.py','--snapshot',str(snapshot),'--manifest',str(manifest),'--output',str(root/name)]
            for x in arguments if selected is None else selected:cmd+=['--cohort',x]
            return subprocess.run(cmd,capture_output=True,text=True)
        good=run('valid');assert good.returncode==0,good.stderr
        r=json.loads((root/'valid/receipt.json').read_text());assert r['models']==2 and r['residues']==4 and r['entry_receipts']==expected_entries
        assert not (root/'valid/excluded.receipt.json').exists()
        for sid,h in expected_entries.items():assert (root/'valid'/(sid+'.receipt.json')).is_symlink() and sha(root/'valid'/(sid+'.receipt.json'))==h
        failures=[]
        def fail(name,selected=None):
            result=run(name,selected);assert result.returncode!=0,name
            assert not (root/name/'receipt.json').exists();failures.append(name)
        fail('missing_cohort',arguments[:1]);fail('duplicate_cohort',arguments+arguments[:1])
        cohorts[0]['excluded_model_ids']=[];write_manifest();fail('undeclared_exclusion');cohorts[0]['excluded_model_ids']=['excluded'];write_manifest()
        original=tables[0].read_bytes();tables[0].write_bytes(b'changed');fail('changed_residue_table');tables[0].write_bytes(original)
        for field,value in [('source_record_id','wrong'),('sequence_sha256','wrong')]:
            saved=models[0][field];models[0][field]=value;save(snapshot/'model_provenance.json',models);snapshot_receipt();write_manifest();fail('changed_'+field);models[0][field]=saved
        save(snapshot/'model_provenance.json',models);snapshot_receipt();write_manifest()
        ap=Path(cohorts[1]['audit'])/'receipt.json';ar=json.loads(ap.read_text());ar['status']='incomplete';save(ap,ar);fail('incomplete_audit')
    r=dict(status='passed_recovered_accessibility_union_cli_fixtures',models_selected=2,models_explicitly_excluded=1,rejected_cases=failures,merger_sha256=sha('scripts/merge_recovered_accessibility.py'),checker_sha256=sha(__file__),scope='Synthetic provenance, subset, symlink and rejection behavior only; not production ASA or coordinate validation.')
    Path('metadata/recovered_afdb_accessibility_merge_fixture_checks_20260927.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
