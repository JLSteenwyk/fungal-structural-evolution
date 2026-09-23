#!/usr/bin/env python3
"""Exercise full ASA merge CLI with explicitly synthetic provenance fixtures."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from catalog_whole_proteome_structures import sha


def save(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);snapshot=root/'snapshot';snapshot.mkdir();inventory=root/'inventory';inventory.mkdir()
        (inventory/'inventory.jsonl').write_text('fixture-only\n')
        save(inventory/'receipt.json',{'status':'complete_disjoint_prediction_inventory_union_with_batch_provenance','artifacts':{'inventory.jsonl':sha(inventory/'inventory.jsonl')}})
        models=[];cohorts=[];entries={};targets=[]
        for i in range(2):
            source=root/('source'+str(i));asa=root/('asa'+str(i));audit=root/('audit'+str(i))
            for p in [source,asa,audit]:p.mkdir()
            original=root/('old'+str(i)+'.cif');relocated=root/('new'+str(i)+'.cif');original.write_text('synthetic fixture coordinate bytes '+str(i));relocated.write_bytes(original.read_bytes());targets.append(relocated)
            sid='MODEL'+str(i);m={'model_id':sid,'length':2,'sequence_sha256':'sequence'+str(i),'sha256':sha(original),'path':str(original)}
            models.append(dict(m,path=str(relocated)))
            save(source/'model_provenance.json',[m]);save(source/'receipt.json',{'status':'fixture_source','artifacts':{'model_provenance.json':sha(source/'model_provenance.json')}})
            config={'snapshot_receipt_sha256':sha(source/'receipt.json'),'models':1,'residues':2,'workers':1,'method':'synthetic-fixture-only'};save(asa/'config.json',config)
            table=asa/(sid+'.residues.tsv.gz');table.write_bytes(b'fixture table bytes')
            entry={'config_sha256':sha(asa/'config.json'),'model_sha256':m['sha256'],'sequence_sha256':m['sequence_sha256'],'table_sha256':sha(table)};rp=asa/(sid+'.receipt.json');save(rp,entry);entries[sid]=sha(rp)
            save(asa/'receipt.json',{'status':'complete_snapshot_predicted_accessibility','config_sha256':sha(asa/'config.json'),'entry_receipts':{sid:sha(rp)}})
            at=audit/'audited_models.tsv';at.write_text('model_id\tresidues\tentry_receipt_sha256\n'+sid+'\t2\t'+sha(rp)+'\n')
            save(audit/'receipt.json',{'status':'passed_full_accessibility_snapshot','assessment_receipt_sha256':sha(asa/'receipt.json'),'snapshot_receipt_sha256':sha(source/'receipt.json'),'models_audited':1,'residues_audited':2,'artifacts':{at.name:sha(at)}})
            cohorts.append({'mapping':str(source),'mapping_receipt_sha256':sha(source/'receipt.json'),'accessibility':str(asa),'audit':str(audit)})
        save(snapshot/'model_provenance.json',models);save(snapshot/'receipt.json',{'source_inventory_sha256':sha(inventory/'inventory.jsonl'),'distinct_models':2,'artifacts':{'model_provenance.json':sha(snapshot/'model_provenance.json')}})
        manifest=root/'cohorts.json';save(manifest,{'snapshot_receipt_sha256':sha(snapshot/'receipt.json'),'cohorts':cohorts})
        base=[sys.executable,'scripts/merge_completed_accessibility.py','--snapshot',str(snapshot),'--inventory',str(inventory),'--cohort-manifest',str(manifest)]
        args=['='.join(c[k] for k in ['mapping','accessibility','audit']) for c in cohorts]
        def run(name,selected):
            cmd=base+['--output',str(root/name)]
            for argument in selected:cmd+=['--cohort',argument]
            return subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        valid=run('valid',args)
        if valid.returncode:raise RuntimeError(valid.stderr)
        result=json.loads((root/'valid/receipt.json').read_text());assert result['models']==2 and result['residues']==4 and result['entry_receipts']==entries
        for sid,digest in entries.items():
            link=root/'valid'/(sid+'.receipt.json');assert link.is_symlink() and sha(link)==digest
            assert (root/'valid'/(sid+'.residues.tsv.gz')).is_symlink()
        failures=[]
        for name,selected in [('missing_cohort',args[:1]),('duplicate_cohort',args+args[:1])]:
            assert run(name,selected).returncode!=0;assert not (root/name/'receipt.json').exists();failures.append(name)
        saved=targets[0].read_bytes();targets[0].write_bytes(b'changed')
        assert run('changed_target_bytes',args).returncode!=0;failures.append('changed_target_bytes');targets[0].write_bytes(saved)
        audit=Path(cohorts[0]['audit'])/'receipt.json';r=json.loads(audit.read_text());r['status']='incomplete';save(audit,r)
        assert run('incomplete_source_audit',args).returncode!=0;failures.append('incomplete_source_audit')
    receipt={'status':'passed_complete_accessibility_merge_cli_fixtures','fixture_models':2,'fixture_cohorts':2,'unchanged_relative_symlink_entries':True,'rejected_cases':failures,'merger_sha256':sha('scripts/merge_completed_accessibility.py'),'script_sha256':sha(__file__),'scope':'Synthetic metadata/file-lifecycle fixture, not physical coordinate or ASA validation. Full CLI verifies relocated byte identity, disjoint source coverage, completion gates and unchanged symlinked entry receipts. No production inputs modified.'}
    Path('metadata/completed_accessibility_merge_cli_fixture_checks.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
