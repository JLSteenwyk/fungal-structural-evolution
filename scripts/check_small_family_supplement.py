#!/usr/bin/env python3
"""Validate small-family supplements against native fixtures and corruptions."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from supplement_small_family_orthologs import expected_pairs, scan_native, supplement

ROOT=Path(__file__).resolve().parents[1]


def rejected(call):
    try:call()
    except (ValueError,FileNotFoundError):return
    raise AssertionError('Corruption was accepted')


def main():
    checked=[]
    with tempfile.TemporaryDirectory(prefix='small-family-supplement-') as temporary:
        root=Path(temporary)
        for tag,fixture_name,n_missing in [
            ('sorted','native-small-family-output-fixture-v2',2),
            ('early','native-small-family-early-singleton-fixture-v1',12)]:
            fixture=ROOT/'results/orthology'/fixture_name
            source=fixture/'Source/WorkingDirectory';result=fixture/'Results_small_fixture'
            native_receipt=json.loads((fixture/'receipt.json').read_text())
            receipt=supplement(source,result,root/tag,native_receipt)
            rows=list(csv.DictReader((root/tag/'supplemental_directed_pairs.tsv').open(),delimiter='\t'))
            actual={(r['family'],r['source_protein'],r['target_protein']) for r in rows}
            assert len(rows)==len(actual)==n_missing==receipt['supplemental_directed_pairs']
            assert actual==set(map(tuple,native_receipt['missing_directed_pairs']))
            expected,sizes,names,species,reverse=expected_pairs(source)
            observed,_,_=scan_native(result,sizes,species,reverse,expected)
            supplement_keys={(r['family'],r['source_native_gene'],r['target_native_gene']) for r in rows}
            assert not observed&supplement_keys and observed|supplement_keys==expected
            assert all((family,right,left) in expected for family,left,right in expected)
            checked.append(dict(fixture=tag,supplemental_pairs=n_missing,expected_small_family_pairs=len(expected)))
        fixture=ROOT/'results/orthology/native-small-family-output-fixture-v2'
        source=fixture/'Source/WorkingDirectory';result=fixture/'Results_small_fixture'
        expected,sizes,names,species,reverse=expected_pairs(source)
        clone=root/'corrupted';shutil.copytree(result,clone)
        path=clone/'Orthologues/Taxon0.tsv';original=path.read_text()
        smallrow=next(line for line in original.splitlines() if line.startswith('OG0000001\t'))
        path.write_text(original+smallrow+'\n')
        rejected(lambda:scan_native(clone,sizes,species,reverse,expected))
        path.write_text(original+'OG0000006\tTaxon2\tprotein_0_1\tprotein_2_1\n')
        rejected(lambda:scan_native(clone,sizes,species,reverse,expected))
        path.unlink()
        rejected(lambda:scan_native(clone,sizes,species,reverse,expected))
        binding=root/'binding';shutil.copytree(fixture,binding)
        command=[sys.executable,str(ROOT/'scripts/supplement_small_family_orthologs.py'),
                 '--source',str(binding/'Source/WorkingDirectory'),
                 '--results',str(binding/'Results_small_fixture'),
                 '--completion',str(binding/'receipt.json'),'--output',str(root/'rejected')]
        ids=binding/'Source/WorkingDirectory/SequenceIDs.txt'
        ids.write_text(ids.read_text().replace('protein_0_0','changed_0_0'))
        failed=subprocess.run(command,capture_output=True,text=True)
        assert failed.returncode and 'Fixture source changed' in failed.stderr
        receipt_path=binding/'receipt.json';r=json.loads(receipt_path.read_text())
        r['status']='running_native_reconciliation';receipt_path.write_text(json.dumps(r))
        failed=subprocess.run(command,capture_output=True,text=True)
        assert failed.returncode and 'Native execution is not complete' in failed.stderr
        print(json.dumps(dict(status='passed_small_family_supplement_fixture_checks',fixtures=checked,
              corruption_checks=['duplicate_pair','unexpected_singleton_pair','missing_taxon_table','changed_bound_source','incomplete_native_execution']),indent=2))


if __name__=='__main__':main()
