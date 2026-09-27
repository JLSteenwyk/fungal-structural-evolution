#!/usr/bin/env python3
"""Independently compare every combined table field with the selected source audit rows."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def rows(path):
    with Path(path).open() as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        return reader.fieldnames,list(reader)


def verify_collection(collection,inventory_path):
    collection=Path(collection);inventory_path=Path(inventory_path)
    receipt=json.loads((collection/'receipt.json').read_text())
    inventory=json.loads(inventory_path.read_text())
    if receipt['status']!='complete_source_preserving_paired_fit_table_collection':raise ValueError('Incomplete collection')
    if receipt['source_sha256'][str(inventory_path)]!=sha(inventory_path):raise ValueError('Inventory binding differs')
    for path,h in receipt['source_sha256'].items():
        if sha(path)!=h:raise ValueError('Changed source')
    for name,h in receipt['artifacts'].items():
        if sha(collection/name)!=h:raise ValueError('Changed collection artifact')
    selected={r['marker']:r for r in inventory['sources']}
    if len(selected)!=len(inventory['sources']) or len(selected)!=receipt['markers']:raise ValueError('Marker universe differs')
    proofs={}
    for name in ['fit_summary.tsv','paired_branches.tsv','warnings.tsv']:
        # Audit locations are read from producer provenance, then their exact
        # selected fit-directory memberships are independently checked below.
        output=collection/name
        actual_fields,actual=rows(output) if output.exists() else (None,[])
        audits={Path(p).parent for p in receipt['source_sha256']
                if Path(p).name=='receipt.json' and Path(p).parent.name.startswith('paired-')
                and json.loads(Path(p).read_text()).get('status')=='complete_paired_fit_audit'}
        if len(audits)!=2:raise ValueError('Expected two audited sources')
        expected=[];fields=None
        for auditdir in sorted(audits):
            ar=json.loads((auditdir/'receipt.json').read_text())
            sourcefile=auditdir/name
            if not sourcefile.exists():
                if name!='warnings.tsv':raise ValueError('Missing required source table')
                continue
            if sha(sourcefile)!=ar['artifacts'][name]:raise ValueError('Changed source audit table')
            header,source_rows=rows(sourcefile)
            if fields is not None and fields!=header:raise ValueError('Source schemas differ')
            fields=header
            # The fit receipt hash determines which native run this audit represents.
            source_roots={Path(r['fit_directory']).parent for r in selected.values()}
            matched=[p for p in source_roots if sha(p/'receipt.json')==ar['fit_receipt_sha256']]
            if len(matched)!=1:raise ValueError('Audit fit binding ambiguous')
            root=matched[0]
            for row in source_rows:
                marker=row['marker'];source=selected.get(marker)
                if source is None or Path(source['fit_directory']).parent!=root:continue
                expected.append(dict(row,source_run=source['source'],source_fit_directory=str(root/marker),source_audit_directory=str(auditdir)))
        allfields=(fields or [])+['source_run','source_fit_directory','source_audit_directory']
        if actual and actual_fields!=allfields:raise ValueError('Combined schema differs')
        encode=lambda records:Counter(tuple(r[k] for k in allfields) for r in records)
        if encode(actual)!=encode(expected):raise ValueError('Combined field or row multiplicity differs: '+name)
        proofs[name]=dict(rows=len(actual),sha256=sha(output) if output.exists() else None)
    if proofs['fit_summary.tsv']['rows']!=4*len(selected) or proofs['paired_branches.tsv']['rows']!=receipt['paired_branches'] or proofs['warnings.tsv']['rows']!=receipt['warning_rows']:raise ValueError('Combined totals differ')
    return dict(status='passed_complete_combined_paired_fit_source_readback',producer_receipt_sha256=sha(collection/'receipt.json'),markers=len(selected),tables=proofs,scope='Every combined table field and row multiplicity independently compared with the correct selected source-audit rows and provenance. Source receipts/artifacts and exact marker source choices checked. Does not rerun native fits or replace their numeric audits; no uncertainty or biological acceleration inference.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--collection',type=Path,required=True);p.add_argument('--inventory',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=verify_collection(a.collection,a.inventory)
    with a.output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')


if __name__=='__main__':main()
