#!/usr/bin/env python3
"""Bind all unchanged qualified encodings to the audited complete model mapping."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import time
import psutil
from merge_qualified_structure_encodings import rows, fingerprints
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args()
    plan=json.loads(a.plan.read_text());plan_sha=sha(a.plan)
    def verify():
        if sha(a.plan)!=plan_sha:
            raise ValueError('Plan changed')
        for name,digest in plan['pins'].items():
            if sha(ROOT/name)!=digest:
                raise ValueError('Changed dependency: '+name)
    verify()
    out=ROOT/plan['output'];control=ROOT/plan['control']
    if out.exists() or control.exists():
        raise FileExistsError('Use fresh immutable outputs')
    control.mkdir(parents=True)
    start=time.time()
    def state(status,**fields):
        p=control/'state.tmp'
        p.write_text(json.dumps(dict(status=status,elapsed_seconds=time.time()-start,**fields),indent=2)+'\n')
        p.replace(control/'state.json')
    state('waiting_for_exact_predecessors')
    for wait in plan['predecessors']:
        while True:
            try:
                process=psutil.Process(wait['pid'])
                live=process.create_time()==wait['create_time'] and process.status()!=psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:
                live=False
            if not live:
                break
            time.sleep(20)
        receipt=json.loads((ROOT/wait['receipt']).read_text())
        if receipt['status']!=wait['status'] or receipt['plan_sha256']!=sha(ROOT/wait['plan']):
            raise ValueError('Incomplete or changed predecessor')
    verify()
    if shutil.disk_usage(ROOT).free<plan['resources']['minimum_free_disk_gib']*2**30:
        raise ValueError('Insufficient disk')
    snapshot=ROOT/plan['mapping'];inventory=ROOT/plan['inventory']
    mr=checked_receipt(snapshot);ir=checked_receipt(inventory)
    mapped_audit=checked_receipt(ROOT/plan['mapping_readback'])
    if (mapped_audit['status']!='passed_full_local_residue_mapping_readback'
            or mapped_audit['mapping_receipt_sha256']!=sha(snapshot/'receipt.json')
            or mapped_audit['models']!=plan['models']
            or mr['source_inventory_sha256']!=sha(inventory/'inventory.jsonl')
            or ir['status']!='complete_disjoint_prediction_inventory_union_with_batch_provenance'):
        raise ValueError('Verified complete mapping required')
    inventory_rows=[json.loads(line) for line in (inventory/'inventory.jsonl').open()]
    expected_models={row['models'][0]['model_id']:dict(row['models'][0],source_record_id=row['record_id']) for row in inventory_rows}
    union_provenance={}; summaries={};sequences=set();totals=Counter();sources=[]
    tables={name:Counter() for name in ['marker_structure_links.tsv','matrix_to_structure_residues.tsv.gz']}
    dynamic_pins={}
    for cohort in plan['cohorts']:
        state('checking_qualified_cohort',cohort=cohort['name'])
        mapping=ROOT/cohort['mapping'];encoding=ROOT/cohort['encodings']
        source_mr=checked_receipt(mapping);er=checked_receipt(encoding)
        mapping_sha=sha(mapping/'receipt.json');encoding_sha=sha(encoding/'receipt.json')
        confidence_path=ROOT/cohort['confidence_readback']
        confidence=json.loads(confidence_path.read_text())
        if (source_mr['matrix_receipt_sha256']!=mr['matrix_receipt_sha256']
                or er['status']!='complete_native_3di_feature_audit' or er['confidence_stage']!='plddt_and_pae'
                or er['mapping_receipt_sha256']!=mapping_sha or er['models']!=cohort['models']
                or confidence['status']!='passed_full_original_pae_context_readback'
                or confidence['models']!=cohort['models']
                or confidence['source_receipts']['snapshot']!=mapping_sha
                or confidence['source_receipts']['qualified']!=encoding_sha):
            raise ValueError('Missing matched full confidence qualification')
        provenance=json.loads((mapping/'model_provenance.json').read_text())
        by_name={Path(m['path']).stem:m for m in provenance}
        if len(by_name)!=cohort['models'] or len(provenance)!=cohort['models']:
            raise ValueError('Source provenance scope differs')
        for model in provenance:
            mid=model['model_id']
            if mid in union_provenance or model!=expected_models.get(mid):
                raise ValueError('Overlapping or changed source provenance')
            union_provenance[mid]=model
        found=set();local_totals=Counter()
        for row in rows(encoding/'model_summary.tsv'):
            name=row['model_name'];sid=row['sequence_sha256']
            if name in summaries or name in found or name not in by_name or sid in sequences:
                raise ValueError('Overlapping or unexpected encoding identity')
            if sid!=by_name[name]['sequence_sha256'] or sha(ROOT/row['encoding_path'])!=row['encoding_sha256']:
                raise ValueError('Changed encoding or sequence')
            summaries[name]=row;found.add(name);sequences.add(sid)
            for key in er['totals']:
                local_totals[key]+=int(row[key])
        if found!=set(by_name) or dict(local_totals)!=er['totals']:
            raise ValueError('Encoding count or summary totals differ')
        totals.update(local_totals)
        for filename in tables:
            tables[filename].update(fingerprints(mapping/filename))
        for path in [mapping/'receipt.json',encoding/'receipt.json',confidence_path]:
            dynamic_pins[str(path)]=sha(path)
        sources.append(dict(name=cohort['name'],mapping=cohort['mapping'],mapping_receipt_sha256=mapping_sha,
                            encodings=cohort['encodings'],encoding_receipt_sha256=encoding_sha,models=len(found),
                            confidence_readback_sha256=sha(confidence_path),native_config_sha256=er['native_config_sha256'],
                            coordinate_audit_receipt_sha256=er['coordinate_audit_receipt_sha256'],pae_receipt_sha256=er['pae_receipt_sha256']))
    state('checking_complete_mapping_union')
    combined_provenance=json.loads((snapshot/'model_provenance.json').read_text())
    if (len(combined_provenance)!=plan['models'] or len(summaries)!=plan['models']
            or {m['model_id']:m for m in combined_provenance}!=union_provenance
            or set(union_provenance)!=set(expected_models)):
        raise ValueError('Combined model/provenance scope differs')
    table_counts={}
    for name,expected in tables.items():
        actual=fingerprints(snapshot/name)
        if actual!=expected or any(n!=1 for n in actual.values()):
            raise ValueError('Combined mapping differs from disjoint cohort union')
        table_counts[name]=sum(actual.values())
    out.mkdir()
    path=out/'model_summary.tsv'
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(next(iter(summaries.values()))),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(summaries[k] for k in sorted(summaries))
    verify()
    for path,digest in dynamic_pins.items():
        if sha(Path(path))!=digest:
            raise ValueError('Source receipt changed during integration')
    result=dict(status='complete_union_of_previously_qualified_encodings',models=len(summaries),totals=dict(totals),
                confidence_stage='plddt_and_pae',mapping_receipt_sha256=sha(snapshot/'receipt.json'),
                inventory_receipt_sha256=sha(inventory/'receipt.json'),source_cohorts=sources,
                verified_mapping_union_rows=table_counts,plan_sha256=plan_sha,script_sha256=sha(Path(__file__)),
                artifacts={'model_summary.tsv':sha(out/'model_summary.tsv')},
                interpretation='Unchanged encodings, complete disjoint source mapping union, and matched full directional confidence readbacks. Seven original prediction configurations retained by the source inventory. No new prediction, coordinate averaging, homology claim or evolutionary test.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    state(result['status'])


if __name__=='__main__':
    main()
