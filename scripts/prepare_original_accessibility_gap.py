#!/usr/bin/env python3
"""Derive the exact unannotated remainder of the completed original prediction cohort."""
import argparse
import json
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from catalog_whole_proteome_structures import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ['full','previous','previous-audit','output']:ap.add_argument('--'+k,type=Path,required=True)
    a=ap.parse_args();full=checked_receipt(a.full);previous=checked_receipt(a.previous);audit=checked_receipt(a.previous_audit)
    if audit['status']!='passed_full_accessibility_snapshot' or audit['snapshot_receipt_sha256']!=sha(a.previous/'receipt.json'):raise ValueError('Prior audit not bound to previous snapshot')
    old=json.loads((a.previous/'model_provenance.json').read_text());complete=json.loads((a.full/'model_provenance.json').read_text())
    before={m['model_id']:m for m in old};after={m['model_id']:m for m in complete}
    if len(before)!=len(old) or len(after)!=len(complete) or not before.keys()<=after.keys():raise ValueError('Invalid model grid')
    if audit['models_audited']!=len(old) or audit['residues_audited']!=sum(m['length'] for m in old):raise ValueError('Prior audit scope differs')
    relocated=0
    for key,m in before.items():
        target=after[key]
        if {k:v for k,v in m.items() if k!='path'}!={k:v for k,v in target.items() if k!='path'}:raise ValueError('Changed common model identity or coordinate checksum')
        relocated+=m['path']!=target['path']
    missing=[m for m in complete if m['model_id'] not in before]
    if not missing:raise ValueError('No remaining cohort')
    a.output.mkdir(parents=True,exist_ok=False);target=a.output/'model_provenance.json';target.write_text(json.dumps(missing,indent=2)+'\n')
    r={'status':'complete_exact_accessibility_gap_model_manifest','models':len(missing),'residues':sum(m['length'] for m in missing),'full_models':len(complete),'previously_audited_models':len(old),'common_models_with_relocated_paths':relocated,'source_receipts':{str(p/'receipt.json'):sha(p/'receipt.json') for p in [a.full,a.previous,a.previous_audit]},'script_sha256':sha(__file__),'artifacts':{target.name:sha(target)},'scope':'Exact model-ID complement, requiring unchanged common provenance fields and coordinate checksums except relocation of paths. Input for missing accessibility only, not a marker mapping or a new coordinate-byte audit. Production ASA and its source audit must verify actual coordinate bytes.'}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
