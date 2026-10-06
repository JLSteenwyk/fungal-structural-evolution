#!/usr/bin/env python3
"""Bind the OG0000017 audit to a required downstream aggregate sensitivity."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', required=True, type=Path)
    parser.add_argument('--audit-check', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Fresh policy receipt required')
    audit = json.loads(args.audit.read_text())
    check = json.loads(args.audit_check.read_text())
    if audit['status'] != 'complete_og0000017_source_annotation_repeat_risk_audit':
        raise ValueError('Complete source-context audit required')
    if check['status'] != 'passed_og0000017_annotation_repeat_context_output_check':
        raise ValueError('Independent audit readback required')
    if audit['family'] != 'OG0000017' or audit['focal_taxon'] != 'F181123':
        raise ValueError('Unexpected family or focal taxon')
    if (audit['source_products'], audit['focal_family_products'], audit['other_focal_products']) != (133973, 46787, 87186):
        raise ValueError('Unexpected complete product census')
    cohorts = audit['cohorts']
    if (cohorts['family']['repeat_keyword_pfam_proteins'],
            cohorts['other']['repeat_keyword_pfam_proteins']) != (19764, 6046):
        raise ValueError('Unexpected descriptive Pfam-risk counts')
    if check['source_products'] != audit['source_products'] or check['cohorts'] != {'family': 46787, 'other': 87186}:
        raise ValueError('Audit/readback census mismatch')
    result = {
        'status': 'required_og0000017_aggregate_sensitivity_policy',
        'checked_utc': datetime.now(timezone.utc).isoformat(),
        'family': audit['family'], 'focal_taxon': audit['focal_taxon'],
        'audit_sha256': sha(args.audit), 'audit_check_sha256': sha(args.audit_check),
        'complete_source_products_retained': audit['source_products'],
        'family_products_retained': audit['focal_family_products'],
        'paired_aggregate_analyses': [
            'complete eligible family universe',
            'same eligible universe with OG0000017 contributions withheld only from aggregate summary'
        ],
        'invariants': [
            'Do not alter OG0000017 membership, sequences, annotations, family tree, or reconciliation inputs.',
            'Use identical eligibility, branch, confidence, missing-data, multiplicity, and phylogenetic controls in both aggregates.',
            'Report the family-withheld change as a dependency sensitivity, not as repeat classification or a biological duplication conclusion.',
            'Do not make family-specific duplication, homology, or structural-evolution claims without separate validated gene-tree and annotation evidence.'
        ],
        'scientific_eligibility': False,
        'scope': 'Binds the validated complete source-context audit to a predeclared downstream aggregate sensitivity. It creates no family exclusion, no tree/reconciliation change, and no biological classification.'
    }
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
