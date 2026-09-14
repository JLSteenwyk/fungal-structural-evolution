#!/usr/bin/env python3
"""Replace affected-marker frame rows after a completed FCS rate workflow."""
import argparse
import json
import math
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['baseline', 'changed', 'inputs', 'accessibility', 'controller', 'output']:
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipts = {k: checked_receipt(getattr(a,k)) for k in ['baseline','changed','inputs','accessibility']}
    br, cr, ir, ar = (receipts[k] for k in ['baseline','changed','inputs','accessibility'])
    if any(r['status'] != 'complete_site_rate_exposure_analysis_frame' for r in [br,cr]):
        raise ValueError('Incomplete source frames')
    if br['source_receipts']['inputs'] != ir['source_input_receipt_sha256'] or cr['source_receipts']['inputs'] != sha(a.inputs/'receipt.json'):
        raise ValueError('Frame input lineage differs')
    if ar['source_receipts']['inputs'] != sha(a.inputs/'receipt.json') or ar['source_receipts']['baseline_inputs'] != ir['source_input_receipt_sha256']:
        raise ValueError('Full accessibility input lineage differs')
    controller = json.loads(a.controller.read_text())
    if controller['status'] != 'complete_fcs_rate_and_exposure_frame_handoff':
        raise ValueError('Incomplete controller')
    # Require the completed frame's exact checkpoint, and every upstream stage.
    stages = controller['stages']
    for stage in stages:
        path = Path(stage['receipt'])
        if sha(path) != stage['receipt_sha256']:
            raise ValueError('Controller stage changed')
    if len(stages) != 11 or not any(stage['receipt_sha256'] == sha(a.changed/'receipt.json') for stage in stages):
        raise ValueError('Complete eleven-stage frame handoff required')
    disposition = read_table(a.inputs/'marker_disposition.tsv')
    ready = {r['marker']:r for r in disposition if r['baseline_status']=='ready_for_inference'}
    changed = {m for m,r in ready.items() if r['sensitivity_status']=='ready_for_inference'}
    unchanged = {m for m,r in ready.items() if r['sensitivity_status']=='unchanged_reuse_baseline_fit'}
    if changed|unchanged != set(ready) or len(changed)!=cr['markers'] or len(ready)!=br['markers']:
        raise ValueError('Unsupported marker universe')
    if any(int(r['all_missing_columns_removed']) for r in ready.values()):
        raise ValueError('Coordinate changes require separate handling')
    sources = {k: read_table(getattr(a,k)/'site_rate_exposure.tsv') for k in ['baseline','changed']}
    key = lambda r:(r['marker'],r['paired_column_1based'])
    indexed = {k:{key(r):r for r in rows} for k,rows in sources.items()}
    if any(len(indexed[k])!=len(sources[k]) for k in sources):
        raise ValueError('Duplicate site identities')
    if set(indexed['changed']) != {k for k in indexed['baseline'] if k[0] in changed}:
        raise ValueError('Changed site grid differs')
    output=[]; provenance=[]
    for row in sources['baseline']:
        marker=row['marker'];replacement=indexed['changed'].get(key(row),row)
        if marker not in ready or set(row)!=set(replacement) or row['matrix_column_1based']!=replacement['matrix_column_1based']:
            raise ValueError('Unexpected frame schema or coordinate')
        if int(replacement['marker_taxa'])!=int(ready[marker]['sensitivity_eligible_taxa']):
            raise ValueError('Retained taxon count differs')
        rates=[v for k,v in replacement.items() if k.endswith('_gamma_rate') or k.endswith('_freerate_rate')]
        if len(rates)!=8 or any(not math.isfinite(float(v)) or float(v)<0 for v in rates):
            raise ValueError('Incomplete valid rate grid')
        output.append(replacement)
    counts={m:sum(int(r['observed_taxa']) for r in output if r['marker']==m) for m in ready}
    if counts != ar['marker_row_counts'] or sum(counts.values())!=ar['rows']:
        raise ValueError('Full exposure observation grid differs')
    diagnostics=[]
    for label in ['baseline','changed']:
        diagnostics.extend(r for r in read_table(getattr(a,label)/'fit_diagnostics.tsv')
                           if r['marker'] in (unchanged if label=='baseline' else changed))
    if len(diagnostics)!=4*len(ready) or len({(r['marker'],r['fit']) for r in diagnostics})!=len(diagnostics):
        raise ValueError('Fit diagnostic grid differs')
    for m in sorted(ready):
        provenance.append({'marker':m,'frame_source':'changed' if m in changed else 'baseline',
                           'source_receipt_sha256':sha((a.changed if m in changed else a.baseline)/'receipt.json')})
    a.output.mkdir(parents=True)
    write_table(a.output/'site_rate_exposure.tsv',output)
    write_table(a.output/'fit_diagnostics.tsv',diagnostics)
    write_table(a.output/'marker_provenance.tsv',provenance)
    result={'status':'complete_full_cohort_fcs_site_rate_exposure_frame','markers':len(ready),'sites':len(output),
            'site_rate_values':8*len(output),'taxon_site_observations':sum(counts.values()),
            'changed_markers':len(changed),'unchanged_markers':len(unchanged),
            'source_receipts':{k:sha(getattr(a,k)/'receipt.json') for k in receipts},
            'controller_receipt_sha256':sha(a.controller),'script_sha256':sha(Path(__file__)),
            'artifacts':{x.name:sha(x) for x in a.output.iterdir()},
            'interpretation':'Full baseline marker universe, with affected marker rows replaced only from the complete audited FCS workflow. Unchanged rows and all selected-row fields copied exactly. Counts match full retained normalized exposure. Not a completed coupling test; rates remain conditional on fitted trees/models and extant exposure is not ancestral or phylogenetically weighted.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
