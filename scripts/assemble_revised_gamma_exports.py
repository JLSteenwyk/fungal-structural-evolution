#!/usr/bin/env python3
"""Assemble an immutable full rate cohort with one explicitly reviewed replacement."""
import argparse
import json
from pathlib import Path
import shutil
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import read_table, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['source','source-audit','refits','readback','plan','output']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--marker', required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    source = json.loads((a.source/'receipt.json').read_text())
    audit = checked_receipt(a.source_audit); refits = checked_receipt(a.refits)
    readback = json.loads(a.readback.read_text()); plan = json.loads(a.plan.read_text())
    if (audit['rate_receipt_sha256'] != sha(a.source/'receipt.json')
        or audit['status'] != 'passed_full_site_rate_output_audit'
        or refits['plan_sha256'] != sha(a.plan)
        or readback['receipt_sha256'] != sha(a.refits/'receipt.json')
        or readback['status'] != 'passed_independent_gamma_refit_rate_topology_readback'):
        raise ValueError('Reviewed source/refit lineage required')
    choices = read_table(a.refits/'summary.tsv')
    chosen = max(choices, key=lambda r: float(r['log_likelihood']))
    if float(chosen['minus_export']) <= 0 or a.marker not in Path(plan['alignment']).parts:
        raise ValueError('Replacement does not improve the specified marker')
    config = json.loads((a.source/'config.json').read_text())
    if sha(a.source/'config.json') != source['config_sha256'] or config['heterogeneity'] != 'G4':
        raise ValueError('Changed source or non-Gamma source')
    config['assembly_script_sha256'] = sha(Path(__file__))
    config['revision'] = dict(marker=a.marker, fit='aa', chosen_start=chosen['start'],
        source_receipt_sha256=sha(a.source/'receipt.json'), refit_receipt_sha256=sha(a.refits/'receipt.json'),
        readback_sha256=sha(a.readback), interpretation='Only this AA fit replaced by highest observed reviewed Gamma likelihood; all other numerical outputs copied byte-for-byte. No global-optimum claim.')
    for path in [a.source/'receipt.json', a.refits/'receipt.json', a.readback, a.plan]:
        config['pinned_files'][str(path)] = sha(path)
    a.output.mkdir(parents=True)
    (a.output/'config.json').write_text(json.dumps(config,indent=2)+'\n'); config_hash=sha(a.output/'config.json')
    receipts = {}; replacement_count = 0
    for key, digest in sorted(source['fit_receipts'].items()):
        marker, label = key.split('/'); folder=a.source/marker; dest=a.output/marker; dest.mkdir(exist_ok=True)
        rp=folder/(label+'.receipt.json'); r=json.loads(rp.read_text())
        if sha(rp)!=digest:
            raise ValueError('Source fit receipt changed')
        for name,h in r['artifacts'].items():
            if sha(folder/name)!=h:
                raise ValueError('Source fit artifact changed')
        request=json.loads((folder/(label+'.config.json')).read_text())
        replaced=marker==a.marker and label=='aa'
        if replaced:
            replacement_count+=1
            request['command']=next(x['command'] for x in plan['requests'] if x['name']==chosen['start'])
            start=Path(request['command'][request['command'].index('-te')+1])
            request['reference_topology_sha256']=request['topology_sha256']
            request['topology_sha256']=sha(start)
            request['reviewed_refit_receipt_sha256']=sha(a.refits/'receipt.json')
        request['parent_config_sha256']=config_hash
        request['source_request_sha256']=sha(folder/(label+'.config.json'))
        for name in r['artifacts']:
            if name==label+'.config.json':
                continue
            src=(a.refits/(chosen['start']+name[len(label):])) if replaced else folder/name
            shutil.copyfile(src,dest/name)
        (dest/(label+'.config.json')).write_text(json.dumps(request,indent=2)+'\n')
        r['parent_config_sha256']=config_hash
        r['source_receipt_sha256']=digest
        r['assembly_action']='reviewed_refit_replacement' if replaced else 'unchanged_numerical_outputs'
        if replaced:
            r['elapsed_seconds']=None
            r['reviewed_refit_receipt_sha256']=sha(a.refits/'receipt.json')
        r['artifacts']={name:sha(dest/name) for name in r['artifacts']}
        (dest/(label+'.receipt.json')).write_text(json.dumps(r,indent=2)+'\n')
        receipts[key]=sha(dest/(label+'.receipt.json'))
    if replacement_count!=1:
        raise ValueError('Expected exactly one replacement')
    result=dict(source, config_sha256=config_hash, fit_receipts=receipts,
        revision=config['revision'], interpretation=config['revision']['interpretation'])
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Assembled',len(receipts),'fits with one reviewed replacement')


if __name__=='__main__':
    main()
