#!/usr/bin/env python3
"""Join audited domain controls to complete duplicate/reference triads."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from prepare_duplication_sequence_covariates import sha

POLICIES = ('alignment_evalue', 'alignment_bitscore', 'envelope_evalue', 'envelope_bitscore')


def classify(classes, conservative):
    if 'identical_model' in classes:
        return 'identical_model_in_triad'
    if any(x in ('neither_annotated', 'one_unannotated') for x in classes):
        return 'incomplete_annotation'
    if any(x in ('different_annotation_content', 'same_content_different_order') for x in classes):
        return 'annotation_difference'
    if set(classes) != {'same_ordered_annotations'}:
        raise ValueError('Unknown architecture class')
    return 'conservative_same_ordered' if all(conservative) else 'same_ordered_not_conservative'


def rows(path):
    with Path(path).open() as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('controls', 'control-readback', 'references', 'reference-readback', 'queue', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args(); pins = {}
    def receipt(root, required):
        path = root / 'receipt.json'; r = json.loads(path.read_text())
        if r['status'] != required: raise ValueError('Wrong producer status')
        pins[str(path)] = sha(path)
        return r
    c = receipt(a.controls, 'complete_duplication_domain_control_inventory_pending_readback')
    r = receipt(a.references, 'complete_provisional_reference_comparison_inventory')
    q = receipt(a.queue, 'complete_reviewed_duplication_model_pair_queue')
    for path, root, status in [(a.control_readback, a.controls, 'passed_full_duplication_domain_control_readback'),
                               (a.reference_readback, a.references, 'passed_full_reference_comparison_ledger_readback')]:
        audit = json.loads(path.read_text())
        if audit['status'] != status or audit['producer_receipt_sha256'] != sha(root/'receipt.json'):
            raise ValueError('Readback not bound to producer')
        pins[str(path)] = sha(path)
    for root, rec, name in [(a.controls,c,'pair_architecture_controls.tsv'),
                            (a.references,r,'event_reference_comparisons.tsv'),
                            (a.queue,q,'event_model_pair_links.tsv')]:
        path = root / name
        if sha(path) != rec['artifacts'][name]: raise ValueError('Changed source table')
        pins[str(path)] = sha(path)
    controls = {}
    for row in rows(a.controls/'pair_architecture_controls.tsv'):
        key = row['pair_key'], row['policy']
        if key in controls or row['policy'] not in POLICIES: raise ValueError('Repeated/unknown policy')
        controls[key] = row
    if len(controls) != c['pair_policy_rows']: raise ValueError('Control scope mismatch')
    events = {}
    event_fields = ('guide','family','gene_node','gene_a','gene_b')
    for row in rows(a.queue/'event_model_pair_links.tsv'):
        key = tuple(row[k] for k in event_fields)
        if key in events: raise ValueError('Repeated event')
        events[key] = row
    triads = defaultdict(dict)
    for row in rows(a.references/'event_reference_comparisons.tsv'):
        key = tuple(row[k] for k in event_fields) + (row['reference_gene'],)
        side = row['focal_side']
        if side not in ('a','b') or side in triads[key]: raise ValueError('Repeated/bad focal side')
        triads[key][side] = row
    a.output.mkdir(parents=True, exist_ok=False)
    counts = Counter(); stability = Counter()
    fields = list(event_fields) + ['reference_gene','lexical_representative','model_a','version_a','model_b','version_b',
        'reference_model','reference_version','policy','duplicate_pair_key','a_reference_pair_key','b_reference_pair_key',
        'duplicate_annotation_class','a_reference_annotation_class','b_reference_annotation_class',
        'duplicate_both_conservative','a_reference_both_conservative','b_reference_both_conservative','triad_class']
    with (a.output/'triad_architecture.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
        for key, sides in sorted(triads.items()):
            if set(sides) != {'a','b'}: raise ValueError('Incomplete reference triad')
            left,right = sides['a'],sides['b']; event = events[key[:-1]]
            for k in ('reference_model','reference_version','lexical_representative'):
                if left[k] != right[k]: raise ValueError('Reference side disagreement')
            for side, row in sides.items():
                if (row['focal_gene'],row['focal_model'],row['focal_version']) != (event['gene_'+side],event['model_'+side],event['version_'+side]):
                    raise ValueError('Focal source mismatch')
            ids = [(event['model_a'],event['version_a']), (event['model_b'],event['version_b']), (left['reference_model'],left['reference_version'])]
            pairs = [event['pair_key'],left['pair_key'],right['pair_key']]
            equal = [ids[0]==ids[1],ids[0]==ids[2],ids[1]==ids[2]]
            statuses=[]
            for policy in POLICIES:
                annotations=[];flags=[]
                for pair, same in zip(pairs,equal):
                    control = None if same else controls[pair,policy]
                    annotations.append('identical_model' if same else control['annotation_class'])
                    flags.append('' if same else control['both_conservative'])
                status = classify(annotations,[x=='1' for x in flags]); statuses.append(status)
                record = dict(zip(event_fields,key[:-1]))
                record.update(reference_gene=key[-1],lexical_representative=left['lexical_representative'],
                    model_a=ids[0][0],version_a=ids[0][1],model_b=ids[1][0],version_b=ids[1][1],
                    reference_model=ids[2][0],reference_version=ids[2][1],policy=policy,triad_class=status)
                for label,pair,annotation,flag in zip(('duplicate','a_reference','b_reference'),pairs,annotations,flags):
                    record[label+'_pair_key']=pair;record[label+'_annotation_class']=annotation;record[label+'_both_conservative']=flag
                writer.writerow(record);counts[key[0]+':'+policy+':'+status]+=1
            stability['same_class_all_policies' if len(set(statuses))==1 else 'policy_sensitive_class']+=1
    for path,h in pins.items():
        if sha(path)!=h: raise ValueError('Source changed during join')
    expected = r['counts']['profile_event_reference_links'] + r['counts']['mafft_event_reference_links']
    if len(triads)!=expected: raise ValueError('Triad count mismatch')
    result=dict(status='complete_duplication_triad_architecture_controls',source_sha256=pins,
        script_sha256=sha(__file__),triads=len(triads),policy_rows=4*len(triads),counts=dict(counts),policy_stability=dict(stability),
        artifacts={'triad_architecture.tsv':sha(a.output/'triad_architecture.tsv')},
        scope='All provisional event/reference ties retained. Ordered Pfam annotation compatibility, not validated domain homology, gain/loss, PAE/orientation, structure divergence or biological orthology. Identical-model triads explicit; priority classes retain all individual pair fields. Policies/guides/ties are dependent sensitivity analyses.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','triads','policy_rows','counts','policy_stability')},indent=2))

if __name__=='__main__':main()
