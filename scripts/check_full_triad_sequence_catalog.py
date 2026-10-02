#!/usr/bin/env python3
"""Check source dispositions, physical reuse and changed sequence catalog exports."""
import argparse
import copy
import json
from pathlib import Path
import prepare_full_triad_sequence_catalog as producer
import readback_full_triad_sequence_catalog as reader
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output;out.mkdir(parents=True,exist_ok=False)
    catalog=[dict(triad_id='first',models=[['A',1],['B',1],['C',1]],source_design_ready_links=2),
             dict(triad_id='swapped',models=[['B',1],['A',1],['C',1]],source_design_ready_links=1),
             dict(triad_id='unmodeled',models=[['A',1],['B',1],['D',1]],source_design_ready_links=0)]
    inputs={('A',1,'full'):dict(sequence='ACDE',original_length=4),
            ('B',1,'full'):dict(sequence='ACE',original_length=3),
            ('C',1,'full'):dict(sequence='AACDE',original_length=5)}
    mapping=dict(target_contexts=9,reference_tie_records=7,duplicate_reference_links=14)
    def sources(*args):return copy.deepcopy(catalog),catalog[:2],copy.deepcopy(inputs),mapping,{}
    producer.load_original=sources;reader.load_original=sources
    plan=out/'plan.json';plan.write_text(json.dumps(dict(output=str(out/'catalog'),resources=dict(minimum_free_disk_gib=0),scope='software fixture')))
    producer.run(plan);reader.inspect(plan,out/'readback.json')
    r=json.loads((out/'readback.json').read_text());assert r['source_ready_triads']==2 and r['unique_sequence_model_sets']==1 and r['native_alignment_states']==12
    root=out/'catalog';rejected=0
    mutations=[('sequence_sets.jsonl','sequences'),('sequence_sets.jsonl','sequence_sha256'),
        ('sequence_sets.jsonl','original_lengths'),('sequence_sets.jsonl','models'),
        ('sequence_sets.jsonl','sequence_set_id'),('triad_sequence_links.jsonl','role_to_sequence_indices'),
        ('triad_sequence_links.jsonl','triad_id'),('triad_sequence_links.jsonl','sequence_control_disposition')]
    rp=root/'receipt.json';saved_receipt=rp.read_bytes()
    for name,field in mutations:
        path=root/name;original=path.read_bytes();rows=[json.loads(l) for l in path.open()]
        rows[0][field]=['altered'] if isinstance(rows[0][field],list) else 'altered'
        path.write_text(''.join(json.dumps(x)+'\n' for x in rows))
        receipt=json.loads(saved_receipt);receipt['artifacts'][name]=sha(path);rp.write_text(json.dumps(receipt))
        try:reader.inspect(plan,out/('changed-'+str(rejected)+'.json'))
        except (AssertionError,KeyError,TypeError,IndexError):rejected+=1
        else:raise AssertionError('Changed sequence catalog accepted '+field)
        finally:path.write_bytes(original);rp.write_bytes(saved_receipt)
    # Missing or additional original triads must fail even after rehashing.
    for extra in [False,True]:
        path=root/'triad_sequence_links.jsonl';original=path.read_bytes();lines=path.read_text().splitlines()
        lines=(lines+[lines[0]]) if extra else lines[:-1];path.write_text('\n'.join(lines)+'\n')
        receipt=json.loads(saved_receipt);receipt['artifacts'][path.name]=sha(path);rp.write_text(json.dumps(receipt))
        try:reader.inspect(plan,out/('changed-'+str(rejected)+'.json'))
        except AssertionError:rejected+=1
        else:raise AssertionError('Changed triad coverage accepted')
        finally:path.write_bytes(original);rp.write_bytes(saved_receipt)
    assert rejected==10
    files=['scripts/full_triad_sequence_sources.py','scripts/prepare_full_triad_sequence_catalog.py',
           'scripts/readback_full_triad_sequence_catalog.py',str(Path(__file__))]
    result=dict(status='passed_full_triad_sequence_catalog_software_contracts',altered_exports_rejected=rejected,
        checked_reused_physical_set=True,checked_unscheduled_original_triads=True,scientific_eligibility=False,
        source_hashes={s:sha(s) for s in files})
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
