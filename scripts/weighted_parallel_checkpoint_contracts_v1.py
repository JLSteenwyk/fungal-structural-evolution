"""Require complete and semantically consistent parallel numerical checkpoints."""
import json
from ancestral_chain_attempt import sha


def checkpoint_pair(primary, independent, entry, index):
    for record in [primary, independent]:
        assert record['cohort_index']==index and record['cohort_id']==entry['cohort_id']
        assert record['entry']==entry and record['audits']==entry['audits']==1200
        assert record['links']==entry['links'] and record['designs']==30
        assert record['counts']==entry['audit_status_counts']
        assert record['link_counts']==entry['link_status_counts']
        assert sum(record['counts'].values())==record['audits']
        assert sum(record['link_counts'].values())==record['links']
        assert sum(record['policy_counts'].values())==record['audits']
        assert record['cached_numeric_inputs_preserved'] is True
        assert record['cached_input_array_bindings_checked']>0
        assert record['worker_address_space_limit_bytes']==12*2**30
        assert record['scientific_eligibility'] is False
        assert record['worker']['pid']>0 and record['worker']['created']>0
        assert isinstance(record['worker']['cmdline'],list) and record['worker']['cmdline']
    for field in ['counts','link_counts','policy_counts','audits','links','designs']:
        assert primary[field]==independent[field]


def checkpoint_artifacts(root,manifest,producer,reader):
    names={'checkpoints/'+str(index).zfill(5)+'.producer.json' for index in range(len(manifest))}
    rnames={'checkpoints/'+str(index).zfill(5)+'.reader.json' for index in range(len(manifest))}
    assert set(reader['reader_checkpoint_artifacts'])==rnames
    assert producer['cached_numeric_inputs_preserved_for_all_cohorts'] is reader['cached_numeric_inputs_preserved_for_all_cohorts'] is True
    assert producer['parallel_workers']==reader['parallel_workers']
    assert producer['worker_address_space_gib']==reader['worker_address_space_gib']==12
    for index,entry in enumerate(manifest):
        p='checkpoints/'+str(index).zfill(5)+'.producer.json'
        r='checkpoints/'+str(index).zfill(5)+'.reader.json'
        assert producer['artifacts'][p]==sha(root/p)
        assert reader['reader_checkpoint_artifacts'][r]==sha(root/r)
        checkpoint_pair(json.loads((root/p).read_text()),json.loads((root/r).read_text()),entry,index)
    return names
