#!/usr/bin/env python3
"""Export all completed categorical groups and every certified numerical review."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import gzip
import json
from pathlib import Path

from ancestral_chain_attempt import sha


def write_table(path,rows):
    assert rows and all(set(r)==set(rows[0]) for r in rows)
    with path.open('x',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--table',type=Path,required=True)
    parser.add_argument('--reviews',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();path=Path('metadata/independent_baliphy_category_v2_completed_20261002.json')
    closed=json.loads(path.read_text());assert closed['status']=='complete_verified_full_independent_baliphy_categorical_comparison'
    archive=Path(closed['full_hash_archive']);assert sha(archive)==closed['full_hash_archive_sha256']
    proof=json.loads(archive.read_text());assert len(proof['services'])==2
    root=archive.parent;chains_path=Path('results/ancestral/full-baliphy-horizon-resources-20261003-v1/chains.jsonl')
    chain_rows=[json.loads(x) for x in chains_path.read_text().splitlines()];lookup={}
    for chain in chain_rows:
        group=chain['model_input_identity'];value=(chain['family'],chain['prior_label'],chain['effective_input_group'],
            json.dumps(chain['original_configuration_ids'],separators=(',', ':')))
        assert group not in lookup or lookup[group]==value;lookup[group]=value
    census_path=Path('metadata/baliphy_horizon_resources_completed_20261003.json');census=json.loads(census_path.read_text())
    census_archive=Path(census['full_hash_archive']);assert sha(census_archive)==census['full_hash_archive_sha256']
    assert sha(chains_path)==json.loads(census_archive.read_text())['source_hashes'][str(chains_path)]
    bindings={str(p):sha(p) for p in [Path(__file__),path,archive,chains_path,census_path,census_archive]}
    tables=[];reviews=[];numeric=Counter();pattern=Counter();coordinates=Counter();group_count=0;failed=0;indicators=0
    for checkpoint in sorted((root/'groups').glob('*.json')):
        assert sha(checkpoint)==proof['source_hashes'][str(checkpoint)];bindings[str(checkpoint)]=sha(checkpoint)
        document=json.loads(checkpoint.read_text());group=document['group'];result=document['result']
        assert group==checkpoint.stem and result['scientific_eligibility'] is False
        family,prior,input_group,aliases=lookup[group];group_count+=1
        is_failed=result['status']=='unresolved_failed_native_chain_retained';failed+=is_failed
        assert set(result['cutoffs'])==set() if is_failed else set(result['cutoffs'])=={'250','500'}
        for cut in ['250','500']:
            value=result['cutoffs'].get(cut)
            row=dict(group=group,family=family,prior=prior,effective_input_group=input_group,original_configuration_aliases_json=aliases,
                discard_through=cut,group_status=result['status'],patterns=None,coordinates=None,declared_indicators=None,
                no_variation_patterns=None,mixing_review_patterns=None,no_variation_coordinates=None,mixing_review_coordinates=None,
                defined_indicator_rows_compared=None,undefined_indicator_rows=None,binary_zero_pair_indicator_rows=None,
                singular_indicator_rows=None,posterior_qualified=False)
            if value is not None:
                row.update(patterns=value['patterns'],coordinates=value['coordinates'],declared_indicators=value['declared_indicator_rows'],
                    no_variation_patterns=value['pattern_status_counts'].get('no_observed_state_variation_requires_review',0),
                    mixing_review_patterns=value['pattern_status_counts'].get('categorical_mixing_requires_review',0),
                    no_variation_coordinates=value['coordinate_status_counts'].get('no_observed_state_variation_requires_review',0),
                    mixing_review_coordinates=value['coordinate_status_counts'].get('categorical_mixing_requires_review',0),
                    defined_indicator_rows_compared=value['numeric_comparison_status_counts'].get('all_defined_metrics_compared',0),
                    undefined_indicator_rows=value['numeric_comparison_status_counts'].get('original_disposition_only_no_defined_metrics',0),
                    binary_zero_pair_indicator_rows=value['binary_zero_pair_indicator_rows'],singular_indicator_rows=value['singular_indicator_rows'])
                numeric.update(value['numeric_comparison_status_counts']);indicators+=value['declared_indicator_rows']
                pattern.update({cut+':'+k:v for k,v in value['pattern_status_counts'].items()})
                coordinates.update({cut+':'+k:v for k,v in value['coordinate_status_counts'].items()})
                if value['binary_zero_pair_indicator_rows']:
                    comparison=root/value['serialized_comparison'];assert sha(comparison)==value['serialized_comparison_sha256']
                    bindings[str(comparison)]=sha(comparison);found=set()
                    with gzip.open(comparison,'rt') as handle:
                        for line in handle:
                            record=json.loads(line)
                            for flag in record['unresolved_numeric_metrics']:
                                assert flag['reason']=='exact_binary_zero_pair_truncation_requires_review'
                                reviews.append(dict(group=group,family=family,prior=prior,discard_through=cut,
                                    pattern_id=record['pattern_id'],coordinate_multiplicity=record['coordinate_multiplicity'],
                                    state=flag['state'],metric=flag['metric'],original_value=flag['original'],
                                    independent_value=flag['independent'],reason=flag['reason'],
                                    admissible_roundoff_outcomes_json=json.dumps(flag['admissible_roundoff_outcomes']),
                                    certificates_json=json.dumps(flag['certificates'],sort_keys=True,separators=(',', ':')),
                                    serialized_comparison=str(comparison),serialized_comparison_sha256=sha(comparison),
                                    posterior_qualified=False));found.add((record['pattern_id'],flag['state']))
                    assert len(found)==value['binary_zero_pair_indicator_rows']
            tables.append(row)
    assert group_count==len(lookup)==405 and failed==2 and len(tables)==810
    indicator_keys={(r['group'],r['discard_through'],r['pattern_id'],r['state']) for r in reviews}
    assert len(indicator_keys)==closed['binary_zero_pair_indicator_rows']==15
    assert indicators==closed['declared_indicator_rows'] and sum(r['patterns'] or 0 for r in tables)==closed['pattern_cutoff_rows']
    assert dict(numeric)==closed['numeric_comparison_status_counts'] and dict(pattern)==closed['pattern_status_counts']
    assert dict(coordinates)==closed['coordinate_status_counts']
    assert Counter(key[1] for key in indicator_keys)=={'250':2,'500':13}
    write_table(args.table,tables);write_table(args.reviews,reviews)
    result=dict(status='exported_full_closed_ancestral_category_review_tables',checked_utc=datetime.now(timezone.utc).isoformat(),
        model_quartets=405,cutoff_rows=810,unresolved_model_quartets=2,explicit_missing_cutoff_rows=4,
        declared_indicator_rows=indicators,certified_numerical_review_indicator_rows=15,
        certified_numerical_review_metric_rows=len(reviews),review_indicator_rows_by_cutoff={'250':2,'500':13},
        source_hashes=bindings,artifacts={str(p):sha(p) for p in [args.table,args.reviews]},scientific_eligibility=False,
        scope='All405closedmodelgroups and bothcutoffs exported, including failedquartets with blank measurements. Every metric flag for all15certifiedreview indicators extracted from every comparison file whose full closed summary contains a review; complete production already checked every indicator. Multiple metric flags for one indicator remain separate rows. No numerical agreement claimed for those reviews; multiplicities/cutoffs overlap and are not independent biological counts. No posterior or model acceptance.')
    with args.receipt.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True)


if __name__=='__main__':main()
