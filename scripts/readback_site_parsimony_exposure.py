#!/usr/bin/env python3
"""Independently rescore all reference-tree sites and read back exposure summaries."""
import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
from Bio import Phylo, SeqIO
from audit_site_parsimony_topologies import set_scores
from audit_joint_path_uncertainty import checked, rows, sha, quantile


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['diagnostic', 'inputs', 'fits', 'exposure', 'output']:
        ap.add_argument('--'+name, type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    r = checked(a.diagnostic); checked(a.inputs); er = checked(a.exposure)
    if r['status'] != 'complete_site_parsimony_exposure_diagnostic':
        raise ValueError('Incomplete diagnostic')
    for name in ['inputs', 'fits', 'exposure']:
        if r['source_receipts'][name] != sha(getattr(a, name)/'receipt.json'):
            raise ValueError('Changed source lineage')
    table = list(rows(a.diagnostic/'site_parsimony_exposure.tsv'))
    bykey = {(x['marker'], int(x['paired_column_1based'])): x for x in table}
    if len(bykey) != len(table) or len(table) != r['sites']:
        raise ValueError('Duplicate or incomplete site rows')
    observations = defaultdict(dict); count = 0
    with gzip.open(a.exposure/'normalized_paired_sites.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['marker'], int(row['paired_column_1based'])
            taxon = row['taxon_id']
            if key not in bykey or taxon in observations[key]:
                raise ValueError('Unexpected or duplicated residue observation')
            observations[key][taxon] = (row['amino_acid'], row['3di_state'], row['model_id'], row['rsa_tien2013_theoretical'], row['rsa_miller1987'])
            count += 1
    if count != er['rows'] or count != r['taxon_site_observations'] or set(observations) != set(bykey):
        raise ValueError('Exposure universe differs')
    ready = [x for x in rows(a.inputs/'marker_summary.tsv') if x['status']=='ready_for_inference']
    seen = set(); scores_checked = 0; quantiles_checked = 0
    for markerrow in ready:
        marker = markerrow['marker']; treepath = a.fits/marker/'aa.treefile'
        if sha(treepath) != r['tree_sha256'][marker]:
            raise ValueError('Changed reference tree')
        tree = Phylo.read(treepath, 'newick')
        aa = {x.id: str(x.seq) for x in SeqIO.parse(a.inputs/marker/'aa.faa','fasta')}
        di = {x.id: str(x.seq) for x in SeqIO.parse(a.inputs/marker/'3di.faa','fasta')}
        tips = [x.name for x in tree.get_terminals()]
        if len(tips)!=len(set(tips)) or set(tips)!=set(aa) or set(aa)!=set(di):
            raise ValueError('Taxon identity mismatch')
        columns = list(rows(a.inputs/marker/'columns.tsv'))
        independent = [set_scores(tree, seq) for seq in [aa,di]]
        if any(len(x)!=len(columns) for x in independent):
            raise ValueError('Column count mismatch')
        for j, col in enumerate(columns):
            key = marker, j+1; original = bykey[key]; obs = observations.pop(key)
            if key in seen or original['matrix_column_1based']!=col['matrix_column_1based']:
                raise ValueError('Site identity mismatch')
            seen.add(key); expected = {t for t,s in aa.items() if s[j]!='?'}
            if set(obs)!=expected or expected!={t for t,s in di.items() if s[j]!='?'}:
                raise ValueError('Observed taxa differ')
            if any(v[0]!=aa[t][j] or v[1]!=di[t][j] for t,v in obs.items()):
                raise ValueError('Observed states differ')
            if int(original['observed_taxa'])!=len(obs) or int(original['unique_models'])!=len({v[2] for v in obs.values()}):
                raise ValueError('Observation/model counts differ')
            for k,label in enumerate(['aa','3di']):
                if int(original[label+'_minimum_changes'])!=int(independent[k][j]) or int(original[label+'_distinct_states'])!=len({v[k] for v in obs.values()}):
                    raise ValueError('Independent character score differs')
                scores_checked += 1
            for label,k in [('tien2013',3),('miller1987',4)]:
                values=[float(v[k]) for v in obs.values() if v[k]!='']
                if int(original[label+'_normalized_taxa'])!=len(values):
                    raise ValueError('Normalized observation count differs')
                for stat,q in [('q25',.25),('median',.5),('q75',.75)]:
                    value=original[label+'_rsa_'+stat]
                    if (not values and value!='') or (values and not math.isclose(float(value),quantile(values,q),rel_tol=1e-12,abs_tol=1e-12)):
                        raise ValueError('Independent exposure quantile differs')
                    quantiles_checked += 1
        print(marker, 'full site readback passed', flush=True)
    if observations or seen!=set(bykey) or len(ready)!=r['markers']:
        raise ValueError('Incomplete readback')
    result={'status':'passed_full_site_parsimony_exposure_readback','markers':len(ready),'sites':len(seen),'taxon_site_observations':count,'independent_character_scores':scores_checked,'exposure_quantiles_checked':quantiles_checked,'source_receipt_sha256':sha(a.diagnostic/'receipt.json'),'script_sha256':sha(Path(__file__)),'score_helper_sha256':sha(Path(__file__).with_name('audit_site_parsimony_topologies.py')),'scope':'Every AA/3Di score independently recomputed by set-membership recurrence; all observed identities, states, model counts and exposure quantiles read back. No ancestral exposure, rate or coupling inference.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
