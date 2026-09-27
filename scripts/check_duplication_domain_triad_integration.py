#!/usr/bin/env python3
"""Reconstruct triad-domain intersections in Python and check every SQLite table/view row."""
import argparse,csv,hashlib,itertools,json,sqlite3
from collections import defaultdict,Counter
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text())
    assert r['status']=='complete_domain_triad_integration_pending_independent_readback' and r['plan_sha256']==sha(a.plan)
    assert r['biological_inference_eligible'] is False
    for p,h in plan['pins'].items():assert sha(p)==h,p
    assert sha(root/'domain_triads.sqlite')==r['artifacts']['domain_triads.sqlite']
    db=sqlite3.connect(f'file:{root}/domain_triads.sqlite?mode=ro',uri=True);db.row_factory=sqlite3.Row
    assert dict(db.execute('SELECT path,sha256 FROM provenance'))==plan['pins']
    tables=['triads','domain_links','coverage','sequence_covariates'];raw_counts={}
    # Compare full records directly to CSV, not just hashes or aggregate row counts.
    for table,source in zip(tables,plan['sources']):
        def original():
            for p in source['files']:
                with Path(p).open() as f:
                    for row in csv.DictReader(f,delimiter='\t'):
                        yield ({'guide':Path(p).name.split('_')[0],**row} if table=='sequence_covariates' else row)
        n=0
        for actual,expected in itertools.zip_longest(db.execute('SELECT * FROM '+table+' ORDER BY rowid'),original()):
            assert actual is not None and expected is not None and dict(actual)==expected,table
            n+=1
        assert n==r['counts'][table];raw_counts[table]=n
    keyfields=['guide','family','gene_node','gene_a','gene_b','reference_gene','policy']
    def key(row):return tuple(row[k] for k in keyfields)
    links=defaultdict(dict)
    for row in db.execute('SELECT * FROM domain_links'):
        d=dict(row);k=d['pair_key'],d['policy'],d['boundary'];pfam=d['pfam_accession']
        assert pfam not in links[k];links[k][pfam]=d
    expected={};availability={}
    for row in db.execute('SELECT * FROM triads'):
        t=dict(row);k=key(t)
        for boundary in ['alignment','envelope']:
            groups=[links.get((t[e+'_pair_key'],t['policy'],boundary),{}) for e in ['duplicate','a_reference','b_reference']]
            common=set(groups[0]) & set(groups[1]) & set(groups[2])
            status='matched_single_copy_domain' if common else ('identical_model_in_triad' if t['triad_class']=='identical_model_in_triad' else 'no_common_single_copy_domain_candidate')
            availability[k+(boundary,)]=dict(t,boundary=boundary,common_domain_count=len(common),domain_status=status)
            for pfam in common:
                domains=[g[pfam] for g in groups]
                maps=[{(d['model_'+side],d['version_'+side]):d['interval_'+side] for side in ['left','right']} for d in domains]
                ma=t['model_a'],t['version_a'];mb=t['model_b'],t['version_b'];mr=t['reference_model'],t['reference_version']
                assert set(maps[0])=={ma,mb} and set(maps[1])=={ma,mr} and set(maps[2])=={mb,mr}
                assert maps[0][ma]==maps[1][ma] and maps[0][mb]==maps[2][mb] and maps[1][mr]==maps[2][mr]
                x=dict(t,boundary=boundary,pfam_accession=pfam,duplicate_domain_pair_key=domains[0]['domain_pair_key'],a_reference_domain_pair_key=domains[1]['domain_pair_key'],b_reference_domain_pair_key=domains[2]['domain_pair_key'],interval_a=maps[0][ma],interval_b=maps[0][mb],interval_reference=maps[1][mr])
                expected[k+(boundary,pfam)]=x
    matched_seen=set()
    for row in db.execute('SELECT * FROM matched_domains'):
        x=dict(row);k=key(x)+(x['boundary'],x['pfam_accession']);assert k not in matched_seen and x==expected[k];matched_seen.add(k)
    assert matched_seen==set(expected)
    avail_seen=set();statuses=Counter()
    for row in db.execute('SELECT * FROM domain_availability'):
        x=dict(row);k=key(x)+(x['boundary'],);assert k not in avail_seen and x==availability[k];avail_seen.add(k)
        statuses[x['guide'],x['domain_status']]+=1
    assert avail_seen==set(availability)
    coverage={(x['pair_key'],x['mask']):dict(x) for x in db.execute('SELECT * FROM coverage')}
    seq={tuple(x[k] for k in keyfields[:5]):dict(x) for x in db.execute('SELECT * FROM sequence_covariates')}
    seqfields=['duplicate_pair_sequence_distance','sequence_tip_a','sequence_tip_b','sequence_signed_difference','sequence_normalized_contrast','sequence_direction','sequence_covariate_status']
    seen=set();totals=Counter();passes=Counter();refstatus=Counter();unique_events=set();screen_events=defaultdict(set)
    for row in db.execute('SELECT * FROM triad_domain_coverage'):
        actual=dict(row);mk=key(actual)+(actual['boundary'],actual['pfam_accession']);mask=actual['mask'];k=mk+(mask,)
        assert k not in seen;seen.add(k);t=expected[mk];s=seq[mk[:5]];wanted=dict(t,mask=mask)
        controls=[coverage[t[e+'_domain_pair_key'],mask] for e in ['duplicate','a_reference','b_reference']]
        for sid in plan['screens']:
            flag=int(all(c[sid+'_pass']=='1' for c in controls));wanted['all_three_'+sid+'_pass']=flag;passes[sid,t['guide'],mask]+=flag
            if flag:screen_events[sid,t['guide'],mask].add(mk[:5])
        wanted.update({f:s[f] for f in seqfields});exact=s['chosen_reference_gene']==t['reference_gene']
        wanted.update(sequence_covariate_reference=s['chosen_reference_gene'],sequence_reference_status='exact_reference' if exact else 'different_tied_reference',distance_a_to_this_reference=s['distance_a_to_reference'] if exact else None,distance_b_to_this_reference=s['distance_b_to_reference'] if exact else None)
        assert actual==wanted
        totals[t['guide'],mask]+=1;refstatus[wanted['sequence_reference_status']]+=1;unique_events.add(mk[:5])
    assert seen=={k+(mask,) for k in expected for mask in ['full','plddt70']}
    for item in r['screen_counts']:
        assert item['domain_triad_rows']==totals[item['guide'],item['mask']]
        assert item['all_three_pass']==passes[item['screen'],item['guide'],item['mask']]
    assert {(x['guide'],x['domain_status']):x['rows'] for x in r['availability_statuses']}==dict(statuses)
    assert len(expected)==r['counts']['matched_domains'] and len(availability)==r['counts']['domain_availability'] and len(seen)==r['counts']['triad_domain_coverage']
    assert db.execute('PRAGMA integrity_check').fetchone()[0]=='ok';db.close()
    assert sha(root/'domain_triads.sqlite')==r['artifacts']['domain_triads.sqlite']
    result=dict(status='passed_full_domain_triad_integration_readback',producer_receipt_sha256=sha(root/'receipt.json'),plan_sha256=sha(a.plan),checker_sha256=sha(__file__),raw_rows_checked=raw_counts,matched_domain_rows=len(expected),availability_rows=len(availability),domain_mask_rows=len(seen),guide_event_keys_with_matched_domains=len(unique_events),reference_status_counts=dict(refstatus),screen_event_counts=[dict(screen=sid,guide=g,mask=m,distinct_guide_events=len(v)) for (sid,g,m),v in sorted(screen_events.items())],biological_inference_eligible=False,scope='All imported CSV values checked. All three-way Pfam intersections and interval identities reconstructed independently in Python, exact availability and matched sets verified, every sequence field/reference identity/coverage decision checked. Guide events are not independent cross-guide observations; no significance or effect inference.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
