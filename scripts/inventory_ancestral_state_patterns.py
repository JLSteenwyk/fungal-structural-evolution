#!/usr/bin/env python3
"""Inventory lossless single-chain trace patterns; no quartet diagnostics."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from ancestral_chain_attempt import sha,write_json
from ancestral_state_patterns import group_patterns,reconstruct_patterns


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--completion',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    completion=json.loads(args.completion.read_text());rp=Path(completion['receipt']);assert sha(rp)==completion['receipt_sha256'];receipt=json.loads(rp.read_text())
    pins={str(p):sha(p) for p in [args.completion,rp,Path(__file__),Path('scripts/ancestral_state_patterns.py')]}
    rows=[];started=time.monotonic()
    for summary in receipt['summaries']:
        for p,h in summary['artifacts'].items():assert sha(p)==h
        path=next(Path(p) for p in summary['artifacts'] if p.endswith('/states.npz'))
        with np.load(path,allow_pickle=False) as saved:
            states=saved['states'];iterations=saved['iterations']
        for cutoff in [250,500]:
            x=states[iterations>cutoff][None,...]
            patterns,mapping=group_patterns(x)
            restored=reconstruct_patterns(patterns,mapping)
            assert np.array_equal(restored,x)
            assert np.array_equal(np.unique(mapping),np.arange(len(patterns)))
            rows.append(dict(chain=summary['chain'],discard_through=cutoff,draws=int(x.shape[1]),coordinates=int(mapping.size),
                unique_single_chain_patterns=len(patterns),raw_state_bytes=x.nbytes,grouped_state_and_map_bytes=patterns.nbytes+mapping.nbytes,
                independently_reconstructed_state_values=int(x.size),source=str(path),source_sha256=sha(path)))
    for p,h in pins.items():assert sha(p)==h
    write_json(args.output,dict(status='lossless_single_chain_pattern_inventory_complete',pins=pins,rows=rows,elapsed_seconds=time.monotonic()-started,
        scope='Only seven verified terminal chains, each analyzed separately. No synthetic quartet or production diagnostic. '
        'Full quartet unique-pattern counts may be larger; observed single-chain compression is not a full-grid runtime estimate.'))
    print(json.dumps(dict(rows=len(rows),coordinates=sum(r['coordinates'] for r in rows),patterns=sum(r['unique_single_chain_patterns'] for r in rows),elapsed_seconds=time.monotonic()-started)))


if __name__=='__main__':main()
