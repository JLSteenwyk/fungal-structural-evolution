#!/usr/bin/env python3
"""Run native reconciliation with an explicit task-stall allowance and tracing."""
import argparse
import functools
import json
import os
from pathlib import Path
import sys
import time
from orthofinder.gene_tree_inference import trees2ologs_of
from orthofinder.run.main import main as native_main


def install(stall_seconds,progress):
    if not 120 <= stall_seconds <= 86400:
        raise ValueError('Task-stall allowance must be between 120 seconds and one day')
    progress.mkdir(exist_ok=False)
    original_parallel=trees2ologs_of.RunOrthologsParallel
    original_analysis=trees2ologs_of.TreeAnalyser.AnalyseTree

    @functools.wraps(original_parallel)
    def parallel(*args,**kwargs):
        kwargs['STALL_TIMEOUT']=stall_seconds
        return original_parallel(*args,**kwargs)

    @functools.wraps(original_analysis)
    def analysis(self,iog):
        path=progress/f'worker-{os.getpid()}.jsonl'
        start=time.time();cpu=time.process_time()
        def write(status,**fields):
            with path.open('a') as f:
                f.write(json.dumps(dict(family=f'OG{iog:07d}',pid=os.getpid(),status=status,
                                       time_unix=time.time(),cpu_seconds=time.process_time()-cpu,**fields))+'\n')
        write('started')
        try:
            result=original_analysis(self,iog)
        except BaseException as exc:
            write('raised',error=repr(exc))
            raise
        write('analysis_returned',elapsed_seconds=time.time()-start,returned_none=result is None)
        return result

    trees2ologs_of.RunOrthologsParallel=parallel
    trees2ologs_of.TreeAnalyser.AnalyseTree=analysis


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stall-timeout-seconds',type=float,required=True)
    parser.add_argument('native_arguments',nargs=argparse.REMAINDER)
    args=parser.parse_args()
    native=args.native_arguments
    if native and native[0]=='--':native=native[1:]
    install(args.stall_timeout_seconds,Path.cwd()/'native_task_progress')
    print(json.dumps(dict(wrapper='task_progress_and_explicit_stall_allowance',stall_timeout_seconds=args.stall_timeout_seconds,
                          interpretation='Native inference code and data unchanged; per-family return does not prove valid biological inference or completed native output.')),flush=True)
    sys.argv=['orthofinder']+native
    native_main()
