#!/usr/bin/env python3
"""Actual native write and normal-exit controls on two artificial words."""
import argparse
import ctypes
import os
from pathlib import Path
import signal

import psutil
from ancestral_chain_attempt import write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--case', choices=['native_write', 'normal_exit'], required=True)
    a = p.parse_args(); a.metadata.parent.mkdir(exist_ok=False)
    words = (ctypes.c_uint64*2)(0x1122334455667788, 0x8877665544332211)
    original = list(words); base = ctypes.addressof(words); inferior = psutil.Process()
    write_json(a.metadata, dict(targets={name: dict(address=base+i*8, bytes=8)
        for i, name in enumerate(['artificial_word_a', 'artificial_word_b'])},
        inferior_pid=inferior.pid, inferior_created=inferior.create_time(),
        process_local_read_only_protection=False, artificial_control=True))
    os.kill(os.getpid(), signal.SIGSTOP)
    if a.case == 'native_write':
        ctypes.memset(base, 0, 8)
        raise AssertionError('Hardware watchpoint must stop at the actual native byte write')
    assert list(words) == original
    write_json(a.receipt, dict(status='completed_artificial_watchpoint_normal_exit_target',
        words_preserved=True, artificial_control=True, scientific_eligibility=False))


if __name__ == '__main__': main()
