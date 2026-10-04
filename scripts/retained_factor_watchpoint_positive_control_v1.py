#!/usr/bin/env python3
"""Artificial native-byte write control for the project's hardware watchpoints."""
import argparse
import ctypes
import os
from pathlib import Path
import signal

import psutil
from ancestral_chain_attempt import write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata', type=Path, required=True)
    args = parser.parse_args(); args.metadata.parent.mkdir(exist_ok=False)
    words = (ctypes.c_uint64*2)(0x1122334455667788, 0x8877665544332211)
    base = ctypes.addressof(words); inferior = psutil.Process()
    write_json(args.metadata, dict(targets={name: dict(address=base+i*8, bytes=8)
        for i, name in enumerate(['artificial_word_a', 'artificial_word_b'])},
        inferior_pid=inferior.pid, inferior_created=inferior.create_time(),
        process_local_read_only_protection=False, artificial_control=True))
    os.kill(os.getpid(), signal.SIGSTOP)
    ctypes.memset(base, 0, 8)
    raise AssertionError('Watchpoint control should stop at the native byte write')


if __name__ == '__main__':
    main()
