#!/usr/bin/env python3
"""Run the unchanged identity audit with a 16-MiB CSV field allowance.

Large native gene groups exceed Python's default 128-KiB field limit. This
entry point preserves the original producer used by the earlier run.
"""
import csv

from audit_grouped_ortholog_identities import main


if __name__ == '__main__':
    csv.field_size_limit(16 * 1024 * 1024)
    main()
