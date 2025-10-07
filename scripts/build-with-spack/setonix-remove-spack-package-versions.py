#!/usr/bin/env python3
"""Remove version constraints from Spack packages.yaml while preserving structure."""

import sys
import re

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.yaml> <output.yaml>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, 'r') as f:
        lines = f.readlines()

    output_lines = []
    skip_next = False

    for i, line in enumerate(lines):
        # Skip lines with "version:" and the array values that follow
        if re.match(r'^\s+version:\s*$', line):
            skip_next = True
            continue

        # Skip array elements after version: (lines starting with whitespace and -)
        if skip_next and re.match(r'^\s+-\s+', line):
            continue
        else:
            skip_next = False

        output_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    print(f"Removed version constraints from {input_file} -> {output_file}")
