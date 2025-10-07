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
        # Skip lines with "version:" on separate line
        if re.match(r'^\s+version:\s*$', line):
            skip_next = True
            continue

        # Skip inline version: [...]
        if re.search(r'\s+version:\s*\[', line):
            continue

        # Skip array elements after version: (lines starting with whitespace and -)
        if skip_next and re.match(r'^\s+-\s+', line):
            continue
        else:
            skip_next = False

        output_lines.append(line)

    # Fix empty package entries - add {} to lines that are just "package_name:"
    # Skip special entries like "all:" which have nested content
    for i, line in enumerate(output_lines):
        if re.match(r'^  \w[\w-]*:\s*$', line) and not re.match(r'^  all:\s*$', line):
            output_lines[i] = line.rstrip() + ' {}\n'

    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    print(f"Removed version constraints from {input_file} -> {output_file}")
