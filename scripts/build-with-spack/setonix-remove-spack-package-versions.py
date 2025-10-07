#!/usr/bin/env python3
"""Remove version constraints from Spack packages.yaml while preserving structure."""

import sys
import yaml

def remove_version_constraints(data):
    """Recursively remove 'version' keys from nested dictionaries."""
    if isinstance(data, dict):
        # Remove 'version' key if it exists
        if 'version' in data:
            del data['version']
        # Recursively process all values
        for key, value in data.items():
            remove_version_constraints(value)
    elif isinstance(data, list):
        # Recursively process list items
        for item in data:
            remove_version_constraints(item)
    return data

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.yaml> <output.yaml>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read input YAML
    with open(input_file, 'r') as f:
        data = yaml.safe_load(f)

    # Remove version constraints
    data = remove_version_constraints(data)

    # Write output YAML
    with open(output_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    print(f"Removed version constraints from {input_file} -> {output_file}")
