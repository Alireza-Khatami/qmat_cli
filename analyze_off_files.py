#!/usr/bin/env python3
"""
Analyze OFF files and generate CSV statistics.
Usage: python analyze_off_files.py <directory> [output.csv]
"""

import os
import sys
import csv
from pathlib import Path
from statistics import mean, median


def parse_off_file(filepath):
    """Parse an OFF file and return (nvertices, nedges, nfaces). Returns (None, None, None) on error."""
    try:
        with open(filepath, 'r') as f:
            header = f.readline().strip()
            if header != "OFF":
                print(f"Warning: {filepath} does not start with 'OFF'", file=sys.stderr)
                return None, None, None

            counts_line = f.readline().strip()
            if not counts_line:
                print(f"Warning: {filepath} has no counts line", file=sys.stderr)
                return None, None, None

            parts = counts_line.split()
            if len(parts) < 3:
                print(f"Warning: {filepath} counts line has fewer than 3 values", file=sys.stderr)
                return None, None, None

            nvertices = int(parts[0])
            nfaces = int(parts[1])
            nedges = int(parts[2])

            return nvertices, nedges, nfaces
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return None, None, None


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_off_files.py <directory> [output.csv]")
        sys.exit(1)

    directory = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "off_analysis.csv"

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a directory", file=sys.stderr)
        sys.exit(1)

    off_files = sorted(Path(directory).rglob("*_mat_simplified.off"))

    if not off_files:
        print(f"No *_mat_simplified.off files found in {directory}", file=sys.stderr)
        sys.exit(1)

    data = []
    vertices_list = []
    edges_list = []
    faces_list = []

    for filepath in off_files:
        nvertices, nedges, nfaces = parse_off_file(filepath)
        if nvertices is not None:
            data.append({
                'filename': filepath.name,
                'vertices': nvertices,
                'edges': nedges,
                'faces': nfaces
            })
            vertices_list.append(nvertices)
            edges_list.append(nedges)
            faces_list.append(nfaces)

    if not data:
        print("No valid OFF files were parsed", file=sys.stderr)
        sys.exit(1)

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['filename', '#vertices', '#edges', '#faces'])

        for row in data:
            writer.writerow([row['filename'], row['vertices'], row['edges'], row['faces']])

        writer.writerow(['MEAN',
                         f"{mean(vertices_list):.2f}",
                         f"{mean(edges_list):.2f}",
                         f"{mean(faces_list):.2f}"])
        writer.writerow(['MEDIAN',
                         f"{median(vertices_list):.2f}",
                         f"{median(edges_list):.2f}",
                         f"{median(faces_list):.2f}"])
        writer.writerow(['MIN', min(vertices_list), min(edges_list), min(faces_list)])
        writer.writerow(['MAX', max(vertices_list), max(edges_list), max(faces_list)])

    print(f"Analysis written to {output_csv} ({len(data)} files processed)")


if __name__ == '__main__':
    main()
