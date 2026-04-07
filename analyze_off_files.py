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


def parse_ma_file(filepath):
    """Parse a .ma file and return (nvertices, nedges, nfaces) from the first line. Returns (None, None, None) on error."""
    try:
        with open(filepath, 'r') as f:
            line = f.readline().strip()
            if not line:
                print(f"Warning: {filepath} is empty", file=sys.stderr)
                return None, None, None

            parts = line.split()
            if len(parts) < 3:
                print(f"Warning: {filepath} first line has fewer than 3 values", file=sys.stderr)
                return None, None, None

            nvertices = int(parts[0])
            nfaces = int(parts[1])
            nedges = int(parts[2])

            return nvertices, nedges, nfaces
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return None, None, None


def find_init_ma(off_path):
    """Given a *_mat_simplified.off path, find the corresponding init_vor_mat .ma file."""
    basename = off_path.name[: -len("_mat_simplified.off")]
    ma_path = off_path.parent / "init_vor_mat" / (basename + ".ma")
    return ma_path if ma_path.exists() else None


def ratio_str(simplified, initial):
    """Return simplified/initial as a formatted string, or 'N/A' if initial is missing or zero."""
    if initial is None or initial == 0:
        return "N/A"
    return f"{simplified / initial:.4f}"


def main():
    # if len(sys.argv) < 2:
    #     print("Usage: python analyze_off_files.py <directory> [output.csv]")
    #     sys.exit(1)

    directory = sys.argv[1] if len(sys.argv) > 1 else "output/abc_dataset"
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "off_analysis.csv"

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a directory", file=sys.stderr)
        sys.exit(1)

    off_files = sorted(Path(directory).rglob("*_mat_simplified.off"))

    if not off_files:
        print(f"No *_mat_simplified.off files found in {directory}", file=sys.stderr)
        sys.exit(1)

    data = []
    simp_vertices_list = []
    simp_edges_list = []
    simp_faces_list = []
    init_vertices_list = []
    init_edges_list = []
    init_faces_list = []

    for filepath in off_files:
        nv, ne, nf = parse_off_file(filepath)
        if nv is None:
            continue

        ma_path = find_init_ma(filepath)
        iv, ie, if_ = parse_ma_file(ma_path) if ma_path else (None, None, None)
        if ma_path is None:
            print(f"Info: no init_vor_mat .ma found for {filepath.name}", file=sys.stderr)

        data.append({
            'filename': filepath.name,
            'simp_v': nv, 'simp_e': ne, 'simp_f': nf,
            'init_v': iv, 'init_e': ie, 'init_f': if_,
        })
        simp_vertices_list.append(nv)
        simp_edges_list.append(ne)
        simp_faces_list.append(nf)
        if iv is not None:
            init_vertices_list.append(iv)
            init_edges_list.append(ie)
            init_faces_list.append(if_)

    if not data:
        print("No valid OFF files were parsed", file=sys.stderr)
        sys.exit(1)

    def fmt(v):
        return v if v is not None else "N/A"

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'filename',
            'simp_#vertices', 'simp_#edges', 'simp_#faces',
            'init_#vertices', 'init_#edges', 'init_#faces',
            'ratio_vertices', 'ratio_edges', 'ratio_faces',
        ])

        for row in data:
            writer.writerow([
                row['filename'],
                row['simp_v'], row['simp_e'], row['simp_f'],
                fmt(row['init_v']), fmt(row['init_e']), fmt(row['init_f']),
                ratio_str(row['simp_v'], row['init_v']),
                ratio_str(row['simp_e'], row['init_e']),
                ratio_str(row['simp_f'], row['init_f']),
            ])

        # Summary rows (simplified only — always available)
        writer.writerow([])
        writer.writerow(['--- simplified stats ---'])
        writer.writerow(['MEAN',
                         f"{mean(simp_vertices_list):.2f}", f"{mean(simp_edges_list):.2f}", f"{mean(simp_faces_list):.2f}",
                         f"{mean(init_vertices_list):.2f}" if init_vertices_list else 'N/A',
                         f"{mean(init_edges_list):.2f}" if init_edges_list else 'N/A',
                         f"{mean(init_faces_list):.2f}" if init_faces_list else 'N/A',
                         '', '', ''])
        writer.writerow(['MEDIAN',
                         f"{median(simp_vertices_list):.2f}", f"{median(simp_edges_list):.2f}", f"{median(simp_faces_list):.2f}",
                         f"{median(init_vertices_list):.2f}" if init_vertices_list else 'N/A',
                         f"{median(init_edges_list):.2f}" if init_edges_list else 'N/A',
                         f"{median(init_faces_list):.2f}" if init_faces_list else 'N/A',
                         '', '', ''])
        writer.writerow(['MIN',
                         min(simp_vertices_list), min(simp_edges_list), min(simp_faces_list),
                         min(init_vertices_list) if init_vertices_list else 'N/A',
                         min(init_edges_list) if init_edges_list else 'N/A',
                         min(init_faces_list) if init_faces_list else 'N/A',
                         '', '', ''])
        writer.writerow(['MAX',
                         max(simp_vertices_list), max(simp_edges_list), max(simp_faces_list),
                         max(init_vertices_list) if init_vertices_list else 'N/A',
                         max(init_edges_list) if init_edges_list else 'N/A',
                         max(init_faces_list) if init_faces_list else 'N/A',
                         '', '', ''])

    init_count = len(init_vertices_list)
    print(f"Analysis written to {output_csv} ({len(data)} files processed, {init_count} with init_vor_mat data)")


if __name__ == '__main__':
    main()
