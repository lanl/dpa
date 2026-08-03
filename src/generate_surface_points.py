#!/usr/bin/env python3
"""
Generate surface points around a protein for DPA analysis

Creates a shell of points at a specified distance from the protein surface.
Simple approach: points that are 2-4 Å from the nearest atom.

Usage:
    python generate_surface_points.py input.pdb output.hetcrd
"""

import numpy as np
from Bio.PDB import PDBParser
import argparse
import sys
from scipy.spatial import cKDTree

def load_structure(pdb_file, chain_id=None):
    """Load protein structure and extract heavy atom coordinates"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    coords = []
    for atom in structure.get_atoms():
        if atom.element != 'H':  # Skip hydrogens
            if chain_id is None or atom.get_parent().get_parent().id == chain_id:
                coords.append(atom.coord)

    coords = np.array(coords)
    print(f"Loaded {len(coords)} heavy atoms from {pdb_file}")
    return coords

def generate_grid(coords, boundary=6.0, spacing=1.0):
    """Generate 3D grid around protein"""
    mins = coords.min(axis=0) - boundary
    maxs = coords.max(axis=0) + boundary

    x = np.arange(mins[0], maxs[0], spacing)
    y = np.arange(mins[1], maxs[1], spacing)
    z = np.arange(mins[2], maxs[2], spacing)

    xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')
    grid_points = np.column_stack([xv.ravel(), yv.ravel(), zv.ravel()])

    print(f"Generated {len(grid_points):,} grid points (spacing={spacing} Å)")
    return grid_points

def find_surface_shell(grid_points, atom_coords, min_dist=2.0, max_dist=4.0):
    """Find points in shell around protein surface"""
    print(f"Finding surface shell ({min_dist}-{max_dist} Å)...")

    tree = cKDTree(atom_coords)

    surface_points = []
    batch_size = 100000

    for i in range(0, len(grid_points), batch_size):
        batch = grid_points[i:i+batch_size]
        distances, _ = tree.query(batch)
        mask = (distances >= min_dist) & (distances <= max_dist)
        surface_points.append(batch[mask])

        if (i // batch_size + 1) % 5 == 0:
            print(f"  {100*i//len(grid_points)}%...")

    return np.vstack(surface_points)

def subsample_uniform(points, max_points):
    """Uniform spatial subsampling (optional)"""
    if max_points is None or len(points) <= max_points:
        return points

    print(f"Subsampling to {max_points} points...")

    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    volume = np.prod(maxs - mins)
    voxel_size = (volume / max_points) ** (1/3)

    voxel_indices = ((points - mins) / voxel_size).astype(int)
    unique_voxels = {}

    for i, voxel in enumerate(voxel_indices):
        key = tuple(voxel)
        if key not in unique_voxels:
            unique_voxels[key] = i

    selected = list(unique_voxels.values())[:max_points]
    return points[selected]

def write_hetcrd(points, output_file):
    """Write .hetcrd format for DPA"""
    with open(output_file, 'w') as f:
        for i, point in enumerate(points, start=1):
            f.write(f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  {i:4d} SURF\n")
    print(f"Wrote {len(points)} points to {output_file}")

def write_pdb(points, output_file):
    """Write PDB for visualization"""
    with open(output_file, 'w') as f:
        for i, point in enumerate(points, start=1):
            f.write(f"HETATM{i:5d}  CA  SRF A{i:4d}    "
                   f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}"
                   f"  1.00 50.00           C\n")
    print(f"Wrote visualization to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate surface points for DPA')
    parser.add_argument('input_pdb', help='Input PDB file')
    parser.add_argument('output_hetcrd', help='Output .hetcrd file')
    parser.add_argument('--chain', help='Specific chain ID')
    parser.add_argument('--spacing', type=float, default=0.75, help='Grid spacing (Å, default: 0.75)')
    parser.add_argument('--min-distance', type=float, default=3.45, help='Min distance from atoms (Å, default: 3.45)')
    parser.add_argument('--max-distance', type=float, default=3.55, help='Max distance from atoms (Å, default: 3.55)')
    parser.add_argument('--max-points', type=int, default=None, help='Max output points (default: no limit)')
    parser.add_argument('--output-pdb', help='Output PDB for visualization')

    args = parser.parse_args()

    print("="*60)
    print("Surface Point Generator")
    print("="*60)

    # Load structure
    coords = load_structure(args.input_pdb, args.chain)

    # Generate grid
    grid = generate_grid(coords, boundary=6.0, spacing=args.spacing)

    # Find surface shell
    surface = find_surface_shell(grid, coords, args.min_distance, args.max_distance)
    print(f"Found {len(surface):,} surface points")

    # Subsample if requested
    if args.max_points is not None and len(surface) > args.max_points:
        surface = subsample_uniform(surface, args.max_points)

    # Write outputs
    write_hetcrd(surface, args.output_hetcrd)
    if args.output_pdb:
        write_pdb(surface, args.output_pdb)

    print("="*60)
    print(f"Done! {len(surface)} points written")
    if args.max_points:
        print(f"(max limit: {args.max_points})")
    else:
        print("(no subsampling)")
    print("="*60)

if __name__ == '__main__':
    main()
