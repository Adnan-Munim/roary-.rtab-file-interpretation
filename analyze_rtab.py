#!/usr/bin/env python3
"""
Analyze a Roary rarefaction .Rtab file
(number_of_new_genes.Rtab, number_of_conserved_genes.Rtab,
 number_of_genes_in_pan_genome.Rtab, number_of_unique_genes.Rtab)

Format: rows = permutation replicates (random genome orderings),
        columns = position when each genome was added (1st, 2nd, ...)

Usage:
    python3 analyze_rtab.py number_of_new_genes.Rtab
    python3 analyze_rtab.py number_of_new_genes.Rtab --out summary.csv --plot curve.png
"""
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser()
    p.add_argument('rtab_file', help='Path to a Roary .Rtab rarefaction file')
    p.add_argument('--out', default=None, help='Write per-position summary stats to this CSV')
    p.add_argument('--plot', default=None, help='Save a mean/median line plot to this file (e.g. curve.png)')
    args = p.parse_args()

    df = pd.read_csv(args.rtab_file, sep='\t', header=None)
    n_perms, n_genomes = df.shape
    print(f"File: {args.rtab_file}")
    print(f"Permutations (rows): {n_perms}")
    print(f"Genomes per permutation (columns): {n_genomes}")
    print()

    row_sums = df.sum(axis=1)
    print("Row sums (should be ~constant for number_of_new_genes.Rtab,")
    print("since order doesn't change the final pan-genome total):")
    print(f"  min={row_sums.min()}  max={row_sums.max()}  mean={row_sums.mean():.0f}")
    print()

    summary = pd.DataFrame({
        'genome_position': range(1, n_genomes + 1),
        'mean': df.mean(axis=0).round(1).values,
        'median': df.median(axis=0).values,
        'min': df.min(axis=0).values,
        'max': df.max(axis=0).values,
        'std': df.std(axis=0).round(1).values,
    })
    print(summary.to_string(index=False))

    if args.out:
        summary.to_csv(args.out, index=False)
        print(f"\nSummary written to {args.out}")

    if args.plot:
        plt.figure(figsize=(9, 5))
        plt.plot(summary['genome_position'], summary['mean'], label='Mean', marker='o', markersize=3)
        plt.plot(summary['genome_position'], summary['median'], label='Median', marker='o', markersize=3)
        plt.fill_between(summary['genome_position'], summary['min'], summary['max'],
                          alpha=0.15, label='Min-Max range')
        plt.xlabel('Genome added (position in random order)')
        plt.ylabel('Value')
        plt.title(args.rtab_file)
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.plot, dpi=300)
        print(f"Plot saved to {args.plot}")


if __name__ == '__main__':
    main()
