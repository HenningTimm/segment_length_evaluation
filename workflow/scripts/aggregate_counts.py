"""Aggregate several count files into one.

Usage:

    python aggregate_counts.py counts/synth/MRLs_synth0_non-canonical_no-mixing.dat [...]  counts/synth/MRLs_synth10_no-canonical_no-mixing.dat --output_path "foo.dat"
"""
import click
from collections import Counter


@click.command()
@click.argument("data_files", nargs=-1, type=click.File(mode="r"))
@click.option("--output_path", type=click.Path(exists=False, dir_okay=False),
              required=True)
def aggregate_files(data_files, output_path):
    """Collect counts from all input files and return them as a combined counter."""
    c = Counter()
    for data_file in data_files:
        for line in data_file:
            length, count = line.strip().split()
            c[int(length)] += int(count)
    with open(output_path, "w") as outfile:
        for mrl, count in sorted(c.items()):
            outfile.write(f"{mrl} {count}\n")


if __name__ == "__main__":
    aggregate_files()
