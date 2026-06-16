#!/usr/bin/env python3
"""
gc_content.py — Calculate GC content for each sequence in a FASTA file.

GC content is the percentage of bases in a DNA sequence that are
either Guanine (G) or Cytosine (C). It is a basic but useful metric
in genomics (e.g. GC-rich regions often correlate with gene density).

Usage:
    python gc_content.py sequences.fasta
"""

import sys


def parse_fasta(filepath):
    """
    Read a FASTA file and yield (header, sequence) pairs.

    A FASTA file looks like:
        >seq1 description
        GGCCATGC
        ATGC
        >seq2
        AATT

    Lines starting with '>' are headers; the lines below (until the
    next '>') are the sequence, which may be split across multiple lines.
    """
    header = None
    sequence_parts = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()          # remove trailing newline/spaces
            if not line:                 # skip blank lines
                continue
            if line.startswith(">"):     # this line is a header
                # If we were already building a sequence, yield it first
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:]        # drop the '>' character
                sequence_parts = []      # reset for the new sequence
            else:
                sequence_parts.append(line)  # part of the current sequence

    # Yield the final sequence after the loop ends
    if header is not None:
        yield header, "".join(sequence_parts)


def gc_content(sequence):
    """
    Return the GC percentage of a sequence.

    Counting is case-insensitive so 'g' and 'G' both count.
    Returns 0.0 for an empty sequence to avoid dividing by zero.
    """
    if len(sequence) == 0:
        return 0.0
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def main():
    # Expect exactly one argument: the path to the FASTA file
    if len(sys.argv) != 2:
        print("Usage: python gc_content.py <fasta_file>")
        sys.exit(1)

    filepath = sys.argv[1]

    found_any = False
    for header, sequence in parse_fasta(filepath):
        found_any = True
        gc = gc_content(sequence)
        # Print the sequence name and its GC% to 2 decimal places
        print(f"{header}\t{gc:.2f}%")

    # Handle an empty (or header-less) file gracefully
    if not found_any:
        print("No sequences found. Is the file empty or not in FASTA format?")


if __name__ == "__main__":
    main()
