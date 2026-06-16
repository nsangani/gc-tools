"""
gc_content.py
-------------
Reads one or more FASTA files and prints the GC% for every sequence.

Usage:
    python gc_content.py file1.fasta [file2.fasta ...]

Output format:
    <sequence_name>  GC: 52.34%
"""

import sys
import os


def parse_fasta(filepath):
    """
    Parse a FASTA file and yield (name, sequence) tuples.

    Each record starts with a header line beginning with '>'.
    Sequence lines that follow (until the next '>' or EOF) are
    concatenated into a single uppercase string.

    Raises:
        FileNotFoundError  – if the path does not exist.
        ValueError         – if the file is empty or has no valid records.
    """
    # Verify the file exists before opening
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    records = []  # Will hold (name, sequence) pairs

    current_name = None
    current_seq_parts = []  # Accumulate sequence lines here

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()

            # Skip blank lines
            if not line:
                continue

            if line.startswith(">"):
                # Save the previous record before starting a new one
                if current_name is not None:
                    records.append((current_name, "".join(current_seq_parts)))

                # The header may contain description text after a space;
                # use only the first word (the sequence ID) as the name.
                current_name = line[1:].split()[0] if len(line) > 1 else "unnamed"
                current_seq_parts = []
            else:
                # Accumulate sequence characters, normalise to uppercase
                current_seq_parts.append(line.upper())

    # Don't forget the last record in the file
    if current_name is not None:
        records.append((current_name, "".join(current_seq_parts)))

    # Raise an error for completely empty or header-only files
    if not records:
        raise ValueError(f"No valid FASTA records found in: {filepath}")

    return records


def calc_gc(sequence):
    """
    Calculate GC content as a percentage.

    Only G and C nucleotides count toward the numerator.
    Non-standard characters (gaps, ambiguity codes, etc.) are
    included in the total length, so they dilute the GC%.

    Returns:
        float  – GC percentage in [0.0, 100.0], or 0.0 for empty sequences.
    """
    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100.0


def process_file(filepath):
    """
    Parse a FASTA file and print GC% for each sequence.

    Handles errors gracefully so that a bad file does not abort
    processing of subsequent files on the command line.
    """
    print(f"\n=== {filepath} ===")

    try:
        records = parse_fasta(filepath)
    except (FileNotFoundError, ValueError) as exc:
        # Report the problem but keep going with other files
        print(f"  [ERROR] {exc}")
        return

    for name, sequence in records:
        gc = calc_gc(sequence)

        # Warn about sequences that contain no recognised bases
        if not sequence:
            print(f"  {name}  [WARNING: empty sequence]")
        else:
            print(f"  {name}  GC: {gc:.2f}%")


def main():
    """Entry point: expects at least one FASTA file path as a CLI argument."""
    if len(sys.argv) < 2:
        print("Usage: python gc_content.py file1.fasta [file2.fasta ...]")
        sys.exit(1)

    # Process every file provided on the command line
    for filepath in sys.argv[1:]:
        process_file(filepath)


if __name__ == "__main__":
    main()
