#!/usr/bin/env python3
#==================================================================================================
# stepOne -- isolate barcode reads from fastq
#
# Inputs:
#   1. Undetermined.fastq.gz generated after demultiplexing
#   2. Stagger file with format: sample_name<TAB>index<TAB>stagger_sequence
#         The stagger_sequence can be empty (no stagger) or a short nucleotide string.
#
# EIG updated Feb/Mar 2026:
#   fixed regex demultiplexing bug by using index to demux
#   fixed UMI bug to [16:28] from [17:26]
#   order of magnitude more efficiency using stream-process reads, regex.search, unique sets
#   added processing for Vpx barcode (Catch*ReSeq)
#
# Usage:
#   python stepOne.py \
#       --r1 /path/to/R1.fastq.gz \
#       --r2 /path/to/R2.fastq.gz \
#       --staggers /path/to/staggers.txt \
#       --outdir /path/to/output/ \
#       --mode GFP
#==================================================================================================

import argparse
import os
import sys
import regex
from Bio import SeqIO
from gzip import open as gzopen

VALID_BASES = set("ACGTN")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Demultiplex and filter paired FASTQ reads by index + stagger + primer matching."
    )
    parser.add_argument("--r1", required=True, help="Path to Read 1 FASTQ.gz (cell ID read)")
    parser.add_argument("--r2", required=True, help="Path to Read 2 FASTQ.gz (barcode read)")
    parser.add_argument("--staggers", required=True,
                        help="Tab-separated file: sample_name<TAB>index<TAB>stagger_sequence")
    parser.add_argument("--outdir", required=True, help="Root output directory")
    parser.add_argument("--mode", choices=["GFP", "Vpx"], default="GFP",
                        help="Primer mode: GFP (default) or Vpx")
    parser.add_argument("--index-mismatches", type=int, default=1,
                        help="Max allowed mismatches when matching index (default: 1)")
    return parser.parse_args()


def load_staggers(path):
    """Parse stagger file. Returns list of (sample_name, index, stagger_sequence) tuples.

    Each line is: sample_name<TAB>index<TAB>stagger_sequence
    The stagger_sequence can be empty (no stagger) or a short nucleotide string.
    Lines starting with # are comments. Blank lines are skipped.
    """
    samples = []
    seen_names = set()
    with open(path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                # sample_name + index, no stagger
                name = parts[0].strip()
                index_seq = parts[1].strip().upper()
                stagger_seq = ""
            elif len(parts) == 3:
                name = parts[0].strip()
                index_seq = parts[1].strip().upper()
                stagger_seq = parts[2].strip().upper()
            else:
                sys.exit(f"ERROR: staggers.txt line {line_num}: expected 2-3 tab-separated "
                         f"fields (sample_name, index, [stagger_sequence]), got {len(parts)}")

            if name in seen_names:
                sys.exit(f"ERROR: staggers.txt line {line_num}: duplicate sample name '{name}'")
            seen_names.add(name)

            if not index_seq:
                sys.exit(f"ERROR: staggers.txt line {line_num}: index cannot be empty for '{name}'")
            if not set(index_seq).issubset(VALID_BASES):
                bad = set(index_seq) - VALID_BASES
                sys.exit(f"ERROR: staggers.txt line {line_num}: index for '{name}' "
                         f"contains invalid characters: {bad}")
            if stagger_seq and not set(stagger_seq).issubset(VALID_BASES):
                bad = set(stagger_seq) - VALID_BASES
                sys.exit(f"ERROR: staggers.txt line {line_num}: stagger for '{name}' "
                         f"contains invalid characters: {bad}")

            samples.append((name, index_seq + "AT", stagger_seq))

    if not samples:
        sys.exit("ERROR: staggers.txt contains no samples")
    return samples


def extract_index(description):
    """Extract index sequence(s) from a FASTQ record description.
    Returns a tuple of individual index sequences, e.g. ('ATCACG',) or ('ATCACG', 'TTAGGC').
    """
    # New format: index is after the last colon in the second whitespace-delimited field
    parts = description.split()
    if len(parts) >= 2:
        # e.g. "1:N:0:ATCACG+TTAGGC" -> "ATCACG+TTAGGC"
        index_field = parts[1].rsplit(":", 1)[-1]
        return tuple(index_field.split("+"))

    # Old format: index between # and /
    if "#" in description:
        after_hash = description.split("#", 1)[1]
        index_field = after_hash.split("/")[0]
        return tuple(index_field.split("+"))

    return ()


def hamming_distance(s1, s2):
    """Count mismatches between two equal-length strings."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))  # length mismatch = treat as max distance
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def get_primer_seq(mode):
    """Return primer sequence based on mode."""
    if mode == "GFP":
        return "GGACGAGCTGTACAAGTAGG"
    elif mode == "Vpx":
        return "CTAGGGGAAGGACATGGGGCAGG"
    else:
        sys.exit(f"Unknown mode: {mode}")


class SampleState:
    """Per-sample output file handles, counters, and uniqueness sets."""

    def __init__(self, name, index_seq, stagger_seq, primer_seq, outdir, max_primer_errors=4):
        self.name = name
        self.index_seq = index_seq
        self.stagger_seq = stagger_seq
        self.stagger_len = len(stagger_seq)
        self.primer_seq = primer_seq
        self.primer_len = len(primer_seq)
        self.region_len = self.stagger_len + self.primer_len

        # Compile a per-sample regex: exact stagger match + fuzzy primer match.
        if stagger_seq:
            self.demux_pattern = regex.compile(
                rf'^({regex.escape(stagger_seq)})({primer_seq}){{e<={max_primer_errors}}}'
            )
        else:
            self.demux_pattern = regex.compile(
                rf'^({primer_seq}){{e<={max_primer_errors}}}'
            )

        self.barcode_slice = slice(self.region_len, self.region_len + 70)
        self.sample_dir = os.path.join(outdir, name)
        os.makedirs(self.sample_dir, exist_ok=True)

        # Output file handles
        self.f_joined = open(os.path.join(self.sample_dir, "joinedRead1Read2.txt"), "w")
        self.f_bad_qscore = open(os.path.join(self.sample_dir, "badQscore.txt"), "w")
        self.f_bad_barcode = open(os.path.join(self.sample_dir, "badBarcode.txt"), "w")
        self.f_unique_screened = open(os.path.join(self.sample_dir, "uniqueScreenedReads.txt"), "w")
        self.f_unique_shaved = open(os.path.join(self.sample_dir, "uniqueShavedReads.txt"), "w")

        # Counters
        self.total_by_index = 0      # reads assigned to this sample by index match
        self.total_by_primer = 0     # of those, reads with valid stagger+primer
        self.bad_qscore_count = 0
        self.bad_barcode_count = 0
        self.screened_count = 0
        self.shaved_count = 0
        self.missing_primer_count = 0

        # Uniqueness sets
        self.unique_screened = set()
        self.unique_shaved = set()

    def close(self):
        for f in (self.f_joined, self.f_bad_qscore, self.f_bad_barcode,
                  self.f_unique_screened, self.f_unique_shaved):
            f.close()

    def write_summary(self):
        with open(os.path.join(self.sample_dir, "summaryFile.txt"), "w") as summary:
            summary.write(f"sample {self.name}\n")
            summary.write(f"index {self.index_seq}\n")
            summary.write(f"stagger_sequence {self.stagger_seq if self.stagger_seq else '(none)'}\n")
            summary.write(f"stagger_length {self.stagger_len}\n")
            summary.write(f"primer_region_length {self.region_len}\n")
            summary.write(f"total reads by index {self.total_by_index}\n")
            summary.write(f"total assigned by stagger+primer {self.total_by_primer}\n")
            summary.write(f"total missingPrimer reads {self.missing_primer_count}\n")
            summary.write(f"total badQscore reads {self.bad_qscore_count}\n")
            summary.write(f"total badBarcode reads {self.bad_barcode_count}\n")
            summary.write(f"total screenedReads reads {self.screened_count}\n")
            summary.write(f"total screened Unique Reads reads {len(self.unique_screened)}\n")
            summary.write(f"total shavedReads reads {self.shaved_count}\n")
            summary.write(f"total shaved Unique Reads reads {len(self.unique_shaved)}\n")


def build_index_lookup(all_states):
    """Build a dict mapping exact index sequences to SampleState for O(1) lookup.

    Returns (lookup_dict, index_length). The index_length is used to validate
    that the provided index in the staggers file matches the length of the index
    in the FASTQ headers.
    """
    lookup = {}
    lengths = set()
    for state in all_states:
        idx = state.index_seq
        lengths.add(len(idx))
        if idx in lookup:
            sys.exit(f"ERROR: Duplicate index '{idx}' for samples "
                     f"'{lookup[idx].name}' and '{state.name}'")
        lookup[idx] = state
    return lookup, lengths


def match_index(read_indices, index_lookup, index_lengths, max_mismatches):
    """Match a read's index to a sample, allowing up to max_mismatches.

    The read's index field may be dual-indexed (i7+i5). Since the user provides only
    one index, we check each component of the read's index against all sample indices.

    Uses exact lookup first (O(1)), falls back to hamming distance scan only if needed.
    """
    for read_idx in read_indices:
        # Only compare if lengths match — different-length indices are from different
        # index reads (i7 vs i5) and should not be compared.
        if len(read_idx) not in index_lengths:
            continue

        # Fast path: exact match
        if read_idx in index_lookup:
            return index_lookup[read_idx]

        # Slow path: allow mismatches. Only reached for reads with sequencing
        # errors in the index. Iterates all sample indices but this is rare.
        if max_mismatches > 0:
            best_state = None
            best_dist = max_mismatches + 1
            for sample_idx, state in index_lookup.items():
                if len(sample_idx) != len(read_idx):
                    continue
                dist = hamming_distance(read_idx, sample_idx)
                if dist < best_dist:
                    best_dist = dist
                    best_state = state
                elif dist == best_dist and best_state is not None:
                    # Ambiguous: two indices equally close. Do not assign.
                    best_state = None
            if best_state is not None and best_dist <= max_mismatches:
                return best_state

    return None


def main():
    args = parse_args()

    # --- Load stagger definitions ---
    samples = load_staggers(args.staggers)
    primer_seq = get_primer_seq(args.mode)

    if args.mode == "Vpx":
        overridden = [(name, s) for name, _, s in samples if s]
        if overridden:
            for name, s in overridden:
                print(f"NOTE: stagger sequence for '{name}' ('{s}') ignored in Vpx mode.")
        samples = [(name, idx, "") for name, idx, _ in samples]

    # --- Build SampleState objects ---
    all_states = []
    for name, index_seq, stagger_seq in samples:
        state = SampleState(name, index_seq, stagger_seq, primer_seq, args.outdir)
        all_states.append(state)
    # Sort by stagger length descending for primer matching priority
    all_states.sort(key=lambda s: s.stagger_len, reverse=True)

    # --- Build index lookup ---
    index_lookup, index_lengths = build_index_lookup(all_states)

    # --- Auto-detect which index component (i7 vs i5) matches ---
    # Peek at the first read to show the user what indices are in the FASTQ
    detected_format = None
    with gzopen(args.r1, "rt") as peek_handle:
        first_record = next(SeqIO.parse(peek_handle, format="fastq"))
        read_indices = extract_index(first_record.description)
        if read_indices:
            detected_format = "+".join(read_indices)
            # Check which component matches any sample index
            matching_components = []
            for i, ri in enumerate(read_indices):
                if len(ri) in index_lengths:
                    matching_components.append(i)
            if not matching_components:
                print(f"WARNING: Index in FASTQ ({detected_format}) has no component matching "
                      f"the length of provided sample indices. Check your staggers file.")

    # --- Constants ---
    CELL_ID_SLICE = slice(0, 16)
    UMI_SLICE = slice(16, 28)
    MIN_QSCORE = 14
    MAX_LOW_QSCORE_BASES = 5
    bad_barcode_pattern = regex.compile(r'AAAA|TTTT|GGGG|CCCC|NN')

    # --- Global counters ---
    total_reads = 0
    no_index_reads = 0
    unmatched_index_reads = 0

    print(f"Mode: {args.mode} | Primer: {primer_seq} ({len(primer_seq)}bp)")
    print(f"Index mismatches allowed: {args.index_mismatches}")
    if detected_format:
        print(f"Index format in FASTQ (first read): {detected_format}")
    print(f"Samples: {len(all_states)}")
    for state in all_states:
        stag_label = state.stagger_seq if state.stagger_seq else "(none)"
        print(f"  {state.name}: index={state.index_seq}, stagger={stag_label} "
              f"({state.stagger_len}bp), region={state.region_len}bp")

    # --- Single-pass: read once, demultiplex by index then stagger+primer ---
    try:
        with gzopen(args.r1, "rt") as r1_handle, \
             gzopen(args.r2, "rt") as r2_handle:

            records_r1 = SeqIO.parse(r1_handle, format="fastq")
            records_r2 = SeqIO.parse(r2_handle, format="fastq")

            print("Processing reads (single pass, demultiplexing by index + stagger)...")

            for record1, record2 in zip(records_r1, records_r2):
                if record1.id != record2.id:
                    continue

                total_reads += 1

                # --- Step 1: Assign to sample by index ---
                read_indices = extract_index(record1.description)
                if not read_indices:
                    no_index_reads += 1
                    continue

                matched_state = match_index(
                    read_indices, index_lookup, index_lengths, args.index_mismatches
                )
                if matched_state is None:
                    unmatched_index_reads += 1
                    continue

                matched_state.total_by_index += 1

                seq1 = str(record1.seq)
                seq2 = str(record2.seq)
                r2_qscores = record2.letter_annotations["phred_quality"]

                # --- Step 2: Verify stagger+primer in read2 ---
                if not matched_state.demux_pattern.search(seq2):
                    matched_state.missing_primer_count += 1
                    continue

                matched_state.total_by_primer += 1
                matched_state.f_joined.write(f"{seq1} {seq2}\n")

                # --- Step 3: Quality and barcode filters ---

                # Filter 1: Quality score in stagger+primer region
                low_q_count = sum(
                    1 for q in r2_qscores[:matched_state.region_len] if q <= MIN_QSCORE
                )
                if low_q_count > MAX_LOW_QSCORE_BASES:
                    matched_state.bad_qscore_count += 1
                    matched_state.f_bad_qscore.write(f"{seq1} {seq2}\n")
                    continue

                # Filter 2: Homopolymer runs or ambiguous bases in read2 (GFP only)
                if args.mode != "Vpx" and bad_barcode_pattern.search(seq2):
                    matched_state.bad_barcode_count += 1
                    matched_state.f_bad_barcode.write(f"{seq1} {seq2}\n")
                    continue

                # Passed all filters
                matched_state.screened_count += 1
                screened_key = (seq1, seq2)
                if screened_key not in matched_state.unique_screened:
                    matched_state.unique_screened.add(screened_key)
                    matched_state.f_unique_screened.write(f"{seq1} {seq2}\n")

                cell_id = seq1[CELL_ID_SLICE]
                umi = seq1[UMI_SLICE]
                barcode = seq2[matched_state.barcode_slice]
                shaved_key = (cell_id, umi, barcode)
                matched_state.shaved_count += 1
                if shaved_key not in matched_state.unique_shaved:
                    matched_state.unique_shaved.add(shaved_key)
                    matched_state.f_unique_shaved.write(f"{cell_id} {umi} {barcode}\n")

                if total_reads % 1_000_000 == 0:
                    print(f"  Processed {total_reads:,} reads...")

    finally:
        for state in all_states:
            state.close()

    # --- Reporting ---
    print(f"Finished processing {total_reads:,} reads.")
    print(f"No index in header: {no_index_reads:,}")
    print(f"Unmatched index: {unmatched_index_reads:,}")

    for state in all_states:
        state.write_summary()
        print(f"  {state.name}: "
              f"{state.total_by_index:,} by index, "
              f"{state.total_by_primer:,} by stagger+primer, "
              f"{state.screened_count:,} screened, "
              f"{len(state.unique_shaved):,} unique shaved")

    # --- Write global summary ---
    with open(os.path.join(args.outdir, "globalSummary.txt"), "w") as gs:
        gs.write(f"mode {args.mode}\n")
        gs.write(f"primer {primer_seq}\n")
        gs.write(f"primer_length {len(primer_seq)}\n")
        gs.write(f"index_mismatches_allowed {args.index_mismatches}\n")
        gs.write(f"total raw reads {total_reads}\n")
        gs.write(f"no index in header {no_index_reads}\n")
        gs.write(f"unmatched index {unmatched_index_reads}\n")
        for state in all_states:
            stag_label = state.stagger_seq if state.stagger_seq else "(none)"
            gs.write(f"sample {state.name} "
                     f"index {state.index_seq} "
                     f"stagger {stag_label} "
                     f"reads_by_index {state.total_by_index} "
                     f"reads_by_primer {state.total_by_primer} "
                     f"missing_primer {state.missing_primer_count} "
                     f"screened {state.screened_count} "
                     f"unique_shaved {len(state.unique_shaved)}\n")

    print("Done. Output written to: " + args.outdir)


if __name__ == "__main__":
    main()