import pysam
import argparse
import pandas as pd
import json
import random
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path
import math

MAX_SAMPLE_ROWS = 10000

def extract_end_motifs(seq, length=4):
    return seq[:length], seq[-length:]

def compute_methylation_stats(df):
    methyl_fracs = df["cpg_methyl_frac"].dropna()
    if methyl_fracs.empty:
        return {
            "cpg_methylation_mean": None,
            "cpg_methylation_median": None,
            "cpg_methylation_std": None,
            "cpg_methylation_iqr": None,
            "cpg_methylation_entropy": None,
            "cpg_methylation_hist": {},
            "n_fragments_with_cpg": 0,
            "n_fragments_without_cpg": int(len(df))
        }

    hist_bins = np.linspace(0, 1, 11)
    hist_counts, _ = np.histogram(methyl_fracs, bins=hist_bins)
    bin_labels = [f"{round(hist_bins[i], 1)}-{round(hist_bins[i+1], 1)}" for i in range(len(hist_bins) - 1)]
    hist_dict = dict(zip(bin_labels, hist_counts.tolist()))

    probs = hist_counts / hist_counts.sum()
    entropy = -np.sum([p * math.log2(p) for p in probs if p > 0])

    return {
        "cpg_methylation_mean": float(methyl_fracs.mean()),
        "cpg_methylation_median": float(methyl_fracs.median()),
        "cpg_methylation_std": float(methyl_fracs.std()),
        "cpg_methylation_iqr": float(methyl_fracs.quantile(0.75) - methyl_fracs.quantile(0.25)),
        "cpg_methylation_entropy": float(entropy),
        "cpg_methylation_hist": hist_dict,
        "n_fragments_with_cpg": int(methyl_fracs.count()),
        "n_fragments_without_cpg": int(len(df) - methyl_fracs.count())
    }

def load_annotations(annotation_file):
    df = pd.read_csv(annotation_file)
    return df.set_index("Run").to_dict(orient="index")

def analyze_bam(bam_path, annotations=None):
    bam = pysam.AlignmentFile(bam_path, "rb")
    fragments = []
    xm_counter = defaultdict(int)

    for read in bam.fetch():
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if not read.is_proper_pair or not read.is_read1:
            continue

        start = read.reference_start
        end = read.next_reference_start + abs(read.template_length)
        insert_size = abs(read.template_length)

        if insert_size == 0 or insert_size > 1000:
            continue

        seq = read.query_sequence
        if seq is None or len(seq) < 4:
            continue

        motif_5p, motif_3p = extract_end_motifs(seq)

        xm = read.get_tag("XM") if read.has_tag("XM") else None
        z_count = xm.count("Z") if xm else 0
        z_unmeth_count = xm.count("z") if xm else 0
        cpg_total = z_count + z_unmeth_count
        cpg_methyl_frac = z_count / cpg_total if cpg_total > 0 else None

        if xm:
            for char in xm:
                xm_counter[char] += 1

        fragments.append({
            "start": start,
            "end": end,
            "length": insert_size,
            "motif_5p": motif_5p,
            "motif_3p": motif_3p,
            "xm": xm,
            "cpg_methylated": z_count,
            "cpg_unmethylated": z_unmeth_count,
            "cpg_total": cpg_total,
            "cpg_methyl_frac": cpg_methyl_frac
        })

    df = pd.DataFrame(fragments)

    if len(df) > MAX_SAMPLE_ROWS:
        df_sampled = df.sample(n=MAX_SAMPLE_ROWS, random_state=42)
    else:
        df_sampled = df

    summary = {
        "n_fragments": int(len(df)),
        "mean_length": float(df["length"].mean()),
        "median_length": float(df["length"].median()),
        "motif_5p_counts": df["motif_5p"].value_counts().to_dict(),
        "motif_3p_counts": df["motif_3p"].value_counts().to_dict(),
        "xm_counts": dict(xm_counter)
    }

    summary.update(compute_methylation_stats(df))

    # Add annotation if available
    sample_id = bam_path.name.split('.')[0]
    if annotations:
        if sample_id in annotations:
                    summary.update(annotations[sample_id])
        else:
            print(f"[WARNING] No annotation found for sample: {sample_id}")

    return df_sampled, summary

def main():
    parser = argparse.ArgumentParser(description="Analyze a bisulfite BAM file")
    parser.add_argument("bam", type=str, help="Input BAM file")
    parser.add_argument("--prefix", type=str, default=None, help="Output file prefix")
    parser.add_argument("--annotations", type=str, help="CSV file with Run metadata (Run column must match BAM filename)")
    args = parser.parse_args()

    bam_path = Path(args.bam)
    prefix = Path(args.prefix) if args.prefix else bam_path.with_suffix("")

    annotations = load_annotations(args.annotations) if args.annotations else None
    df_sampled, summary = analyze_bam(bam_path, annotations)

    df_sampled.to_csv(f"{prefix}.sampled.csv", index=False)
    with open(f"{prefix}.summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote {len(df_sampled)} sampled fragments to {prefix}.sampled.csv")
    print(f"Wrote summary statistics to {prefix}.summary.json")

if __name__ == "__main__":
    main()
