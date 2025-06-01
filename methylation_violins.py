import argparse
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set(style="whitegrid")


def expand_histograms(summary_dir):
    rows = []
    for path in Path(summary_dir).glob("*.summary.json"):
        with open(path) as f:
            data = json.load(f)
            sample = path.stem.replace(".summary", "")
            label = data.get("disease_status", "unknown")
            hist = data.get("cpg_methylation_hist", {})

            for bin_label, count in hist.items():
                try:
                    low, high = map(float, bin_label.split("-"))
                    center = (low + high) / 2
                    rows.extend([{
                        "sample_id": sample,
                        "disease_status": label,
                        "methylation_frac": center
                    }] * count)
                except Exception as e:
                    print(f"Skipping malformed bin: {bin_label} ({e})")
    return pd.DataFrame(rows)


def plot_violins(df, output_path):
    plt.figure(figsize=(max(10, len(df['sample_id'].unique()) * 0.6), 6))
    sns.violinplot(
        data=df,
        x="sample_id",
        y="methylation_frac",
        hue="disease_status",
        split=True,
        palette="Set2",
        inner="quartile"
    )
    plt.xticks(rotation=45, ha="right")
    plt.title("Per-Sample CpG Methylation Distribution")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summary_dir", help="Directory with .summary.json files")
    parser.add_argument("--output_dir", default=".", help="Directory to save output plot")
    args = parser.parse_args()

    df = expand_histograms(args.summary_dir)
    if df.empty:
        print("No methylation data found.")
        return

    output_path = Path(args.output_dir) / "methylation_violins.png"
    plot_violins(df, output_path)


if __name__ == "__main__":
    main()
