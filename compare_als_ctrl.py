import argparse
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")


def load_summaries(json_dir):
    summaries = []
    for path in Path(json_dir).glob("*.summary.json"):
        with open(path) as f:
            data = json.load(f)
            data["sample_id"] = path.stem.replace(".summary", "")
            summaries.append(data)
    return pd.DataFrame(summaries)


def plot_fragment_lengths(df, output_dir):
    plt.figure(figsize=(8, 5))
    sns.boxplot(data=df, x="disease_status", y="mean_length")
    plt.title("Mean Fragment Length by Disease Status")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "fragment_length_by_status.png")
    plt.close()


def plot_methylation_distribution(df, output_dir):
    plt.figure(figsize=(8, 5))
    sns.boxplot(data=df, x="disease_status", y="cpg_methylation_mean")
    plt.title("Mean CpG Methylation per Sample")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "cpg_methylation_mean_by_status.png")
    plt.close()

    plt.figure(figsize=(8, 5))
    sns.boxplot(data=df, x="disease_status", y="cpg_methylation_entropy")
    plt.title("Methylation Entropy per Sample")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "cpg_methylation_entropy_by_status.png")
    plt.close()


def plot_most_common_motifs(df, motif_col, title, output_dir):
    all_counts = []
    for _, row in df.iterrows():
        counts = row[motif_col]
        if isinstance(counts, dict):
            for motif, count in counts.items():
                all_counts.append({
                    "motif": motif,
                    "count": count,
                    "disease_status": row["disease_status"]
                })
    motif_df = pd.DataFrame(all_counts)
    top_motifs = motif_df.groupby("motif")["count"].sum().nlargest(10).index
    motif_df = motif_df[motif_df["motif"].isin(top_motifs)]

    plt.figure(figsize=(10, 6))
    sns.barplot(data=motif_df, x="motif", y="count", hue="disease_status")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / f"{motif_col}_top_motifs.png")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Compare ALS vs Control samples")
    parser.add_argument("summary_dir", help="Directory containing .summary.json files")
    parser.add_argument("--output_dir", default=".", help="Directory to save output plots")
    args = parser.parse_args()

    df = load_summaries(args.summary_dir)

    plot_fragment_lengths(df, args.output_dir)
    plot_methylation_distribution(df, args.output_dir)
    plot_most_common_motifs(df, "motif_5p_counts", "Top 5' End Motifs by Disease Status", args.output_dir)
    plot_most_common_motifs(df, "motif_3p_counts", "Top 3' End Motifs by Disease Status", args.output_dir)


if __name__ == "__main__":
    main()
