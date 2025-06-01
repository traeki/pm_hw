import argparse
import json
import subprocess
import joblib
import pandas as pd
from pathlib import Path


def extract_features_from_summary(summary_path):
    with open(summary_path) as f:
        data = json.load(f)

    row = {}
    for feature in ["fragment_length_hist", "cpg_methylation_hist", "start_position_hist"]:
        hist = data.get(feature, {})
        total = sum(hist.values())
        for bin_label, count in hist.items():
            row[f"{feature}_{bin_label}"] = count / total if total > 0 else 0.0

    return pd.DataFrame([row])


def main():
    parser = argparse.ArgumentParser(description="Predict ALS vs CTRL for a new BAM file")
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("--model", default="als_ctrl_classifier.joblib", help="Trained model path")
    parser.add_argument("--analyze_script", default="analyze_bam.py", help="Path to analyze_bam.py script")
    args = parser.parse_args()

    prefix = Path(args.bam).with_suffix("")
    summary_file = Path(f"{prefix}.summary.json")

    # Run analyze_bam.py on the input BAM
    cmd = ["python", args.analyze_script, args.bam, "--prefix", str(prefix)]
    subprocess.run(cmd, check=True)

    if not summary_file.exists():
        print(f"Error: expected summary JSON not found: {summary_file}")
        return

    # Load model and predict
    clf = joblib.load(args.model)
    features = extract_features_from_summary(summary_file)
    probs = clf.predict_proba(features)[0]
    prediction = clf.predict(features)[0]
    pred_label = "ALS" if prediction == 1 else "CTRL"
    print(f"Prediction for {args.bam}: {pred_label}")
    print(f"Probabilities: CTRL={probs[0]:.3f}, ALS={probs[1]:.3f}")


if __name__ == "__main__":
    main()
