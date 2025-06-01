import argparse
import json
from pathlib import Path
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import make_pipeline
import joblib
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns


def load_features(json_dir):
    rows = []
    for path in Path(json_dir).glob("*.summary.json"):
        with open(path) as f:
            data = json.load(f)
        if data.get("disease_status") not in {"als", "ctrl"}:
            continue
        row = {"sample_id": path.stem.replace(".summary", ""), "disease_status": data["disease_status"]}

        # Normalize histograms and flatten into features
        for feature in ["fragment_length_hist", "cpg_methylation_hist", "start_position_hist"]:
            hist = data.get(feature, {})
            total = sum(hist.values())
            for bin_label, count in hist.items():
                row[f"{feature}_{bin_label}"] = count / total if total > 0 else 0.0

        rows.append(row)
    return pd.DataFrame(rows)


def train_and_evaluate(df, model_output_path):
    df = df.copy()
    df["label"] = df["disease_status"].map({"als": 1, "ctrl": 0})
    X = df.drop(columns=["sample_id", "disease_status", "label"])
    y = df["label"]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    clf = RandomForestClassifier(random_state=42)
    f1_scores = cross_val_score(clf, X_scaled, y, cv=5, scoring="f1")
    precision_scores = cross_val_score(clf, X_scaled, y, cv=5, scoring="precision")
    recall_scores = cross_val_score(clf, X_scaled, y, cv=5, scoring="recall")
    print(f"Cross-validated F1 scores: {f1_scores}")
    print(f"Cross-validated Precision scores: {precision_scores}")
    print(f"Cross-validated Recall (Sensitivity) scores: {recall_scores}")
    print(f"Mean F1: {f1_scores.mean():.3f} ± {f1_scores.std():.3f}")
    print(f"Mean Precision: {precision_scores.mean():.3f} ± {precision_scores.std():.3f}")
    print(f"Mean Recall: {recall_scores.mean():.3f} ± {recall_scores.std():.3f}")

    # Fit model on full data and save it
    clf_pipeline = make_pipeline(StandardScaler(), RandomForestClassifier(random_state=42))
    clf_pipeline.fit(X, y)
    joblib.dump(clf_pipeline, model_output_path)
    print(f"Saved trained classifier to {model_output_path}")


def main():
    parser = argparse.ArgumentParser(description="Train ALS/CTRL classifier and save model")
    parser.add_argument("summary_dir", help="Directory of .summary.json files")
    parser.add_argument("--model_out", default="als_ctrl_classifier.joblib", help="Path to save trained model")
    args = parser.parse_args()

    df = load_features(args.summary_dir)
    if df.empty:
        print("No valid summary files found.")
        return

    train_and_evaluate(df, args.model_out)


if __name__ == "__main__":
    main()
