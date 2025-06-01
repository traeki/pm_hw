# ALS vs CTRL Classifier from cfDNA

This repository analyzes bisulfite-treated cfDNA from BAM files to classify ALS vs control samples using fragment size, start position, and CpG methylation features.

## 🔧 Setup

Create the environment:

```bash
conda env create -f environment.yml
conda activate cfdna
```

## 🧪 Training (optional)

To retrain the classifier from `.summary.json` outputs (and see associated statistics like precision/recall):

```bash
python als_ctrl_classifier.py path/to/summary_jsons --model_out als_ctrl_classifier.joblib
```

This trains a `RandomForestClassifier` and saves the model for later use.

## 🔍 Predict on New BAM

To classify a new BAM file:

```bash
python predict_als_ctrl.py sample.bam --model als_ctrl_classifier.joblib
```

This will:

1. Run `analyze_bam.py` on the BAM file
2. Extract features from the resulting `.summary.json`
3. Predict `ALS` or `CTRL`, and print class probabilities

## 📁 Notes

* `analyze_bam.py` outputs:

  * `sample.sampled.csv`: downsampled fragment details
  * `sample.summary.json`: fragment/methylation summaries used by the classifier
* The classifier uses only:

  * Normalized fragment length histogram
  * CpG methylation histogram
  * Start position histogram
