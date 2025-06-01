import argparse
import pandas as pd
import subprocess
from pathlib import Path

def batch_analyze(csv_path, bam_dir, output_dir, script_path="analyze_bam.py"):
    csv = pd.read_csv(csv_path)
    bam_dir = Path(bam_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for run_id in csv["Run"]:
        # Match files like SRR13404368*.bam
        bam_files = list(bam_dir.glob(f"{run_id}*.bam"))
        if not bam_files:
            print(f"[WARNING] BAM not found for {run_id}")
            continue

        bam_file = bam_files[0]
        prefix = output_dir / run_id
        cmd = [
            "python", script_path,
            str(bam_file),
            "--prefix", str(prefix),
            "--annotations", str(csv_path)
        ]
        print(f"[INFO] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="Metadata CSV with 'Run' column")
    parser.add_argument("bam_dir", help="Directory containing BAM files")
    parser.add_argument("output_dir", help="Directory for outputs")
    args = parser.parse_args()

    batch_analyze(args.csv, args.bam_dir, args.output_dir)
