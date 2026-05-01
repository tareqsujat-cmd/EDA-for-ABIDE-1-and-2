"""
Download ABIDE 1 and ABIDE 2 datasets from FCP-INDI S3 (public, no credentials needed).

What this script downloads:
  ABIDE 1:
    - Phenotypic CSVs (consolidated + per-site) → ABIDE_1/raw/
    - CPAC preprocessed CC200 ROI time series (nofilt_noglobal) → ABIDE_1/preprocessed/rois_cc200/
      (~437 MB, 1102 subjects, standard atlas for autism FC research)

  ABIDE 2:
    - DCAN Gordon2014 parcellated resting-state time series → ABIDE_2/preprocessed/dcan_gordon2014/
      (~967 MB, 1009 sessions, 333+subcortical ROIs)
    - NOTE: ABIDE 2 phenotypic CSV is not hosted on S3. Download it manually from:
      http://fcon_1000.projects.nitrc.org/indi/abide/abide_II.html
      (file: ABIDEII_Composite_Phenotypic.csv) and place in ABIDE_2/raw/

Requires: pip install boto3 botocore
"""

import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3
from botocore import UNSIGNED
from botocore.config import Config

BASE_DIR = Path(__file__).parent
BUCKET = "fcp-indi"

s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))


# ── helpers ────────────────────────────────────────────────────────────────────

def list_s3_keys(prefix, suffix_filter=None):
    keys = []
    token = None
    while True:
        kwargs = {"Bucket": BUCKET, "Prefix": prefix, "MaxKeys": 1000}
        if token:
            kwargs["ContinuationToken"] = token
        resp = s3.list_objects_v2(**kwargs)
        for obj in resp.get("Contents", []):
            key = obj["Key"]
            if key.endswith("/"):
                continue
            if suffix_filter and not any(key.endswith(s) for s in suffix_filter):
                continue
            keys.append((key, obj["Size"]))
        if not resp.get("IsTruncated"):
            break
        token = resp.get("NextContinuationToken")
    return keys


def download_one(s3_key, local_path):
    local_path.parent.mkdir(parents=True, exist_ok=True)
    if local_path.exists() and local_path.stat().st_size > 0:
        return "skip"
    try:
        s3.download_file(BUCKET, s3_key, str(local_path))
        return "ok"
    except Exception as e:
        print(f"\n  ERROR {s3_key}: {e}")
        return "fail"


def download_batch(tasks, label, workers=8):
    """tasks: list of (s3_key, local_path)"""
    total = len(tasks)
    counts = {"ok": 0, "skip": 0, "fail": 0}
    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"  Files: {total}  |  Workers: {workers}")
    print(f"{'='*60}")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_one, k, p): (k, p) for k, p in tasks}
        for i, fut in enumerate(as_completed(futures), 1):
            result = fut.result()
            counts[result] += 1
            done = counts["ok"] + counts["skip"] + counts["fail"]
            print(
                f"  [{done:4d}/{total}]  ok={counts['ok']}  skip={counts['skip']}  fail={counts['fail']}",
                end="\r",
            )
    print(f"\n  Done — ok={counts['ok']}  skip={counts['skip']}  fail={counts['fail']}")
    return counts


# ── ABIDE 1 ────────────────────────────────────────────────────────────────────

def download_abide1_phenotypic():
    raw_dir = BASE_DIR / "ABIDE_1" / "raw"
    tasks = []

    # consolidated phenotypic CSVs at bucket root
    for key in [
        "data/Projects/ABIDE/Phenotypic_V1_0b.csv",
        "data/Projects/ABIDE/Phenotypic_V1_0b_preprocessed.csv",
        "data/Projects/ABIDE/Phenotypic_V1_0b_preprocessed1.csv",
    ]:
        tasks.append((key, raw_dir / Path(key).name))

    # per-site phenotypic CSVs
    site_dir = raw_dir / "phenotypic_sites"
    for key, _ in list_s3_keys("data/Projects/ABIDE/PhenotypicData/", suffix_filter=[".csv"]):
        tasks.append((key, site_dir / Path(key).name))

    download_batch(tasks, "ABIDE 1 — Phenotypic CSVs", workers=4)


def download_abide1_rois_cc200():
    out_dir = BASE_DIR / "ABIDE_1" / "preprocessed" / "rois_cc200"
    prefix = "data/Projects/ABIDE/Outputs/cpac/nofilt_noglobal/rois_cc200/"
    print(f"\nListing ABIDE 1 CC200 files…")
    items = list_s3_keys(prefix, suffix_filter=[".1D"])
    tasks = [(k, out_dir / Path(k).name) for k, _ in items]
    total_mb = sum(s for _, s in items) / 1e6
    print(f"  Found {len(tasks)} files  (~{total_mb:.0f} MB)")
    download_batch(tasks, "ABIDE 1 — CC200 ROI time series (nofilt_noglobal, CPAC)", workers=8)


# ── ABIDE 2 ────────────────────────────────────────────────────────────────────

def download_abide2_dcan_gordon():
    out_dir = BASE_DIR / "ABIDE_2" / "preprocessed" / "dcan_gordon2014"
    prefix = "data/Projects/ABIDE2/Derivatives/DCAN/"
    suffix = ["atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii"]
    print(f"\nListing ABIDE 2 DCAN Gordon2014 files…")
    items = list_s3_keys(prefix, suffix_filter=suffix)
    tasks = [(k, out_dir / Path(k).name) for k, _ in items]
    total_mb = sum(s for _, s in items) / 1e6
    print(f"  Found {len(tasks)} files  (~{total_mb:.0f} MB)")
    download_batch(tasks, "ABIDE 2 — DCAN Gordon2014 parcellated time series", workers=8)


def note_abide2_phenotypic():
    print("""
NOTE — ABIDE 2 Phenotypic data
  The composite phenotypic CSV for ABIDE 2 is not hosted on the FCP-INDI S3 bucket.
  Download ABIDEII_Composite_Phenotypic.csv manually from:
    http://fcon_1000.projects.nitrc.org/indi/abide/abide_II.html
  Then place it at:
    ABIDE_2/raw/ABIDEII_Composite_Phenotypic.csv
""")


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    print("ABIDE EDA Dataset Downloader")
    print("Downloads phenotypic data + preprocessed ROI time series for EDA.\n")

    print("What will be downloaded:")
    print("  [1] ABIDE 1 — Phenotypic CSVs             (~1.2 MB)")
    print("  [2] ABIDE 1 — CC200 ROI time series       (~437 MB, 1102 subjects)")
    print("  [3] ABIDE 2 — DCAN Gordon2014 time series (~967 MB, 1009 sessions)")
    print("  [4] All of the above (recommended for EDA)")
    print("  [0] Exit")

    choice = input("\nEnter choice: ").strip()

    if choice == "0":
        sys.exit(0)
    elif choice == "1":
        download_abide1_phenotypic()
        note_abide2_phenotypic()
    elif choice == "2":
        download_abide1_rois_cc200()
    elif choice == "3":
        download_abide2_dcan_gordon()
        note_abide2_phenotypic()
    elif choice == "4":
        download_abide1_phenotypic()
        download_abide1_rois_cc200()
        download_abide2_dcan_gordon()
        note_abide2_phenotypic()
    else:
        print("Invalid choice.")
        sys.exit(1)

    print("\nAll downloads complete.")


if __name__ == "__main__":
    main()
