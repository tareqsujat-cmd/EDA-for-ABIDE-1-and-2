"""Non-interactive runner — downloads all EDA-relevant ABIDE data."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3
from botocore import UNSIGNED
from botocore.config import Config

BASE_DIR = Path(r"d:\EDA for ABIDE 1 and 2")
BUCKET = "fcp-indi"

s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))


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
    total = len(tasks)
    counts = {"ok": 0, "skip": 0, "fail": 0}
    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"  Files: {total}  |  Workers: {workers}")
    print(f"{'='*60}")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_one, k, p): (k, p) for k, p in tasks}
        for fut in as_completed(futures):
            result = fut.result()
            counts[result] += 1
            done = counts["ok"] + counts["skip"] + counts["fail"]
            print(
                f"  [{done:4d}/{total}]  ok={counts['ok']}  skip={counts['skip']}  fail={counts['fail']}",
                end="\r",
            )
    print(f"\n  Done — ok={counts['ok']}  skip={counts['skip']}  fail={counts['fail']}")
    return counts


# ── ABIDE 1 phenotypic ─────────────────────────────────────────────────────────
print("\n[1/3] ABIDE 1 — Phenotypic data")
raw_dir = BASE_DIR / "ABIDE_1" / "raw"
tasks = []
for key in [
    "data/Projects/ABIDE/Phenotypic_V1_0b.csv",
    "data/Projects/ABIDE/Phenotypic_V1_0b_preprocessed.csv",
    "data/Projects/ABIDE/Phenotypic_V1_0b_preprocessed1.csv",
]:
    tasks.append((key, raw_dir / Path(key).name))
site_dir = raw_dir / "phenotypic_sites"
for key, _ in list_s3_keys("data/Projects/ABIDE/PhenotypicData/", suffix_filter=[".csv"]):
    tasks.append((key, site_dir / Path(key).name))
download_batch(tasks, "ABIDE 1 — Phenotypic CSVs", workers=4)

# ── ABIDE 1 CC200 ROI time series ──────────────────────────────────────────────
print("\n[2/3] ABIDE 1 — CC200 ROI time series (nofilt_noglobal)")
out_dir = BASE_DIR / "ABIDE_1" / "preprocessed" / "rois_cc200"
prefix = "data/Projects/ABIDE/Outputs/cpac/nofilt_noglobal/rois_cc200/"
print("  Listing files…")
items = list_s3_keys(prefix, suffix_filter=[".1D"])
tasks = [(k, out_dir / Path(k).name) for k, _ in items]
total_mb = sum(s for _, s in items) / 1e6
print(f"  Found {len(tasks)} files  (~{total_mb:.0f} MB)")
download_batch(tasks, "ABIDE 1 — CC200 ROI time series", workers=8)

# ── ABIDE 2 DCAN Gordon2014 ────────────────────────────────────────────────────
print("\n[3/3] ABIDE 2 — DCAN Gordon2014 parcellated time series")
out_dir = BASE_DIR / "ABIDE_2" / "preprocessed" / "dcan_gordon2014"
prefix = "data/Projects/ABIDE2/Derivatives/DCAN/"
suffix = ["atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii"]
print("  Listing files…")
items = list_s3_keys(prefix, suffix_filter=suffix)
tasks = [(k, out_dir / Path(k).name) for k, _ in items]
total_mb = sum(s for _, s in items) / 1e6
print(f"  Found {len(tasks)} files  (~{total_mb:.0f} MB)")
download_batch(tasks, "ABIDE 2 — DCAN Gordon2014 parcellated time series", workers=8)

print("""
============================================================
All downloads complete.

NOTE — ABIDE 2 Phenotypic data is NOT on S3.
Download ABIDEII_Composite_Phenotypic.csv from:
  http://fcon_1000.projects.nitrc.org/indi/abide/abide_II.html
Place it in: ABIDE_2/raw/
============================================================
""")
