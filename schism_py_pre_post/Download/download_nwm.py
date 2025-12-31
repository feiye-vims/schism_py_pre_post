#!/usr/bin/env python3
import os
import argparse
from datetime import datetime, timedelta
import boto3
from botocore import UNSIGNED
from botocore.client import Config
from botocore.exceptions import ClientError

BUCKET = "noaa-nwm-retrospective-3-0-pds"
ALL_TYPES = ["CHANOBS", "CHRTOUT", "FORCING", "GWOUT", "LAKEOUT", "LDASOUT", "RTOUT"]

STEM_OVERRIDES = {
    "FORCING": "LDASIN",
    # Add more if you find exceptions later; default uses the folder/type itself as stem.
}


def filename_stem(ftype: str) -> str:
    return STEM_OVERRIDES.get(ftype, ftype)


def key_for(region, fmt, ftype, ts):
    # e.g., CONUS/netcdf/CHRTOUT/2017/201701010000.CHRTOUT_DOMAIN1
    stem = filename_stem(ftype)
    return f"{region}/{fmt}/{ftype}/{ts[:4]}/{ts}.{stem}_DOMAIN1"


def exists(s3, key):
    try:
        s3.head_object(Bucket=BUCKET, Key=key)
        return True
    except ClientError as e:
        if e.response.get("Error", {}).get("Code") in {"404", "NoSuchKey", "NotFound"}:
            return False
        raise


def download(keys, outdir, dry_run=False):
    # Mirror S3 folder structure under outdir
    # e.g., outdir/CONUS/netcdf/CHRTOUT/2017/201701010000.CHRTOUT_DOMAIN1
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    found = downloaded = 0
    missing = []

    for key in keys:
        if exists(s3, key):
            found += 1
            dst = os.path.join(outdir, key)  # preserve tree
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if dry_run:
                print(f"[FOUND] s3://{BUCKET}/{key} -> {dst} (dry-run)")
            else:
                print(f"[DOWN ] s3://{BUCKET}/{key} -> {dst}")
                s3.download_file(BUCKET, key, dst)
                downloaded += 1
        else:
            print(f"[MISS ] s3://{BUCKET}/{key}")
            missing.append(key)

    print("\nSummary:")
    print(f"  Found: {found}")
    print(f"  Downloaded: {downloaded}{' (dry-run)' if dry_run else ''}")
    if missing:
        print("  Missing:")
        for k in missing:
            print(f"    - s3://{BUCKET}/{k}")


def build_timestamps(start: str, end: str, step_hours: int = 1):
    """
    Build an inclusive list of timestamps between start and end at step_hours.
    Timestamps are YYYYMMDDHHMM (e.g., '201701010000').
    """
    t0 = datetime.strptime(start, "%Y%m%d%H%M")
    t1 = datetime.strptime(end, "%Y%m%d%H%M")
    if t1 < t0:
        raise ValueError("end must be >= start")
    if step_hours <= 0:
        raise ValueError("step_hours must be positive")
    out = []
    cur = t0
    step = timedelta(hours=step_hours)
    while cur <= t1:
        out.append(cur.strftime("%Y%m%d%H%M"))
        cur += step
    return out


def normalize_types(file_types):
    if isinstance(file_types, str):
        return ALL_TYPES if file_types.upper() == "ALL" else [file_types]
    # list passed from argparse
    return ALL_TYPES if (len(file_types) == 1 and file_types[0].upper() == "ALL") else file_types


def fetch_nwm_snapshots(
    timestamps, file_types="CHRTOUT", region="CONUS", fmt="netcdf",
    outdir=None, dry_run=False,
):
    """
    Download NWM retrospective files for a list of timestamps.
    """
    file_types = normalize_types(file_types)
    if outdir is None:
        outdir = f"nwm_retro_{timestamps[0]}"

    keys = [key_for(region, fmt, ftype, ts) for ts in timestamps for ftype in file_types]
    print(f"Timestamps: {', '.join(timestamps)}")
    print(f"Types     : {', '.join(file_types)}")
    print(f"Region/fmt: {region}/{fmt}")
    print(f"Output dir: {outdir}\n")
    download(keys, outdir, dry_run=dry_run)


def main():
    """
    Command-line interface to download NWM retrospective snapshots.
    python nwm_retro_dl.py -t 201701010000 201701011200 -k ALL
    python nwm_retro_dl.py --start 201701010000 --end 201701012300 --step 1 -k CHRTOUT RTOUT
    """
    ap = argparse.ArgumentParser(description="Download NWM v3 retrospective snapshots.")
    # Either explicit timestamps OR start/end (with step); at least one option set is required.
    ap.add_argument("-t", "--timestamps", nargs="+",
                    help="Explicit timestamps (YYYYMMDDHHMM), e.g., 201701010000 201701011200")
    ap.add_argument("--start", help="Range start timestamp (YYYYMMDDHHMM)")
    ap.add_argument("--end", help="Range end timestamp (YYYYMMDDHHMM)")
    ap.add_argument("--step", type=int, default=1, help="Range step in hours (default: 1)")

    ap.add_argument("-k", "--types", nargs="+", default=["CHRTOUT"],
                    help="File types (e.g., CHRTOUT RTOUT or ALL)")
    ap.add_argument("-r", "--region", default="CONUS")
    ap.add_argument("-f", "--fmt", default="netcdf")
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("-n", "--dry-run", action="store_true")
    args = ap.parse_args()

    # Decide timestamps source
    if args.start and args.end:
        timestamps = build_timestamps(args.start, args.end, args.step)
    elif args.timestamps:
        timestamps = args.timestamps
    else:
        ap.error("Provide either --timestamps or --start/--end (with optional --step).")

    fetch_nwm_snapshots(
        timestamps=timestamps,
        file_types=args.types,
        region=args.region,
        fmt=args.fmt,
        outdir=args.outdir,
        dry_run=args.dry_run
    )


def sample_download_call():
    """
    Example call to fetch_nwm_snapshots()
    Specify timestamps explicitly;
        timestamps = ["201612010000", "201612010100", "201612010200"]
    , or generate a range:
        timestamps = build_timestamps("201612010000", "201701010400", step_hours=1)
    """

    timestamps = build_timestamps("201601010000", "201701010000", step_hours=1)

    fetch_nwm_snapshots(
        timestamps=timestamps,
        file_types=["CHRTOUT"],  # ["LDASOUT"],  # or "ALL"
        outdir="/sciclone/schism10/feiye/STOFS3D-v8/NWM/",
        dry_run=False
    )


if __name__ == "__main__":
    # see usage of cmd line interface in main()
    # main()

    sample_download_call()

    print("Done.")
