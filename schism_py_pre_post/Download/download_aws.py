"Download files and folders from AWS S3 bucket."

import boto3
import os
from botocore import UNSIGNED
from botocore.client import Config

def download_s3_folder(bucket, prefix, local_dir):
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    paginator = s3.get_paginator('list_objects_v2')
    pages = paginator.paginate(Bucket=bucket, Prefix=prefix)

    for page in pages:
        for obj in page.get('Contents', []):
            key = obj['Key']
            relative_path = key[len(prefix):].lstrip('/')  # remove leading slash here!
            if not relative_path:
                continue
            local_file_path = os.path.join(local_dir, relative_path)
            local_file_dir = os.path.dirname(local_file_path)
            os.makedirs(local_file_dir, exist_ok=True)
            print(f"Downloading s3://{bucket}/{key} to {local_file_path}")
            s3.download_file(bucket, key, local_file_path)

if __name__ == "__main__":
    download_s3_folder(
        bucket="noaa-nos-stofs3d-pds",
        prefix="STOFS-3D-Atl/para/tmp_dir/stofs.v7.2.0_Cycle_20250625/",
        local_dir="/sciclone/schism10/feiye/TEMP/Biased/"
    )
