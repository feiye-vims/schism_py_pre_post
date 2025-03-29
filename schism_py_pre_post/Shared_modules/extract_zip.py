'''
Processing SWOT data
'''

import os
import zipfile
from glob import glob

import pandas as pd
import geopandas as gpd
from tqdm import tqdm


def extract_zip_files(folder_path, extracted_folder):
    """
    Extract all ZIP files in a folder
    """
    os.makedirs(extracted_folder, exist_ok=True)
    zip_files = glob(os.path.join(folder_path, '*.zip'))
    for i, zip_file in enumerate(tqdm(zip_files)):
        if i < 4245:
            continue
        try:
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                zip_ref.extractall(extracted_folder)
        except Exception as e:
            print(f"Error extracting {zip_file}: {e}")


def merge_shapefiles(folder_path):
    """
    Merge all shapefiles in a folder into a single GeoDataFrame
    """
    shapefiles = glob(os.path.join(folder_path, '**', '*.shp'), recursive=True)
    merged_gdf = gpd.GeoDataFrame(
        pd.concat([gpd.read_file(shp) for shp in shapefiles], ignore_index=True))
    merged_gdf.to_file(os.path.join(folder_path, 'merged_shapefile.shp'))


if __name__ == '__main__':
    FOLDER_PATH = '/sciclone/schism10/Hgrid_projects/SWOT/Zip/'
    EXTRACTED_FOLDER = '/sciclone/schism10/Hgrid_projects/SWOT/Extracted/'

    extract_zip_files(FOLDER_PATH, EXTRACTED_FOLDER)
    merge_shapefiles(EXTRACTED_FOLDER)
    print('Done')
