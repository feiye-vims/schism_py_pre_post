import os
import geopandas as gpd
from concurrent.futures import ThreadPoolExecutor

download_dir = '/sciclone/schism10/Hgrid_projects/DEMs/CONED_2022/'
shp_idx = gpd.read_file('/sciclone/schism10/Hgrid_projects/DEMs/CONED_2022/tileindex_CoNED_NGOM2_DEM_2022.shp')

def download_tile(tile_info):
    tile_name = tile_info[1]['location']
    tile_url = tile_info[1]['url']

    # download the tile
    if not os.path.exists(f'{download_dir}{tile_name}'):
        os.system(f'curl -o {download_dir}{tile_name} {tile_url}')

# Use ThreadPoolExecutor to download tiles in parallel
with ThreadPoolExecutor(max_workers=10) as executor:
    executor.map(download_tile, shp_idx.iterrows())  
    pass