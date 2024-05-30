import os
from pathlib import Path
from glob import glob

from concurrent.futures import ProcessPoolExecutor
import multiprocessing

from osgeo import gdal
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling


fnames = sorted(glob('/sciclone/schism10/Hgrid_projects/DEMs/CONED_2022/original/*.tif'))
outdir = '/sciclone/schism10/Hgrid_projects/DEMs/CONED_2022/lonlat/'


def reproject_to_lonlat(input_tif):
    """
    Reproject a GeoTIFF file to the longitude and latitude coordinate system (EPSG:4326).

    Parameters:
    input_tif (str): Path to the input GeoTIFF file.
    """
    output_tif = f"{outdir}{Path(input_tif).stem}_lonlat.tif"

    # Open the source GeoTIFF file
    with rasterio.open(input_tif) as src:
        # Define the target CRS (EPSG:4326 for WGS84)
        dst_crs = 'EPSG:4326'
        
        # Calculate the transform and dimensions of the reprojected raster
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        
        # Define the metadata for the output file
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        
        # Open the destination file and perform the reprojection
        with rasterio.open(output_tif, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest
                )
    
    print(f"Reprojection complete. Output saved to {output_tif}")

# Example usage
# reproject_to_lonlat(fnames[0], f"{outdir}{Path(fnames[0]).stem}_lonlat.tif")
with multiprocessing.Pool(processes=os.cpu_count()) as pool:
    pool.map(reproject_to_lonlat, fnames)
