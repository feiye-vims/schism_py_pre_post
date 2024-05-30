# read a tif file
import rasterio

with rasterio.open('/sciclone/schism10/feiye/Alaska/DEM/akw_Yukon_BaBpEnsCu1_w84_mllw_32m_pu_-180_180_LMSL.tif') as src:
    # Read the raster data as an array
    data = src.read()

    # Create a transformation that subtracts 180 from the longitude
    # assuming the longitude is encoded in a straightforward manner
    # This might need adjustment depending on the CRS and data format
    transform = src.transform
    transform = list(transform)
    transform[2] -= 180
    transform = tuple(transform)

    # Update metadata with new transform
    new_meta = src.meta.copy()
    new_meta.update({"transform": transform})

    # Write the modified raster data to a new file
    with rasterio.open('path_to_your_output_file.tif', 'w', **new_meta) as dst:
        dst.write(data)

pass