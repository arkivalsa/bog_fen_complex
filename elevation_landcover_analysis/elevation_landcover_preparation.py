import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
# from rasterio.merge import merge
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from rasterio.mask import mask
from rasterio.merge import merge
from shapely.geometry import box
import geopandas as gpd

import plotting_functions
###

dem_path = r"C:\your_path\peat_dems\WSL"

# Step 1:
# Reproject the MERIT DEMS to EPSG:3857 from EPSG:4326 -- comment out when completed
# plotting_functions.ReprojectDEM(input_file=os.path.join(dem_path, "n65e040_dem.tif"),output_file=os.path.join(dem_path, "n65e040_dem_reproj.tif"),output_crs="EPSG:3857")

# Step 2:
# For WSL only merged the DEMS to one with the following snippet
# Load the DEMs and name the new merged file
with rasterio.open(os.path.join(dem_path, "merged_dem.tif")) as src:
    merit = src.read(1)  # Read the first band
    
    # Initial assessment of the tif file
    # Get cell size (resolution)
    cell_size = src.res  # (x resolution, y resolution)
    print(f"Cell Size (Resolution): {cell_size} (meters per pixel)")
    
    # Get map extent (bounding box)
    bounds = src.bounds  # (minx, miny, maxx, maxy)
    print(f"Map Extent: {bounds}")
    
    # Get map size (number of rows and columns)
    width = src.width
    height = src.height
    print(f"Map Size: {width} x {height} (columns x rows)")
    
    # Coordinate system
    crs = src.crs
    print(f"Coordinate Reference System (CRS): {crs}")
    
    merit[merit == src.nodata] = None  # Mask NoData values if present
    extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top] 

merged_dem_path = dem_path = r"C:\Users\archo\Documents\NIOZ\Peats\peat_dems\WSL\merged_dem.tif"
# Step 3:
# Comment out when completed
# Define the bounding box (xmin, ymin, xmax, ymax) in EPSG:3857 coordinates
# Adjust to the DEM under examination 
# WSL - replace with appropriate
crop_extents = [box(8400000, 7900000, 8500000, 8000000), 
                box(8400000, 7800000, 8500000, 7900000),
                box(8200000, 8000000, 8300000, 8100000),
                box(8400000, 8100000, 8500000, 8200000),
                box(8100000, 7900000, 8200000, 8000000),
                box(8000000, 8000000, 8100000, 8100000)]

for i, crop_extent in enumerate(crop_extents):
    with rasterio.open(os.path.join(dem_path, "merged_dem.tif")) as src:
        geo = gpd.GeoDataFrame({'geometry': [crop_extent]}, crs=src.crs)  # Convert to geodataframe
        out_image, out_transform = mask(src, geo.geometry, crop=True)
        out_meta = src.meta.copy()
        print(f"Raster bounds: {src.bounds}")
        print(f"Crop bounding box: {crop_extent.bounds}")
        
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],  # Adjust height
            "width": out_image.shape[2],   # Adjust width
            "transform": out_transform
        })
        # Remove empty dimensions (Rasterio returns a 3D array)
        # out_image = out_image[0]
        # print(np.isnan(out_image).sum())  # Counts NaNs in the output

        # out_image = np.ma.masked_equal(out_image, -9999)
        
        # Save the cropped raster
        output_file = os.path.join(dem_path, f"cropped_dem_{i+1}.tif")
        
        with rasterio.open(output_file, "w", **out_meta) as dest:
            dest.write(out_image)

        print(f"Saved cropped DEM {i+1} to {output_file}")

# Step 4.1:
# Prepare the landcover files and comment out when completed
# For WSL we merge landcover files - for rest of locations only reproject and align with elevation
# WSL and HBL landcover merge
landcover_files = sorted(glob.glob(r"C:\your_path\LandCover_Sentinel2\WSL_2024\*.tif"))
merged_landcover_path = r"C:\your_path\LandCover_processed\merged_landcover.tif"
processed_landcover_path = r"C:\your_path\LandCover_processed\final_landcover.tif"
temp_reprojected_dir = r"C:\your_path\LandCover_processed\reprojected"
os.makedirs(temp_reprojected_dir, exist_ok=True)
# We reproject all landcover maps to the DEM CRS:
reprojected_tiles = []
    
for i, lc_path in enumerate(landcover_files):
    with rasterio.open(lc_path) as src:
        dst_crs = "EPSG:3857"
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        dst_path = os.path.join(temp_reprojected_dir, f"reprojected_tile_{i}.tif")
        with rasterio.open(dst_path, 'w', **kwargs) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest  # Important for categorical landcover
            )
        reprojected_tiles.append(dst_path)

# Only for WSL and HBL
# Merge the reprojected landcover tiles 
srcs = [rasterio.open(p) for p in reprojected_tiles]
merged_array, merged_transform = merge(srcs)

merged_meta = srcs[0].meta.copy()
merged_meta.update({
    "height": merged_array.shape[1],
    "width": merged_array.shape[2],
    "transform": merged_transform
})

with rasterio.open(merged_landcover_path, 'w', **merged_meta) as dest:
    dest.write(merged_array)

# Step 4.2
# For ALL landcover files     
# Align the merged landcover with the DEM
# Change the paths accordingly for all locations - for Arkhangelsk there is no merging so it is the reprojected dem file
with rasterio.open(merged_dem_path) as dem_src:
    with rasterio.open(merged_landcover_path) as lc_src:
        profile = lc_src.meta.copy()
        profile.update({
            'crs': dem_src.crs,
            'transform': dem_src.transform,
            'width': dem_src.width,
            'height': dem_src.height,
            'dtype': lc_src.dtypes[0],
            'nodata': lc_src.nodata
        })

        with rasterio.open(processed_landcover_path, 'w', **profile) as dst:
            reproject(
                source=rasterio.band(lc_src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=lc_src.transform,
                src_crs=lc_src.crs,
                dst_transform=dem_src.transform,
                dst_crs=dem_src.crs,
                resampling=Resampling.nearest
            )

print("Landcover merged and aligned to DEM!")

# For all locations
# Create  bog-fen mask - keeping necessary classes
with rasterio.open(processed_landcover_path) as cover:
    print("CRS:", cover.crs)
    print("Transform:", cover.transform)
    print("Bounds:", cover.bounds)
    landuse = cover.read(1)  # Read the raster band
    profile = cover.profile  # Save profile for future writing
    
bog_fen_mask = np.where(landuse == 4, 2,                # bog → 2
                  np.where(landuse == 2, 0,             # fen → 0
                  np.where(landuse == 11, 1,            # margin → 1
                  np.where(landuse == 1, 3, 255))))     # water → 3, else → 255 

# Convert to uint8
bog_fen_mask = bog_fen_mask.astype(np.uint8)

# Optional: Confirm unique values
classes = np.unique(bog_fen_mask)
print("Saved bog-fen mask classes:", classes)

# Update profile for saving
mask_profile = profile.copy()
mask_profile.update({
    'dtype': 'uint8',
    'nodata': 255  # Optional: Set nodata value if needed
})

# Save the new mask
output_mask_path = r"C:\your_path\LandCover_processed\bog_fen_mask.tif"
with rasterio.open(output_mask_path, 'w', **mask_profile) as dst:
    dst.write(bog_fen_mask, 1)

print(f"Mask saved to: {output_mask_path}")