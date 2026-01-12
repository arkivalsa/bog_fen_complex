import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from scipy.ndimage import label, generate_binary_structure, binary_dilation
from collections import Counter
from rasterio.windows import from_bounds
from scipy.stats import gaussian_kde
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

# Function to reproject to EPSG:3857 from EPSG:4326
def ReprojectDEM(input_file,output_file,output_crs,target_resolution=130):

    with rasterio.open(input_file) as src:
        transform, width, height = rasterio.warp.calculate_default_transform(
            src.crs, output_crs, src.width, src.height, *src.bounds,
            resolution=target_resolution)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': output_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
    
        with rasterio.open(output_file, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                rasterio.warp.reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=output_crs,
                    resampling=rasterio.warp.Resampling.nearest)

structure = generate_binary_structure(2, 2)  # 8-connectivity
def reclassify_small_patches_by_size_order(veg_map, min_size):
    output_map = veg_map.copy()
    original_map = veg_map.copy()
    classes = [c for c in np.unique(original_map) if c != 255]

    # Collect all small patches
    patch_info_list = []

    for c in classes:
        binary_map = (original_map == c)
        labeled_array, num_features = label(binary_map, structure=structure)

        for patch_id in range(1, num_features + 1):
            patch_mask = (labeled_array == patch_id)
            patch_size = np.sum(patch_mask)

            if patch_size < min_size:
                patch_info_list.append({
                    'class': c,
                    'size': patch_size,
                    'mask': patch_mask
                })

    # Sort patches by size (largest first)
    patch_info_list.sort(key=lambda x: -x['size'])

    # Apply reclassification from largest to smallest
    for patch_info in patch_info_list:
        patch_mask = patch_info['mask']
        neighbors = get_neighbor_classes(output_map, patch_mask, structure)

        if neighbors:
            majority_class = Counter(neighbors).most_common(1)[0][0]
            if majority_class != patch_info['class']:
                output_map[patch_mask] = majority_class
        else:
            output_map[patch_mask] = 0  # fallback to background

    return output_map

def get_neighbor_classes(arr, patch_mask, structure, max_iterations=5):
    current_mask = patch_mask.copy()
    for _ in range(max_iterations):
        dilated = binary_dilation(current_mask, structure=structure)
        neighbors_mask = dilated & (~current_mask)
        neighbor_classes = arr[neighbors_mask]
        neighbor_classes = neighbor_classes[(neighbor_classes != 255)]

        if neighbor_classes.size > 0:
            return neighbor_classes.tolist()
        current_mask = dilated

    return []

# Track remaining small patches
def count_remaining_small_patches(veg_map, min_size, structure):
    from scipy.ndimage import label
    small_patch_count = 0
    for c in np.unique(veg_map):
        if c == 255:
            continue
        binary = veg_map == c
        labeled, num = label(binary, structure=structure)
        for pid in range(1, num + 1):
            if np.sum(labeled == pid) < min_size:
                small_patch_count += 1
    print(f"Remaining small patches: {small_patch_count}")
   
def calculate_flow_direction_d8(dem, x_res, y_res):
    # Define the 8 possible neighbor directions
    directions = np.array([
        [-1, -1], [-1, 0], [-1, 1],
        [ 0, -1],          [ 0, 1],
        [ 1, -1], [ 1, 0], [ 1, 1]
    ])
    
    # Flow direction matrix
    flow_direction = np.zeros_like(dem, dtype=int)
    
    # Iterate over the internal cells of the DEM (ignore boundaries)
    for i in range(1, dem.shape[0] - 1):
        for j in range(1, dem.shape[1] - 1):
            max_diff = 0  # Initialize with zero, we want to find the max descent
            best_dir = 0  # Initialize flow direction
            
            # Iterate over the 8 neighbors
            for d in range(directions.shape[0]):
                di, dj = directions[d]
                neighbor_elevation = dem[i + di, j + dj]
                elevation_diff = dem[i, j] - neighbor_elevation  # Positive for descent
                
                # Update if this neighbor has a steeper descent
                if elevation_diff > max_diff:
                    max_diff = elevation_diff
                    best_dir = d + 1  # Directions are indexed from 1 to 8
            
            # Assign the best direction (steepest descent)
            flow_direction[i, j] = best_dir
    
    return flow_direction

from collections import deque
def calculate_flow_accumulation(flow_direction): 
    acc = np.ones_like(flow_direction)
    rows, cols = flow_direction.shape
    directions = [
        (-1, -1), (-1, 0), (-1, 1),
        ( 0, -1),          ( 0, 1),
        ( 1, -1), ( 1, 0), ( 1, 1)
    ]
    
    # Count how many cells flow into each cell
    inflow_count = np.zeros_like(flow_direction, dtype=int)
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            d = flow_direction[i, j] - 1
            if 0 <= d < 8:
                di, dj = directions[d]
                ni, nj = i + di, j + dj
                inflow_count[ni, nj] += 1

    # Start with cells that have no inflows
    queue = deque([(i, j) for i in range(rows) for j in range(cols) if inflow_count[i, j] == 0])

    while queue:
        i, j = queue.popleft()
        d = flow_direction[i, j] - 1
        if 0 <= d < 8:
            di, dj = directions[d]
            ni, nj = i + di, j + dj
            acc[ni, nj] += acc[i, j]
            inflow_count[ni, nj] -= 1
            if inflow_count[ni, nj] == 0:
                queue.append((ni, nj))
    
    return acc

def fill_depressions(dem):
    # Get the dimensions of the DEM
    rows, cols = dem.shape

    # Create a copy of the DEM to modify
    filled_dem = np.copy(dem)

    # Define the 8 possible neighbor directions (vertical, horizontal, diagonal)
    neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    # Iterate until no changes are made
    change = True
    while change:
        change = False
        for i in range(1, rows-1):
            for j in range(1, cols-1):
                # Get the current cell's elevation
                current_value = filled_dem[i, j]
                
                # Get the values of the neighbors
                neighbor_values = [filled_dem[i + di, j + dj] for di, dj in neighbors]
                
                # Check if the current cell is a depression
                if current_value < min(neighbor_values):
                    # Fill the depression by raising it to the level of the lowest neighbor
                    filled_dem[i, j] = min(neighbor_values)
                    change = True  # A change was made, so we need to iterate again
    
    return filled_dem

# --- Utility: Extract elevations for each class from DEM + vegetation raster ---
def extract_elevations(dem_path, veg_path, class_map):
    with rasterio.open(dem_path) as dem_src, rasterio.open(veg_path) as veg_src:
        # Ensure CRS match
        if dem_src.crs != veg_src.crs:
            raise ValueError("CRS mismatch! Reproject vegetation raster first.")

        # Ensure same resolution/extent: read veg raster window covering DEM
        window = from_bounds(*dem_src.bounds, transform=veg_src.transform)
        veg = veg_src.read(1, window=window, out_shape=(dem_src.height, dem_src.width))

        # Read DEM
        dem = dem_src.read(1)

        # Handle nodata
        dem = np.where(dem == dem_src.nodata, np.nan, dem)
        veg = np.where(veg == 255, np.nan, veg)  # if 255 used for nodata

    # Extract elevations by class
    elev_dict = {}
    for cls_val, cls_name in class_map.items():
        elevs = dem[veg == cls_val]
        # Remove NaN and infs
        elevs = elevs[np.isfinite(elevs)]
        elev_dict[cls_name] = elevs

    return elev_dict

def density_overlap(arr1, arr2, grid=None):
    arr1 = arr1[np.isfinite(arr1)]
    arr2 = arr2[np.isfinite(arr2)]
    if len(arr1) < 2 or len(arr2) < 2:
        return np.nan

    kde1 = gaussian_kde(arr1)
    kde2 = gaussian_kde(arr2)

    # Evaluation grid
    if grid is None:
        lo = min(arr1.min(), arr2.min())
        hi = max(arr1.max(), arr2.max())
        grid = np.linspace(lo, hi, 500)

    d1 = kde1(grid)
    d2 = kde2(grid)

    overlap = np.trapz(np.minimum(d1, d2), grid)
    return overlap

def clean_array(arr):
    return arr[np.isfinite(arr)]