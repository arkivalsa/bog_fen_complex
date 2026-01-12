# Code to create Figure 4 of main results 
# Analysis with small peat domains
import os
import rasterio
from rasterio.windows import from_bounds

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.lines as mlines

from scipy.ndimage import sobel, gaussian_filter
from scipy.ndimage import generate_binary_structure
from scipy.ndimage import distance_transform_edt

# Load the functions file
import plotting_functions

# Paths for landcover, two cropped area dems and 
save_path = r"C:\your_path\dem_landcover_analysis\small_dems"
dem1_path = r"C:\your_path\peat_dems\WSL\location1\cropped_dem_1.tif"
dem2_path = r"C:\your_path\peat_dems\WSL\location1\cropped_dem_3.tif"
bog_fen_mask = r"C:\your_path\LandCover_processed\bog_fen_mask.tif"

# Step 1:
# Load the first dem 
with rasterio.open(dem1_path) as src_crop:
    cropped_dem = src_crop.read(1)
    print("Min:", np.nanmin(cropped_dem))
    print("Max:", np.nanmax(cropped_dem))
    cropped_res = src_crop.res
    cropped_profile = src_crop.profile
    cropped_transform = src_crop.transform
    cropped_crs = src_crop.crs
    cropped_nodata = src_crop.nodata
    cropped_bounds = src_crop.bounds

    # Build a mask: True where valid
    valid_mask = ~np.isnan(cropped_dem) if cropped_nodata is None else (cropped_dem != cropped_nodata)
    
cropped_dem_masked = np.ma.masked_equal(cropped_dem, src_crop.nodata)
cropped_dem_nodata = np.where(cropped_dem == cropped_nodata, np.nan, cropped_dem)

# Load the second dem 
with rasterio.open(dem2_path) as src_crop2:
    cropped_dem2 = src_crop2.read(1)
    print("Min:", np.nanmin(cropped_dem2))
    print("Max:", np.nanmax(cropped_dem2))
    cropped_res2 = src_crop2.res
    cropped_profile2 = src_crop2.profile
    cropped_transform2 = src_crop2.transform
    cropped_crs2 = src_crop2.crs
    cropped_nodata2 = src_crop2.nodata
    cropped_bounds2 = src_crop2.bounds

    # Build a mask: True where valid
    valid_mask2 = ~np.isnan(cropped_dem2) if cropped_nodata2 is None else (cropped_dem2 != cropped_nodata2)
    
cropped_dem_masked2 = np.ma.masked_equal(cropped_dem2, src_crop2.nodata)
cropped_dem_nodata2 = np.where(cropped_dem2 == cropped_nodata2, np.nan, cropped_dem2)

# Step 2:
# Define vegetation class values and labels
# Create colormap and normalization for vegetation cover maps
class_values = [0, 1, 2, 3]  
class_names = ['Fen', 'Margin', 'Bog', 'Water']
colors = ['forestgreen', 'darkorange', 'magenta', 'aqua']  # Match the order of class_values margin:'mediumseagreen'

cmap = mcolors.ListedColormap(colors)
bounds = [v - 0.5 for v in class_values] + [class_values[-1] + 0.5]  # [−0.5, 0.5, 1.5, 2.5, 3.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)
legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, class_names)]

# Open bog_fen mask for dem1
with rasterio.open(bog_fen_mask) as src_full:
    full_transform = src_full.transform
    full_crs = src_full.crs

    # Ensure the CRS matches (reproject first if they don't)
    assert full_crs == cropped_crs, "CRS mismatch! Reproject first."

    # Get window in full raster corresponding to cropped DEM bounds
    window = from_bounds(*cropped_bounds, transform=full_transform)

    # Read that window
    cropped_mask_from_full = src_full.read(1, window=window)
    full_profile = src_full.profile

print("Full mask CRS:", full_crs)
print("Cropped DEM CRS:", cropped_crs)

# Open bog_fen mask for dem2
with rasterio.open(bog_fen_mask) as src_full2:
    full_transform2 = src_full2.transform
    full_crs2 = src_full2.crs

    # Ensure the CRS matches (reproject first if they don't)
    assert full_crs2 == cropped_crs2, "CRS mismatch! Reproject first."

    # Get window in full raster corresponding to cropped DEM bounds
    window2 = from_bounds(*cropped_bounds2, transform=full_transform2)

    # Read that window
    cropped_mask_from_full2 = src_full2.read(1, window=window2)
    full_profile2 = src_full2.profile

# Step 3:
# Bog-fen class values only where the DEM has valid data
masked_data = np.where(valid_mask, cropped_mask_from_full, np.nan)
masked_data2 = np.where(valid_mask2, cropped_mask_from_full2, np.nan)
###

# Reclassify bog_fen mask so it is not patchy
veg_map_masked = masked_data.copy()
veg_map_masked[veg_map_masked == 255] = np.nan  # 255 is the nan value assigned from the preparation steps

veg_map_masked2 = masked_data2.copy()
veg_map_masked2[veg_map_masked2 == 255] = np.nan  # 255 is the nan value assigned from the preparation steps

structure = generate_binary_structure(2, 2)  # 8-connectivity
neighbors_reclassified = plotting_functions.reclassify_small_patches_by_size_order(veg_map_masked, 4000)
small_patches = plotting_functions.count_remaining_small_patches(neighbors_reclassified, 4000, structure)

neighbors_reclassified2 = plotting_functions.reclassify_small_patches_by_size_order(veg_map_masked2, 4000)
small_patches2 = plotting_functions.count_remaining_small_patches(neighbors_reclassified2, 4000, structure)

# Copy the vegetation classification arrays to use further 
plot_map = neighbors_reclassified.copy()
plot_map[np.isnan(plot_map)] = -1  # use a placeholder

plot_map2 = neighbors_reclassified2.copy()
plot_map2[np.isnan(plot_map2)] = -1  # use a placeholder

# Step 4:
# Calculate bog-fen edge and distance from this boundary
bog_mask = np.isin(plot_map, [2, 1]) # returns a True - False array to calculate edges and distances 
# Detect edges
gy, gx = np.gradient(bog_mask.astype(float))
# Edge mask --> the bog-fen boundary 
edge_mask = np.sqrt(gx**2 + gy**2) > 0

distance_to_bog = distance_transform_edt(~edge_mask, sampling=cropped_res)
distance_to_bog[np.isnan(plot_map)] = np.nan                          # do not use the areas out of the map 
distance_to_bog[bog_mask] *= -1                                       # divide the distance towards the bog and towards the fen

bog_mask2 = np.isin(plot_map2, [2, 1]) # returns a True - False array to calculate edges and distances 
# Detect edges
gy2, gx2 = np.gradient(bog_mask2.astype(float))
# Edge mask --> the bog-fen boundary 
edge_mask2 = np.sqrt(gx2**2 + gy2**2) > 0

distance_to_bog2 = distance_transform_edt(~edge_mask2, sampling=cropped_res2)
distance_to_bog2[np.isnan(plot_map2)] = np.nan                          # do not use the areas out of the map 
distance_to_bog2[bog_mask2] *= -1 

# Step 5: Hydrological analysis 
dx, dy = cropped_res # or from dem2 they should have same resolution from preparation steps

fill_dem = plotting_functions.fill_depressions(cropped_dem_nodata)
flow_direction = plotting_functions.calculate_flow_direction_d8(cropped_dem_nodata, dx, dy)
flow_accumulation = plotting_functions.calculate_flow_accumulation(flow_direction)
flow_acc_log = np.log1p(flow_accumulation)
flow_acc_norm = (flow_acc_log - flow_acc_log.min()) / (flow_acc_log.max() - flow_acc_log.min())

fill_dem2 = plotting_functions.fill_depressions(cropped_dem_nodata2)
flow_direction2 = plotting_functions.calculate_flow_direction_d8(cropped_dem_nodata2, dx, dy)
flow_accumulation2 = plotting_functions.calculate_flow_accumulation(flow_direction2)
flow_acc_log2 = np.log1p(flow_accumulation2)
flow_acc_norm2 = (flow_acc_log2 - flow_acc_log2.min()) / (flow_acc_log2.max() - flow_acc_log2.min())

# Step 6: Wetness proxy
h_max = np.nanmax(cropped_dem_nodata)
topo_relief = h_max - cropped_dem_nodata
topo_norm = (topo_relief - np.nanmin(topo_relief)) / (np.nanmax(topo_relief) - np.nanmin(topo_relief))
wetness_proxy = (topo_norm** 0.8) * (flow_acc_norm** 0.2)

h_max2 = np.nanmax(cropped_dem_nodata2)
topo_relief2 = h_max2 - cropped_dem_nodata2
topo_norm2 = (topo_relief2 - np.nanmin(topo_relief2)) / (np.nanmax(topo_relief2) - np.nanmin(topo_relief2))
wetness_proxy2 = (topo_norm2** 0.8) * (flow_acc_norm2** 0.2)

# Step 7.1:
# Plot distance from the bog-fen boundary with wetness and bog, fen densities
# Prepare for dem1
# Flatten arrays for plotting
distances_flat = distance_to_bog.flatten()
elevation_flat = cropped_dem_nodata.flatten()
wetness_flat = wetness_proxy.flatten()
veg_flat = plot_map.flatten()

# Remove NaN values (if any)
valid_mask = ~np.isnan(distances_flat) & ~np.isnan(elevation_flat) & ~np.isnan(wetness_flat)
# Apply mask
x = distances_flat[valid_mask]
y = elevation_flat[valid_mask]
C = wetness_flat[valid_mask]

## Prepare the vegetation data 
veg_flat_valid = veg_flat[valid_mask]
# Create masks for each class
bog_mask_flat = (veg_flat_valid == 2)
margin_mask_flat = (veg_flat_valid == 1)
fen_mask_flat = (veg_flat_valid == 0)
# Apply masks to get elevation values for each class
elev_bog = y[bog_mask_flat]
elev_margin = y[margin_mask_flat]
elev_fen = y[fen_mask_flat]

# Step 7.2:
# Prepare for dem2 
# Flatten arrays for plotting
distances_flat2 = distance_to_bog2.flatten()
elevation_flat2 = cropped_dem_nodata2.flatten()
wetness_flat2 = wetness_proxy2.flatten()
veg_flat2 = plot_map2.flatten()

# Remove NaN values (if any)
valid_mask2 = ~np.isnan(distances_flat2) & ~np.isnan(elevation_flat2) & ~np.isnan(wetness_flat2)
# Apply mask
# veg = veg_flat[valid_mask]
x2 = distances_flat2[valid_mask2]
y2 = elevation_flat2[valid_mask2]
C2 = wetness_flat2[valid_mask2]

## Prepare the vegetation data 
veg_flat_valid2 = veg_flat2[valid_mask2]
# Create masks for each class
bog_mask_flat2 = (veg_flat_valid2 == 2)
margin_mask_flat2 = (veg_flat_valid2 == 1)
fen_mask_flat2 = (veg_flat_valid2 == 0)
# Apply masks to get elevation values for each class
elev_bog2 = y2[bog_mask_flat2]
elev_margin2 = y2[margin_mask_flat2]
elev_fen2 = y2[fen_mask_flat2]

# Step 8:
# Clean elevation masked arrays from nan or inf before passing to the gaussian kdes for figure 4
# for dem1
elev_fen_clean = plotting_functions.clean_array(elev_fen)
elev_bog_clean = plotting_functions.clean_array(elev_bog)
elev_margin_clean = plotting_functions.clean_array(elev_margin)

# for dem2
elev_fen_clean2 = plotting_functions.clean_array(elev_fen2)
elev_bog_clean2 = plotting_functions.clean_array(elev_bog2)
elev_margin_clean2 = plotting_functions.clean_array(elev_margin2)

# Step 9:
### FIGURE  4 ###
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde

fig = plt.figure(figsize=(12, 7.5))
outer_gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1.2, 2.8], hspace=0.15, wspace=0.12)
# -------------------------------
# TOP ROW - DEM1
# -------------------------------
## Column 1: DEM1 base map + veg
ax_map1 = fig.add_subplot(outer_gs[0, 0])
ax_map1.imshow(cropped_dem_masked, cmap='Greys', interpolation='none')
ax_map1.imshow(plot_map, cmap=cmap, norm=norm, alpha=0.5)
ax_map1.axis('off')

## Column 2: DEM1 hexbin + histograms
inner_gs1 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer_gs[0, 1], width_ratios=[4.5, 0.9, 0.8, 0.15], wspace=0.1)

ax_hex1 = fig.add_subplot(inner_gs1[0])
hb1 = ax_hex1.hexbin(x, y, C=C, reduce_C_function=np.mean, gridsize=200, cmap='RdYlBu_r', mincnt=10)
ax_hex1.axvline(0, color="black", linestyle="--")
# ax_hex1.set_xticks([0, 20000, 40000, 60000])
# ax_hex1.set_xlabel("Distance [m]")
# ax_hex1.set_ylabel("Elevation [m]")
ax_hex1.set_title("Elevation [m]", fontsize=14)
ax_hex1.grid(True, alpha=0.2)

ax_histy1 = fig.add_subplot(inner_gs1[1], sharey=ax_hex1)

# KDE histogram
# Fen KDE
combined_clean = np.concatenate([elev_fen_clean, elev_bog_clean, elev_margin_clean])
y_vals1 = np.linspace(np.min(combined_clean), np.max(combined_clean), 500)
# y_vals1 = np.linspace(min(elevation_flat), max(elevation_flat), 500)
elev_kde_fen = gaussian_kde(elev_fen_clean)
elev_x_fen = elev_kde_fen(y_vals1)
# Bog KDE
elev_kde_bog = gaussian_kde(elev_bog_clean)
elev_x_bog = elev_kde_bog(y_vals1)
# Margin KDE
elev_kde_margin = gaussian_kde(elev_margin_clean)
elev_x_margin = elev_kde_margin(y_vals1)
# Plot on the same axis (horizontal orientation)
ax_histy1.plot(elev_x_bog, y_vals1, color='magenta', lw=2, label="Bog")
ax_histy1.plot(elev_x_margin, y_vals1, color='darkorange', lw=2, label="Margin")
ax_histy1.plot(elev_x_fen, y_vals1, color='forestgreen', lw=2, label="Fen")
ax_histy1.tick_params(labelleft=False)

# -------------------------------
# BOTTOM ROW - DEM2
# -------------------------------

## Column 1: DEM2 base map + veg
ax_map2 = fig.add_subplot(outer_gs[1, 0])
ax_map2.imshow(cropped_dem_masked2, cmap='Greys', interpolation='none')
ax_map2.imshow(plot_map2, cmap=cmap, norm=norm, alpha=0.5)
ax_map2.axis('off')

from matplotlib.patches import Patch
legend_patches = [
    Patch(facecolor='forestgreen', label='Fen'),
    Patch(facecolor='orange', label='Margin'),
    Patch(facecolor='magenta', label='Bog')
]

ax_map2.legend(
    handles=legend_patches,
    loc='upper center',
    bbox_to_anchor=(0.55, -0.05),
    ncol=3,  # You have 3 patches: Fen, Margin, Bog
    frameon=False,
    fontsize=14,
    handletextpad=0.5,   # spacing between patch and label
    columnspacing=0.8    # spacing between columns
)

## Column 2: DEM2 hexbin + histograms
inner_gs2 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer_gs[1, 1], width_ratios=[4.5, 0.9, 0.8, 0.15], wspace=0.1)

ax_hex2 = fig.add_subplot(inner_gs2[0])
hb2 = ax_hex2.hexbin(x2, y2, C=C2, reduce_C_function=np.mean, gridsize=200, cmap='RdYlBu_r', mincnt=10)
ax_hex2.axvline(0, color="black", linestyle="--", label="Bog-Fen Boundary") #, label="Bog-Fen Boundary"
ax_hex2.set_xlabel("Distance from Bog-Fen Boundary [m]", fontsize=14)
ax_hex2.legend(loc='lower left', framealpha=0, fontsize=12)
ax_hex2.set_xticks([-15000, -5000, 0, 5000, 15000])
# ax_hex2.set_ylabel("Elevation [m]")
ax_hex2.grid(True, alpha=0.2)

# Add inset colorbar below legend
cb_ax = inset_axes(
    ax_hex2,
    width="40%",  # width relative to parent axes
    height="3%",  # height in %
    loc='upper right',
    bbox_to_anchor=(-0.38, -0.73, 0.90, 1),  # fine-tune position (x, y, width, height)
    bbox_transform=ax_hex2.transAxes,
    borderpad=0
)

cbar = plt.colorbar(hb2, cax=cb_ax, orientation="horizontal")
cbar.ax.xaxis.set_label_position('top')  # Move label to top
cbar.set_label("Wetness index", fontsize=12)
cbar.ax.tick_params(labelsize=10)

ax_histy2 = fig.add_subplot(inner_gs2[1], sharey=ax_hex2)

# KDE histogram
# Fen KDE
combined_clean2 = np.concatenate([elev_fen_clean2, elev_bog_clean2, elev_margin_clean2])
y_vals2 = np.linspace(np.min(combined_clean2), np.max(combined_clean2), 500)
# y_vals2 = np.linspace(min(elevation_flat2), max(elevation_flat2), 500)
elev_kde_fen2 = gaussian_kde(elev_fen_clean2)
elev_x_fen2 = elev_kde_fen2(y_vals2)
# Bog KDE
elev_kde_bog2 = gaussian_kde(elev_bog_clean2)
elev_x_bog2 = elev_kde_bog2(y_vals2)
# Margin KDE
elev_kde_margin2 = gaussian_kde(elev_margin_clean2)
elev_x_margin2 = elev_kde_margin2(y_vals2)
# Plot on the same axis (horizontal orientation)
ax_histy2.plot(elev_x_bog2, y_vals2, color='magenta', lw=2, label="Bog")
ax_histy2.plot(elev_x_margin2, y_vals2, color='darkorange', lw=2, label="Margin")
ax_histy2.plot(elev_x_fen2, y_vals2, color='forestgreen', lw=2, label="Fen")
ax_histy2.tick_params(labelleft=False)
ax_histy2.set_xlabel("Density", fontsize=12)

for ax in [ax_hex1, ax_hex2]:
    ax.tick_params(axis='both', labelsize=12)

for ax in [ax_histy1, ax_histy2]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.tick_params(axis='both', labelsize=12)
    ax.tick_params(axis='both', which='both', labelleft=False, labelbottom=False)

plt.tight_layout()
plt.savefig(os.path.join(save_path, "name_file.png"), dpi=300, bbox_inches='tight')
plt.show()