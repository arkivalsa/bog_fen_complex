import os
from os import path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from scipy.ndimage import distance_transform_edt, binary_dilation
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.stats import gaussian_kde
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import gridspec

# Colormaps
# Load the existing colormap
YlGn = plt.get_cmap("YlGn")
YlGn_array = YlGn(np.linspace(0, 1, 256))  # Get 256 color values
YlGn_array[:128, -1] = np.linspace(0, 1, 128)  # Gradually increase opacity
new_YlGn = mcolors.ListedColormap(YlGn_array)

OrRd = plt.get_cmap("spring_r")
OrRd_array = OrRd(np.linspace(0, 1, 128))  # Get 256 color values
OrRd_array[:64, -1] = np.linspace(0, 1, 64)  # Gradually increase opacity
new_OrRd = mcolors.ListedColormap(OrRd_array)

def bog_fen_edge_mask(rs, ds, FinalCount):
    # Binary mask to plot borders
    bog_mask = (rs[:,:,FinalCount] > 0.7)
    margin_mask = (rs[:,:,FinalCount] > 0.01) & (ds[:,:,FinalCount] > 0.01)
    fen_mask = ~(bog_mask | margin_mask) & (ds[:,:,FinalCount] > 0.01)
    
    combined_mask = np.zeros_like(rs[:, :, FinalCount], dtype=int)
    combined_mask[bog_mask] = 0
    combined_mask[margin_mask] = 1
    combined_mask[fen_mask] = 2

    bog_margin_mask = bog_mask | margin_mask

    # Detect edges
    gy, gx = np.gradient(bog_margin_mask.astype(float))
    # Edge mask
    edge_mask = np.sqrt(gx**2 + gy**2) > 0
    return edge_mask, combined_mask, bog_margin_mask, bog_mask, margin_mask, fen_mask # bog_fen_edge_mask(rs, ds, FinalCount)[0], [1], [2], [3], [4], [5]

def transect_plot(results, FinalCount, dX, save_path, std=False):
    ss, sos, hs,  rs, ds = results  # Unpack the results tuple
    QR       = 0.1
    channel_mask = (QR/(QR+hs[:,:,FinalCount]) < 0.8)
    edge_mask = bog_fen_edge_mask(rs, ds, FinalCount)[0]
    
    sediment_transect = (sos[:, :, FinalCount] + ss[:, :, FinalCount])[700, 510:610] # [row, col start : col end]
    water_transect = (hs[:, :, FinalCount])[700, 510:610]
    channel_mask_transect = channel_mask[700, 510:610]
    bog_transect = (rs[:, :, FinalCount])[700, 510:610]
    fen_transect = (ds[:, :, FinalCount])[700, 510:610]

    length = np.arange(510, 610) * dX

    # Find boundary along the transect
    boundary_row = edge_mask[700, 510:610]  # edge_mask --> the bog-fen line
    boundary_cols = np.where(boundary_row > 0.5)[0]
    boundary_positions = (510 + boundary_cols) * dX

    # Boolean masks for vegetation zones
    bog_mask_transect = bog_transect > 0.01  # or some other threshold if needed
    fen_mask_transect = fen_transect > 0.85

    # Slightly above sediment line for visibility
    offset = 0.06  # vertical offset for the vegetation ribbons

    plt.figure(figsize=(9.2, 2.8))
    # Bog and fen areas
    plt.fill_between(length, sediment_transect + offset, sediment_transect + offset + 0.04,
                     where=fen_mask_transect, color='forestgreen', alpha=0.6 , label=None)
    plt.fill_between(length, sediment_transect + offset, sediment_transect + offset + 0.04,
                     where=bog_mask_transect, color='magenta', alpha=0.6 , label=None)
    # Sediment and water levels
    plt.plot(length, sediment_transect, linewidth=2.5, color="saddlebrown", label="Peat surface [m]")

    mask = (channel_mask_transect) & (water_transect > sediment_transect)
    plt.fill_between(length,sediment_transect,water_transect,where=mask,interpolate=True,color='skyblue',alpha=0.7,label='Water in Creek')

    # Bog–fen boundary
    if boundary_cols.size > 0:
        boundary_x = (510 + boundary_cols[0]) * dX
        plt.axvline(x=boundary_x, linestyle='--', color='k')
        plt.text(boundary_x - 1.15, 2.3, 'Bog–Fen\nBoundary',
             fontsize=12, va='center', ha='left', fontweight='bold')

    # Rainfall over bog
    plt.annotate('Rainfall',
                 xy=(52, 1.3), xytext=(52, 1.5),
                 arrowprops=dict(arrowstyle='->', color='black'),
                 fontsize=12, color='black', fontweight='bold',
                 ha='center')

    # Rainfall over fen
    plt.annotate('Rainfall',
                 xy=(58, 0.4), xytext=(58, 0.6),
                 arrowprops=dict(arrowstyle='->', color='black'),
                 fontsize=12, color='black', fontweight='bold',
                 ha='center')

    # Overland flow
    plt.annotate('',
                 xy=(56.5, 0.5), xytext=(55.5, 0.8),
                 arrowprops=dict(arrowstyle='->', color='black'))
    plt.text(55.5, 0.95, 'Overland \nflow',
                 fontsize=12, color='black', fontweight='bold',
                 ha='center')

    # Regional flow
    plt.annotate('',
                 xy=(60.4, 0.3), xytext=(60.8, 0.6),
                 arrowprops=dict(arrowstyle='->', color='black',connectionstyle='arc3,rad=0.3'))
    plt.text(59.2, 0.8, 'Regional water',
                 fontsize=12, color='black', fontweight='bold',
                 va='center')

    plt.xlabel("Distance along transect [m]")
    plt.grid(True, alpha=0.2)
    plt.xlim(length[0], length[-1])
    plt.ylim(0, 2.8)

    handles, labels = plt.gca().get_legend_handles_labels()
    filtered = [(h, l) for h, l in zip(handles, labels) if l in ["Peat surface [m]", "Water level [m]"]]
    handles, labels = zip(*filtered)
    plt.legend(handles, labels, loc='upper right', frameon=False)

    plt.tight_layout()
    plt.savefig(path.join(save_path, "FIG2A_model_transect_schematic.svg"), dpi=250)
    plt.show()

def model_summary(input_data, save_path):
    ds, rs, sos, ss, ns, FinalCount, Length_X, Length_Y, dX = input_data
    
    fig, axes = plt.subplots(1, 3, figsize=(10, 6))
    fig.subplots_adjust(wspace=0.10, bottom=0.3)  # Reserve space for colorbars

    # --- First panel: Fen + Bog vegetation ---
    im0 = axes[0].imshow(ds[:, :, FinalCount], cmap=new_YlGn, extent=[0, Length_X, 0, Length_Y], alpha=0.6)
    im1 = axes[0].imshow(rs[:, :, FinalCount], cmap=new_OrRd, extent=[0, Length_X, 0, Length_Y], alpha=0.6)

    # Plot transect line
    transect_y = (1024 - 700) * dX  # Convert row to y-coordinate in meters
    transect_x_start = 510 * dX
    transect_x_end = 610 * dX
    # Plot edge points
    axes[0].plot(transect_x_start, transect_y, marker='o', color='black', markersize=1)
    axes[0].plot(transect_x_end, transect_y, marker='o', color='black', markersize=1)

    # Patches colorbars for vegetation
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='forestgreen', alpha=0.6, label='Fen vegetation'),
        Patch(facecolor='magenta', alpha=0.6, label='Bog vegetation')
    ]
    # Position legend manually 
    axes[0].legend(handles=legend_elements,loc='upper center', bbox_to_anchor=(0.5, 0.0),  ncol=1, fontsize=14,frameon=False, columnspacing=0.5, handletextpad=0.6)

    # --- Second panel: Elevation ---
    im2 = axes[1].imshow(sos[:, :, FinalCount] + ss[:, :, FinalCount],
                         cmap='gist_earth', extent=[0, Length_X, 0, Length_Y], vmin=0, vmax=1.5, alpha=0.9)
    ax1_pos = axes[1].get_position()
    cax2 = fig.add_axes([ax1_pos.x0, ax1_pos.y0 - 0.06, ax1_pos.width, 0.025])
    cb2 = fig.colorbar(im2, cax=cax2, orientation='horizontal')
    cb2.set_label('Elevation [m]')

    # --- Third panel: Residence Time or Velocity ---
    # im3 = axes[2].imshow(residence_time,
    #                      cmap='Greys', extent=[0, Length_X, 0, Length_Y], vmin=0, vmax=8)
    im3 = axes[2].imshow(ns[:, :, FinalCount],
                         cmap='CMRmap', extent=[0, Length_X, 0, Length_Y], vmin=0, vmax=np.max(ns))
    ax2_pos = axes[2].get_position()
    cax3 = fig.add_axes([ax2_pos.x0, ax2_pos.y0 - 0.06, ax2_pos.width, 0.025])
    cb3 = fig.colorbar(im3, cax=cax3, orientation='horizontal')
    cb3.set_label('Flow velocity [m/s]')

    # --- Scale bar on first panel ---
    scale_length = 20
    scale_x_position = 3
    scale_y_start = 80
    scale_y_end = scale_y_start + scale_length

    axes[0].add_line(Line2D([scale_x_position, scale_x_position],
                            [scale_y_start, scale_y_end],
                            color='black', linewidth=2))
    axes[0].text(scale_x_position + 2.0,
                 scale_y_start + scale_length / 2,
                 f'{scale_length} m',
                 horizontalalignment='left',
                 fontsize=17,
                 color='black')

    # --- Hide tick labels ---
    for ax in axes:
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

    plt.savefig(path.join(save_path, "FIG2B_model_variables.svg"), dpi=250)
    
def classes_distance(input_data3, save_path):
    ds, rs, FinalCount, Length_X, Length_Y, dX = input_data3
    edge_mask = bog_fen_edge_mask(rs, ds, FinalCount)[0]
    combined_mask = bog_fen_edge_mask(rs, ds, FinalCount)[1]
    bog_margin_mask = bog_fen_edge_mask(rs, ds, FinalCount)[2]
    
    distance_to_bog = distance_transform_edt(~edge_mask) * dX
    distance_to_bog[bog_margin_mask] *= -1
    
    # Discrete colormap for Fen, Margin, Bog (0, 1, 2)
    cmap_classes = ListedColormap(['magenta', 'darkorange', 'forestgreen'])
    norm_classes = BoundaryNorm([0, 1, 2, 3], cmap_classes.N)  # 3 classes: 0–1, 1–2, 2–3

    fig, axes = plt.subplots(2, 1, figsize=(5, 7))
    fig.subplots_adjust(wspace=0.2)

    # --- Top panel: Fen/Margin/Bog classes ---
    im0 = axes[0].imshow(combined_mask, cmap=cmap_classes, norm=norm_classes, extent=[0, Length_X, 0, Length_Y], alpha=0.6)
    divider0 = make_axes_locatable(axes[0])
    cax0 = divider0.append_axes("right", size="5%", pad=0.05)

    # Separate creation and alpha assignment for the colorbar
    from matplotlib.colorbar import ColorbarBase
    # Create a colormap with alpha
    colors = ['forestgreen', 'darkorange', 'magenta']
    colors_with_alpha = []
    for c in colors:
        rgba = plt.cm.colors.to_rgba(c)
        colors_with_alpha.append((rgba[0], rgba[1], rgba[2], 0.6))  # set alpha here

    cmap_alpha = ListedColormap(colors_with_alpha)
    # Reverse for the colorbar
    cmap_alpha_reversed = cmap_alpha.reversed()

    cb0 = ColorbarBase(cax0, cmap=cmap_alpha_reversed, norm=norm_classes,
                       ticks=[0.5, 1.5, 2.5], orientation='vertical')
    cb0.set_ticklabels(['Bog', 'Margin', 'Fen'])  # Top → bottom

    # --- Bottom panel: Distance from Bog boundary ---
    im1 = axes[1].imshow(distance_to_bog, cmap='gist_yarg', origin='upper')
    divider1 = make_axes_locatable(axes[1])
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)
    cb1 = fig.colorbar(im1, cax=cax1)
    cb1.set_label('Distance from \nBog-Fen Boundary [m]')

    # Contour for boundary
    axes[1].contour(edge_mask, levels=[0.5], colors='black', linestyles='--', linewidths=1)
    
    # --- Scale bar on second panel ---
    scale_length = 200
    scale_x_position = 1003
    scale_y_start = 20
    scale_y_end = scale_y_start + scale_length

    axes[1].add_line(Line2D([scale_x_position, scale_x_position],
                            [scale_y_start, scale_y_end],
                            color='white', linewidth=2))
    axes[1].text(scale_x_position - 9.0,
                 scale_y_start + scale_length / 2,
                 '20 m',
                 horizontalalignment='right',
                 fontsize=14,
                 color='white')

    # Hide ticks if not needed
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    plt.savefig(path.join(save_path, "FIG3A_classes_distance.svg"), dpi=250)
    
## Functions for hydro-topographic controls and Figure 3b ##
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


def variables_distance(input_data4, save_path):
    ds, rs, sos, ss, ns, FinalCount, Length_X, Length_Y, dX, dY = input_data4
    
    # Input
    edge_mask = bog_fen_edge_mask(rs, ds, FinalCount)[0]
    combined_mask = bog_fen_edge_mask(rs, ds, FinalCount)[1]
    bog_margin_mask = bog_fen_edge_mask(rs, ds, FinalCount)[2]
    bog_mask = bog_fen_edge_mask(rs, ds, FinalCount)[3]
    margin_mask = bog_fen_edge_mask(rs, ds, FinalCount)[4]
    fen_mask = bog_fen_edge_mask(rs, ds, FinalCount)[5]
    
    distance_to_bog = distance_transform_edt(~edge_mask) * dX 
    distance_to_bog[bog_margin_mask] *= -1
    elevation = ss[:,:,FinalCount] + sos[:,:,FinalCount]
    velocity = ns[:,:,FinalCount]
    
    # Precalculations for hydrotopographic metrics
    flow_direction = calculate_flow_direction_d8(elevation, dX, dY)
    flow_accumulation = calculate_flow_accumulation(flow_direction)

    flow_acc_log = np.log1p(flow_accumulation)
    flow_acc_norm = (flow_acc_log - flow_acc_log.min()) / (flow_acc_log.max() - flow_acc_log.min())

    h_max = np.max(elevation)
    
    topo_relief = h_max - elevation
    topo_norm = (topo_relief - topo_relief.min()) / (topo_relief.max() - topo_relief.min())
    wetness_proxy = (topo_norm** 0.8) * (flow_acc_norm** 0.2)
    
    # Precalculations for distances
    distance_flat = distance_to_bog.flatten()
    elevation_flat = elevation.flatten()  # Total elevation
    velocity_flat = velocity.flatten()
    wetness_proxy_flat = wetness_proxy.flatten()

    # Remove NaN values (if any)
    valid_mask = ~np.isnan(distance_flat) & ~np.isnan(elevation_flat) & ~np.isnan(velocity_flat)

    distance_flat = distance_flat[valid_mask]
    elevation_flat = elevation_flat[valid_mask]
    vwetness_proxy_flat = wetness_proxy_flat[valid_mask]
    
    # To ensure correct plotting
    bog_mask_flat = bog_mask.flatten()
    elev_bog = elevation_flat[bog_mask_flat]
    vel_bog = velocity_flat[bog_mask_flat]

    margin_mask_flat = margin_mask.flatten()
    elev_margin = elevation_flat[margin_mask_flat]
    vel_margin = velocity_flat[margin_mask_flat]

    fen_mask_flat = fen_mask.flatten()
    elev_fen = elevation_flat[fen_mask_flat]
    vel_fen = velocity_flat[fen_mask_flat]
    
    # Figure 3b 
        # Layout setup: 2 main rows (for two plots), each with extra space for histograms
    fig = plt.figure(figsize=(7.2, 7.5))
    gs = gridspec.GridSpec(4, 5, height_ratios=[0.05, 6, 0.05, 6], width_ratios=[1, 1, 1, 0.7, 0.7], hspace=0.2, wspace=0.2)

    # First plot: elevation
    ax_scatter1 = fig.add_subplot(gs[1, 0:3])

    # Second plot: velocity
    ax_scatter2 = fig.add_subplot(gs[3, 0:3], sharex=ax_scatter1)

    ax_histy1 = fig.add_subplot(gs[1, 3], sharey=ax_scatter1)  # 2 columns wide
    ax_histy2 = fig.add_subplot(gs[3, 3], sharey=ax_scatter2)

    # --- First panel: Elevation ---
    hb = ax_scatter1.hexbin(distance_flat, elevation_flat,
                        C=wetness_proxy_flat, 
                        reduce_C_function=np.mean, 
                        gridsize=200, mincnt=10, cmap="RdYlBu_r") # gridsize=1200, mincnt=1
    ax_scatter1.tick_params(labelbottom=False)
    ax_scatter1.axvline(0, color="black", linestyle="--", label="Bog-Fen Boundary")
    ax_scatter1.set_ylabel("Elevation [m]")
    ax_scatter1.legend(loc='upper right', framealpha=0, fontsize=12)
    ax_scatter1.grid(True, alpha=0.2)

    # Add inset colorbar below legend
    cb_ax = inset_axes(
        ax_scatter1,
        width="40%",  # width relative to parent axes
        height="3%",  # height in %
        loc='upper right',
        bbox_to_anchor=(-0.06, -0.23, 1, 1),  # fine-tune position (x, y, width, height)
        bbox_transform=ax_scatter1.transAxes,
        borderpad=0
    )

    cbar = plt.colorbar(hb, cax=cb_ax, orientation="horizontal")
    cbar.ax.xaxis.set_label_position('top')  # Move label to top
    cbar.set_label("Wetness index", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    # KDE histogram
    # Fen KDE
    y_vals1 = np.linspace(min(elevation_flat), max(elevation_flat), 500)
    elev_kde_fen = gaussian_kde(elev_fen)
    elev_x_fen = elev_kde_fen(y_vals1)
    # Bog KDE
    elev_kde_bog = gaussian_kde(elev_bog)
    elev_x_bog = elev_kde_bog(y_vals1)
    # Margin KDE
    elev_kde_margin = gaussian_kde(elev_margin)
    elev_x_margin = elev_kde_margin(y_vals1)
    # Plot on the same axis (horizontal orientation)
    ax_histy1.plot(elev_x_bog, y_vals1, color='magenta', lw=2, label="Bog")
    ax_histy1.plot(elev_x_margin, y_vals1, color='darkorange', lw=2, label="Margin")
    ax_histy1.plot(elev_x_fen, y_vals1, color='forestgreen', lw=2, label="Fen")

    # --- Second panel: Velocity ---
    ax_scatter2.hexbin(distance_flat, velocity_flat,
                        C=wetness_proxy_flat, 
                        reduce_C_function=np.mean, 
                        gridsize=200, mincnt=10, cmap="RdYlBu_r") # gridsize=1200, mincnt=1
    ax_scatter2.axvline(0, color="black", linestyle="--", label="Bog-Fen Boundary")
    ax_scatter2.set_ylabel("Flow velocity [m/s]")
    ax_scatter2.set_xlabel("Distance from Bog-Fen Boundary [m]")
    ax_scatter2.grid(True, alpha=0.2)

    # KDE histogram
    # Fen KDE
    y_vals2 = np.linspace(min(velocity_flat), max(velocity_flat), 500)
    vel_kde_fen = gaussian_kde(vel_fen)
    vel_x_fen = vel_kde_fen(y_vals2)
    # Bog KDE
    vel_kde_bog = gaussian_kde(vel_bog)
    vel_x_bog = vel_kde_bog(y_vals2)
    # Margin KDE
    vel_kde_margin = gaussian_kde(vel_margin)
    vel_x_margin = vel_kde_margin(y_vals2)
    # Plot on the same axis (horizontal orientation)
    ax_histy2.plot(vel_x_bog, y_vals2, color='magenta', lw=2, label="Bog")
    ax_histy2.plot(vel_x_margin, y_vals2, color='darkorange', lw=2, label="Margin")
    ax_histy2.plot(vel_x_fen, y_vals2, color='forestgreen', lw=2, label="Fen")

    for ax in [ax_histy1, ax_histy2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    #     ax.tick_params(axis='both', labelsize=12)
        ax.tick_params(axis='both', which='both', labelleft=False, labelbottom=False)


    # Final layout and save
    plt.savefig(path.join(save_path, "FIG3B_distance_wetness_with_histograms_082025_kde.png"), dpi=250, bbox_inches='tight')
    plt.show()