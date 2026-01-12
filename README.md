Information for the files and scripts in the repository that belong to the simulations and analysis implemented for the article with title: "Hydrological regime shifts across space as drivers of peatland landscape complexity". The version archived on Zenodo corresponds to the model and analysis used to produce the results presented in the manuscript (https://doi.org/10.5281/zenodo.18063614).

1. Files for the model simulations (model_files):
clPy.bog_fen_peatland_012026.ipynb -- The model script with the sensitivity analysis set-up, and summary plotting of the model variables. The script can run the fen wetland, or the bog-fen complex simulation depending on the settings (total runtime for the model shown in the manuscript is EndTime = 100000). In the script r (variable) is for bog vegetation and d (variable) for fen vegetation.

Hydrofunctions_iPy.cl       -- Definitions of derivatives for diffusion, advection, pressure and flow, and of boundary conditions.

plot_functions.py           -- Help module to plot the model figures shown in Results and do the hydro-topographic analysis

variables_.svg              -- Snapshot of last timestep of the model run for main model variables, for fen wetland

variables_statistics_.svg   -- Summary statistics for main model variables, for fen wetland 

variables_bog_fen_.svg              -- Snapshot of last timestep of the model run for main model variables, for bog-fen complex 

variables_statistics_bog_fen_.svg   -- Summary statistics for main model variables, for bog-fen complex 

FIG2A_model_transect_schematic.svg -- Figure 2a, model transect

FIG2B_model_variables.svg          -- Figure 2b, main model variables at last timestep

FIG3A_classes_distance.svg         -- Figure 3a, vegetation classes and distance from the bog-fen boundary

FIG3B_distance_wetness_with_histograms_082025_kde.png -- Figure 3b, elevation and velocity, plotted over distance from bog-fen boundary, highlighted by wetness index.

Fen_1024_v12b_012026.npz -- last time step state of fen simulation (main variables)

BogFenLandscape_1024_v12b_012026.npz -- last time step state of bog-fen complex simulation (main variables)

Examples of arrays used for the sensitivity analysis: 
Hin_arr = np.array([0.000005, 0.00001, 0.000015, 0.00002, 0.000025, 0.00003, 0.000035, 0.00004, 0.000045, 0.00005])
Oin_arr = np.array([0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014, 0.0016, 0.0018, 0.002])
rR_arr = np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2])
Oout_arr = np.array([0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014, 0.0016, 0.0018, 0.002])
QR_arr = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])

2. In folder elevation_landcover_files; 
   Files for the remote sensing analysis of vegetation cover in West Siberia Lowlands, Arkhangelsk and Hudson Bay Lowlands:
Data files needed for the vegetation cover analysis, are the ESRI landcover files for 2024, openly accessible in this website: https://livingatlas.arcgis.com/landcoverexplorer/. The tiles used for each one of the three locations are the following:
WSL: 42V_20240101-20241231, 43V_20240101-20241231, 44V_20240101-20241231
Arkhangelsk: 38W_20240101-20241231
HBL: 16U_20240101-20241231, 17U_20240101-20241231, 18U_20240101-20241231
   (Bog-fen masks after analysis are in folder bog_fen_masks)

   Files for the remote sensing analysis of elevation changes in West Siberia Lowlands, Arkhangelsk and Hudson Bay Lowlands:
Data files needed for the elevation analysis, are the MERIT DEMS, openly accessible in this website: https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/. The tiles used for each of the three peat complexes are the following: 
WSL: n55e070_dem, n55e075_dem, n60e070_dem, n65e070_dem, n65e075_dem
Arkhangelsk: n60e040_dem
HBL: n50w080_dem 

All examined locations were divided into smaller areas, the bounding boxes ((xmin, ymin, xmax, ymax) in EPSG:3857) of each plotted peat map (with a reference to the figure) are given below (folder, KML_files):
A. West Siberia Lowlands
area 1 : (8400000, 7900000, 8500000, 8000000) - Figure 4  
area 2 : (8400000, 7800000, 8500000, 7900000) - Figure 4
area 3 : (8200000, 8000000, 8300000, 8100000) - Figure S6A
area 4 : (8400000, 8100000, 8500000, 8200000) - Figure S6A
area 5 : (8100000, 7900000, 8200000, 8000000) - Figure S6A
area 6 : (8000000, 8000000, 8100000, 8100000) - Figure S6A
 
B. Arkhangelsk
area 1 : (4900000, 9750000, 5000000, 9850000) - Figure S6B
area 2 : (4900000, 9650000, 5000000, 9750000) - Figure S6B

C. Hudson Bay Lowlands
area 1 : (-9350000, 6700000, -9250000, 6800000) - Figure S6B
area 2 : (-9250000, 6700000, -9150000, 6800000) - Figure S6B

Code files for elevation and landcover analysis:
elevation_landcover_preparation.py -- Preparation of elevation and landcover files. Reproject files, merge where needed and prepare the bog-fen mask with necessary classes. 
elevation_landcover_analysis.py    -- Plotting Figure 4 and panels of S6A, S6B, after loading two cropped DEM areas of each location
plotting_functions.py              -- Help module for the preparation and analysis of Figures 4 and S6A, S6B
