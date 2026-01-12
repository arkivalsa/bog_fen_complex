Files for the remote sensing analysis of vegetation cover in West Siberia Lowlands, Arkhangelsk and Hudson Bay Lowlands:
Data files needed for the vegetation cover analysis, are the ESRI landcover files for 2024, openly accessible in this website: https://livingatlas.arcgis.com/landcoverexplorer/. The tiles used for each one of the three locations are the following:
WSL: 42V_20240101-20241231, 43V_20240101-20241231, 44V_20240101-20241231
Arkhangelsk: 38W_20240101-20241231
HBL: 16U_20240101-20241231, 17U_20240101-20241231, 18U_20240101-20241231
   (Bog-fen masks after analysis are in folder bog_fen_masks --> large model output files are archived on Zenodo due to size constraints
and are not stored in this GitHub repository.)

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
