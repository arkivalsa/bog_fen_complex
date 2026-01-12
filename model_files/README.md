Files for the model simulations (model_files):
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
Large model output files (.npz files) are archived on Zenodo due to size constraints and are not stored in this GitHub repository.

Examples of arrays used for the sensitivity analysis: 
Hin_arr = np.array([0.000005, 0.00001, 0.000015, 0.00002, 0.000025, 0.00003, 0.000035, 0.00004, 0.000045, 0.00005])
Oin_arr = np.array([0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014, 0.0016, 0.0018, 0.002])
rR_arr = np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2])
Oout_arr = np.array([0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0014, 0.0016, 0.0018, 0.002])
QR_arr = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])

