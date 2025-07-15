# Scaling-of-streamflow-inputs-to-lakes-ISIMIP3a
Gridded water fluxes simulated by the global hydrological model WaterGAP 2.2e, as part of the standard model evaluation experiment (obsclim_histsoc_default) under ISIMIP Phase 3a, were used to estimate streamflow inputs to lakes.

Input data (available at: https://doi.org/10.5281/zenodo.15917845):
- WaterGAP 2.2e outputs (for each 0.5° × 0.5° grid cell):
   - Total runoff (qtot) [kg m⁻² s⁻¹] (surface + subsurface)
   - Groundwater runoff (qg) [kg m⁻² s⁻¹]
   - Discharge (dis) [m³ s⁻¹]
- River routing data: Flow direction and slope for each 0.5° × 0.5° grid cell.
- Lake information:
   - Lake area [km²]
   - Lake catchment area [m²]
   - Lake centroid coordinates [latitude, longitude]
 
The streamflow into the lake is calculated based on total (surface + subsurface) runoff, qtot [m s-1], and groundwater runoff, qg [m s-1], scaled by the lake catchment area. 
The lake catchment is delineated using flow direction data. Grid cells area classified into levels, and the selection of grid cells belonging to the catchment is based on the ratio between catchment area and the lake grid area (N).

Approach Ia: For lakes with N<=1, only the grid cells occupied partially by the lake are considered in the calculation.

Approach Ib: For lakes with N>1, both grid cells partially occupied by the lake and upstream grid cells are considered in the calculation. If more upstream grid cells are available than needed, those with the steepest slopes are prioritized. Alternatively, all available upstream grid cells can be included in the calculation.

Approach II: For lakes with lake area bigger than the lake grid cell, streamflow is calculated using both river discharge (dis [m³ s⁻¹]) from upstream grid cells adjacent to the lake the qtot [m s-1] and qg [m s-1] proportional to the land area of the grid cells partially occupied by the lake. 

Script outputs at available at: https://doi.org/10.5281/zenodo.15919494 and include monthly and yearly streamflow inputs to 71 lakes in Sweden, covering the period from January 1, 1901, to December 31, 2019. Validation plots and performance metrics are also included.
