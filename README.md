# Scaling-of-streamflow-inputs-to-lakes-ISIMIP3a
Gridded water fluxes simulated by the global hydrological model WaterGAP 2.2e, as part of the standard model evaluation experiment (obsclim_histsoc_default) under ISIMIP Phase 3a, were used to estimate streamflow inputs to lakes.

Input data (available at: https://doi.org/10.5281/zenodo.17588875):
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

Script outputs at available at: https://doi.org/10.5281/zenodo.17588905 and include monthly and yearly streamflow inputs to 70 lakes in Sweden, covering the period from January 1, 1901, to December 31, 2019. Validation plots and performance metrics are also included.

Comments to the main_script.R:
- Line 21: Load helper functions avalaible at helper_functions_script.R 
- Lines 23-80: Load input data
- Lines 82-373: Approach I(a/b) and validation (HYPE outputs)
   - Lines: 85-99: lakes selection
   - Lines 100-101: method for selection of the upstreamd grids (approach Ib)
   - Lines 112-127: grids occupied by the lake
   - Lines 129-130: upstream grids
   - Lines 132-138: number of grids to pick up
   - Lines 140-240: calculation of qin
        - Lines 141-157: upstream grids <= grids occupied by the lake (approach Ia)
        - Lines 159-233: upstream grids > grids occupied by the lake (approach Ib): Lines 180-201: method==slope, Lines 202-217: method==all.
   - Lines 242-367: validation (HYPE outputs) - plots and performance metrics (KGE and KGE components)
- Lines 374-608: Approach I(a/b) and validation (observations)
   - Lines 377-386: lakes selection
   - Lines 387-388: method for selection of the upstreamd grids (approach Ib)
   - Lines 396-411: grids occupied by the lake
   - Lines 413-414: upstream grids
   - Lines 416-422: number of grids to pick up
   - Lines 424-524: calculation of qin
        - Lines 425-441: upstream grids <= grids occupied by the lake (approach Ia)
        - Lines 443-518: upstream grids > grids occupied by the lake (approach Ib): Lines 464-485: method==slope, Lines 486-501: method==all.
   - Lines 526-642: validation (observations) - plots and performance metrics (KGE and KGE components)
- Lines 644-1348: Approach II and validation (HYPE outputs and observations)
   - Lines 647-879: lake Vänern
        - Lines 653-667: grids occupied by the lake
        - Lines 669-670: upstream grids
        - Lines 672-680: number of grids to pick up
        - Lines 682-691: land/water area of each grid
        - Lines 694-711: qin for boundary grids
        - Lines 713-729: qin for inner grids
        - Lines 731-742: total qin
        - Lines: 750-874: validation (HYPE outputs and observations) - plots and performance metrics (KGE and KGE components)
   - Lines 881-1117: lake Vättern
        - Lines 887-901: grids occupied by the lake
        - Lines 903-904: upstream grids
        - Lines 906-914: number of grids to pick up
        - Lines 916-925: land/water area of each grid
        - Lines 928-941: qin for boundary grids
        - Lines 943-960: qin for inner grids
        - Lines 962-973: total qin
        - Lines 981-1112: validation (HYPE outputs and observations) - plots and performance metrics (KGE and KGE components)
   - Lines 1119-1348: lake Mälaren
        - Lines 1125-1139: grids occupied by the lake
        - Lines 1141-1142: upstream grids
        - Lines 1144-1152: number of grids to pick up
        - Lines 1154-1163: land/water area of each grid
        - Lines 1166-1181: qin for boundary grids
        - Lines 1183-1198: qin for inner grids
        - Lines 1200-1211: total qin
        - Lines 1219-1343: validation (HYPE outputs and observations) - plots and performance metrics (KGE and KGE components)
  
