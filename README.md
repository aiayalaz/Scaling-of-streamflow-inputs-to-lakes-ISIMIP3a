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
   - Lines 377-382: lakes selection
   - Lines 383-384: method for selection of the upstreamd grids (approach Ib)
   - Lines 392-407: grids occupied by the lake
   - Lines 409-410: upstream grids
   - Lines 412-418: number of grids to pick up
   - Lines 420-520: calculation of qin
        - Lines 421-437: upstream grids <= grids occupied by the lake (approach Ia)
        - Lines 439-514: upstream grids > grids occupied by the lake (approach Ib): Lines 460-481: method==slope, Lines 482-497: method==all.
   - Lines 522-608: validation (observations) - plots and performance metrics (KGE and KGE components)
- Lines 609-1313: Approach II and validation (HYPE outputs and observations)
   - Lines 612-844: lake Vänern
        - Lines 618-632: grids occupied by the lake
        - Lines 634-635: upstream grids
        - Lines 637-645: number of grids to pick up
        - Lines 647-656: land/water area of each grid
        - Lines 659-676: qin for boundary grids
        - Lines 678-694: qin for inner grids
        - Lines 696-706: total qin
        - Lines: 715-839: validation (HYPE outputs and observations) - plots and performance metrics (KGE and KGE components)
   - Lines 846-1082: lake Vättern
        - Lines 852-866: grids occupied by the lake
        - Lines 868-869: upstream grids
        - Lines 871-879: number of grids to pick up
        - Lines 881-890: land/water area of each grid
        - Lines 893-906: qin for boundary grids
        - Lines 908-928: qin for inner grids
        - Lines 927-940: total qin
        - Lines 946-1077: validation (HYPE outputs and observations) - plots and performance metrics (KGE and KGE components)
   - Lines 1084-1313: lake Mälaren
        - Lines 1090-1104: grids occupied by the lake
        - Lines 1106-1107: upstream grids
        - Lines 1109-1117: number of grids to pick up
        - Lines 1119-1128: land/water area of each grid
        - Lines 1131-1146: qin for boundary grids
        - Lines 1148-1163: qin for inner grids
        - Lines 1165-1178: total qin
        - Lines 1184-1308: validation (HYPE outputs and observations) - plots and performance metrics (KGE and KGE components)
  
