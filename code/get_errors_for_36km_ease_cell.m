function err = get_errors_for_36km_ease_cell(err_type,day_or_night,...
    sentinel_res,max_vwc,ease_row,ease_col)
%GET_ERRORS_FOR_36KM_EASE_CELL Return error values at a specific grid cell
%for EOF superresolution method from "Low Rank Gap-Filling and Downscaling 
%for SMAP Soil Moisture Datasets" by Beale, Bras, and Romberg. Ecohydrology
%2025.
%   INPUTS:
%       err_type     = 'RMSD','ubRMSD','MAD','R'
%       day_or_night = 'day' or 'night'
%       sentinel_res = 1 or 3, for 1km or 3km resolution results
%       max_vwc      = 3 or 5, for max VWC (kg/m^2) of pixels allowed
%       ease_row     = 36km EASE grid cell row index, in range 1:406
%       ease_col     = 36km EASE grid cell column index, in range 1:964
%   OUTPUT:
%       error = the error value requested
param_str = [num2str(sentinel_res),'km',num2str(max_vwc),'VWC/'];
GLOBAL_RESULTS_DIR = ['./data/Global_Error_Data/',param_str];

% load global results struct
load_struct    = load([GLOBAL_RESULTS_DIR,'global_results']);
global_results = load_struct.global_results;

% read out error
err_type_str = ['avg_',err_type,'_36km'];
err = global_results.images.(day_or_night).SR.(err_type_str)(ease_row,ease_col);
end

