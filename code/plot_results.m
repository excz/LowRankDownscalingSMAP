%plot_results.m
% Plots Figures 3 and 4 from "Low Rank Gap-Filling and Downscaling for SMAP 
% Soil Moisture Datasets" by Beale, Bras, and Romberg. Ecohydrology 2025.
datapath        = './data';
dataset_label   = 'SMAPVEX16I'; % SMAPVEX16I, SMAPVEX16M, SMAPVEX15A
num_components  = -1;           % number of EOFs to use; -1 means use EOFs with R^2 >= R2_threshold
R2_threshold    = 0.1;

%% Load PALS and SMAP Data
pals                  = load_pals_data(dataset_label,datapath);
[smap,smap_composite] = load_smap_data(dataset_label,datapath);
valid_cells           = not(isnan(pals{1}.km1.vwc));
params                = get_params();

%% Load HR Means and EOF Data
[raw_data_hr_means,mc_hr_means,eof,pm] = load_hr_means_and_eof_data(datapath,dataset_label,params);

%% Compute EOF Predictions for all EASE Cells at SMAPVEX Location
pred = get_EOF_predictions_and_baselines(dataset_label,smap,raw_data_hr_means,mc_hr_means,eof,...
    pm,num_components,R2_threshold,params);

%% Plot EOF Predictions Compared to SMAP-Sentinel Data
valid_data = plot_results_against_smap_sentinel(eof,raw_data_hr_means,mc_hr_means,pm,num_components,...
    dataset_label,R2_threshold,datapath,params);

%% Plot EOF Predictions Compared to PALS Data
plot_results_against_pals(datapath,dataset_label,pals,pred);

%% Perform Synthetic Experiment of Gap-Filling Missing Pixels using EOFs
perform_EOF_interpolation_experiment(datapath,dataset_label,eof,valid_data,R2_threshold,params);

%% Ancillary Functions
%==========================================================================
function params = get_params()
%GET_PARAMS
params.MAX_SMAP_SENTINEL_OVERPASS_DIFF_DAYTIME   = 6;  % [hours]
params.MAX_SMAP_SENTINEL_OVERPASS_DIFF_NIGHTTIME = 12; % [hours]
params.MIN_PIXEL_QUALITY = 'uncertain'; % none, uncertain, recommended
%
params.REMOVE_ANOM_HIGH             = 1;    % remove anomalously high soil moisture values
params.REMOVE_DISAGG_FAIL           = 1;    % remove pixels where disaggregation failed
params.MAX_ACCEPTABLE_VWC_SMAP      = 5;    % remove SMAP Enhanced Radiometer pixels with VWC > this threshold [kg/m^2]
params.MAX_ACCEPTABLE_VWC_SENTINEL  = 3;    % remove SMAP-Sentinel pixels with VWC > this threshold [kg/m^2]
params.MAX_ACCEPTABLE_WBF           = 0.05; % remove pixels with water body fraction > this threshold
params.REMOVE_COASTAL               = 1;    % remove pixels marked as being a 'coastal area'
params.REMOVE_URBAN                 = 1;    % remove pixels marked as being 'urban area' (based on SMAP definition -- above 25 percent of pixel is urban)
params.REMOVE_FROZEN                = 1;    % remove pixels marked as being 'frozen ground'
params.REMOVE_MOUNTAINOUS           = 1;    % remove pixels marked as 'mountainous terrain' (based on threshold set by SMAP -- DEM slope STD > 3)
params.REMOVE_PRECIP                = 1;    % remove pixels where it has rained recently (based on SMAP definition -- ???)
params.REMOVE_SCENE_EDGE            = 1;    % remove pixels acquired at the edge of the Sentinel scenewhere disaggregation can be suspect
params.REMOVE_ANOM_SIGMA0           = 1;    % remove pixels were Sentinel sigma0 data in the grid cell were outside nominal expected range for the scene
params.MIN_REASONABLE_SM            = 0.02; 
params.MAX_REASONABLE_SM            = 0.8;    

% matrix completion parameters
params.MATRIX_COMPLETION_METHOD          = 'BPDN'; % 'BP' or 'BPDN'
params.MU                                = 1e-1;   % smoothing parameter for matrix completion (penalizes deviation from X0)
params.EPSILON_MC                        = 0.008;  % max L2 deviation constraint for matrix completion (average per pixel)
params.MAX_EXPECTED_RANK                 = 10;     % maximum rank we expect soil moisture images or blocks to have for any given location
params.MIN_OBSERVATION_FACTOR            = 3;      % factor we multiply by the max expected rank to get the number of observations we require for each pixel
params.MAX_EXPECTED_RANK_BLOCKS          = 5; 
params.MIN_OBSERVATION_FACTOR_BLOCKS     = 3; 
params.MIN_ACCEPTABLE_OBSERVATIONS       = round(params.MAX_EXPECTED_RANK*params.MIN_OBSERVATION_FACTOR);
params.MIN_ACCEPTABLE_OBSERVATIONS_BLOCK = round(params.MAX_EXPECTED_RANK_BLOCKS*params.MIN_OBSERVATION_FACTOR_BLOCKS);

% parameters for EOF method
params.DATASET_FIRST_VALID_DAY  = '01-Apr-2015'; % first date to include measurements from
params.DATASET_LAST_VALID_DAY   = '01-Apr-2022'; % last date to include measurements from
params.TEST_SET_START_DATE      = '29-Oct-2018'; % only measurements occurring after this date will be in test set if params.RANDOMLY_SAMPLE_DATASETS=false
params.TEST_SET_END_DATE        = '29-Oct-2019'; % only measurements occurring before this date will be in test set if params.RANDOMLY_SAMPLE_DATASETS=false
params.VALID_SET_START_DATE     = '01-Jun-2016'; % only measurements occurring after this date will be in validation set if params.RANDOMLY_SAMPLE_DATASETS=false
params.VALID_SET_END_DATE       = '31-Jul-2016'; % only measurements occurring before this date will be in validation set if params.RANDOMLY_SAMPLE_DATASETS=false
params.RANDOMLY_SAMPLE_DATASETS = true;          % true=randomly sample train/test/valid datasets, false=use specific year for test set (specify in above two lines)
params.RNG_SEED                 = 1;             % seed for random number generator (to allow ability to replicate same random sampling)
params.MIN_TRAIN_SET_SIZE       = params.MIN_ACCEPTABLE_OBSERVATIONS_BLOCK; % minimum number of observations we require in training set to even consider applying method
params.MIN_TEST_SET_SIZE        = 5;      % test set for EOF-SR error estimation
params.MIN_VALID_SET_SIZE       = 5;      % minimum validation set size we require
params.SINGULAR_VALUE_THRESHOLD = 1e-3;

% other parameters
params.EASE_GRID_RES = 36;           % EASE cell spatial resolution [km]
params.SMAP_RES      = 9;            % SMAP pixel resolution [km]
params.SENTINEL_RES  = 1;            % Sentinel pixel resolution [km]
%
params.LR_SIDE_PIXELS           = params.EASE_GRID_RES/params.SMAP_RES;
params.HR_SIDE_PIXELS           = params.EASE_GRID_RES/params.SENTINEL_RES;
params.NUM_PIXELS_PER_BLOCK     = (params.HR_SIDE_PIXELS / params.LR_SIDE_PIXELS)^2;
params.MIN_ACCEPTABLE_LR_PIXELS = 1;
params.NUM_LR_PIXELS            = params.LR_SIDE_PIXELS*params.LR_SIDE_PIXELS;
params.NUM_HR_PIXELS            = params.HR_SIDE_PIXELS*params.HR_SIDE_PIXELS;
end
%--------------------------------------------------------------------------

%==========================================================================
function pals = load_pals_data(dataset_label,datapath)
disp(['Loading PALS data for ', dataset_label])
switch dataset_label
    case 'SMAPVEX15A'
        data = load([datapath,'/PALS/pals_az.mat']);
    case 'SMAPVEX16I'
        data = load([datapath,'/PALS/pals_ia.mat']);
    case 'SMAPVEX16M'
        data = load([datapath,'/PALS/pals_ma.mat']);
end
pals = data.pals;
disp('    Done.')
end
%--------------------------------------------------------------------------

%==========================================================================
function [smap,smap_composite] = load_smap_data(dataset_label,datapath)
%LOAD_SMAP_DATA
disp(['Loading SMAP data for ', dataset_label])
switch lower(dataset_label)
    case {'smapvex16i','vex16i','16i','iowa','ia','i'}
        ease_rows = 66:67;
        ease_cols = 232:233;
        measurement_dates_day = {'26-May-2016','28-May-2016',...           % dates are in UTC
            '31-May-2016','02-Jun-2016','03-Jun-2016','05-Jun-2016',...
            '03-Aug-2016','05-Aug-2016','06-Aug-2016','08-Aug-2016',...
            '11-Aug-2016','13-Aug-2016','14-Aug-2016','16-Aug-2016',...
            '19-Aug-2016'};
        measurement_dates_night = {'27-May-2016','28-May-2016',...
            '31-May-2016','02-Jun-2016','04-Jun-2016','05-Jun-2016',...
            '03-Aug-2016','05-Aug-2016','07-Aug-2016','08-Aug-2016',...
            '11-Aug-2016','13-Aug-2016','15-Aug-2016','16-Aug-2016',...
            '19-Aug-2016'};
        rstr      = 'R0';
        cstr      = 'C';
    case {'smapvex16m','vex16m','16m','manitoba','ma','m'}
        ease_rows = 48:49;
        ease_cols = 220:221;
        measurement_dates_day = {'06-Jun-2016','08-Jun-2016',...
            '09-Jun-2016','11-Jun-2016','12-Jun-2016','14-Jun-2016',...
            '16-Jun-2016','19-Jun-2016','20-Jun-2016','22-Jun-2016',...
            '13-Jul-2016','14-Jul-2016','16-Jul-2016','18-Jul-2016',...
            '19-Jul-2016','21-Jul-2016','22-Jul-2016','24-Jul-2016'};
        measurement_dates_night = {'07-Jun-2016','08-Jun-2016',...
            '10-Jun-2016','12-Jun-2016','13-Jun-2016','15-Jun-2016',...
            '16-Jun-2016','18-Jun-2016','20-Jun-2016','21-Jun-2016',...
            '23-Jun-2016','12-Jul-2016','14-Jul-2016','15-Jul-2016',...
            '17-Jul-2016','18-Jul-2016','20-Jul-2016','22-Jul-2016',...
            '23-Jul-2016',};
        rstr      = 'R0';
        cstr      = 'C';
    case {'smapvex15a','smapvex15','vex15','15','arizona','az','a'}
        ease_rows = 96:98;
        ease_cols = 185:188;
        measurement_dates_day = {'31-Jul-2015','02-Aug-2015',...
            '03-Aug-2015','05-Aug-2015','08-Aug-2015','10-Aug-2015',...
            '11-Aug-2015','13-Aug-2015','18-Aug-2015','19-Aug-2015',...
            '21-Aug-2015'};
        measurement_dates_night = {'01-Aug-2015','02-Aug-2015',...
            '04-Aug-2015','07-Aug-2015','09-Aug-2015','10-Aug-2015',...
            '12-Aug-2015','15-Aug-2015','17-Aug-2015','18-Aug-2015',...
            '20-Aug-2015'};
        rstr      = 'R0';
        cstr      = 'C';
    otherwise
        error('invalid dataset label');
end
num_rows       = length(ease_rows);
num_cols       = length(ease_cols);
num_pixels     = 4; % 12/sentinel_res
num_row_pixels = num_rows*num_pixels;
num_col_pixels = num_cols*num_pixels;

% get number of measurements in validation experiment
num_measurements_day   = length(measurement_dates_day);
num_measurements_night = length(measurement_dates_night);

% prepare composite SMAP image struct
smap_composite.day.soil_moisture                   = cell(num_measurements_day,1);
smap_composite.day.soil_moisture_dca               = cell(num_measurements_day,1);
smap_composite.day.time                            = NaN(num_measurements_day,1);
smap_composite.day.measurement_inds                = NaN(num_measurements_day,1);
smap_composite.day.vegetation_opacity              = cell(num_measurements_day,1);
smap_composite.day.vegetation_opacity_dca          = cell(num_measurements_day,1);
smap_composite.day.vegetation_water_content        = cell(num_measurements_day,1);
smap_composite.day.surface_temperature             = cell(num_measurements_day,1);
smap_composite.day.albedo                          = cell(num_measurements_day,1);
smap_composite.day.albedo_dca                      = cell(num_measurements_day,1);
smap_composite.day.roughness                       = cell(num_measurements_day,1);
smap_composite.day.roughness_dca                   = cell(num_measurements_day,1);
smap_composite.day.surface_flags                   = cell(num_measurements_day,1);
smap_composite.day.retrieval_quality_flags         = cell(num_measurements_day,1);
smap_composite.day.retrieval_quality_flags_dca     = cell(num_measurements_day,1);
smap_composite.day.retrieval_quality_flags_merge   = cell(num_measurements_day,1);
%
smap_composite.night.soil_moisture                 = cell(num_measurements_night,1);
smap_composite.night.soil_moisture_dca             = cell(num_measurements_night,1);
smap_composite.night.time                          = NaN(num_measurements_night,1);
smap_composite.night.measurement_inds              = NaN(num_measurements_night,1);
smap_composite.night.vegetation_opacity            = cell(num_measurements_night,1);
smap_composite.night.vegetation_opacity_dca        = cell(num_measurements_night,1);
smap_composite.night.vegetation_water_content      = cell(num_measurements_night,1);
smap_composite.night.surface_temperature           = cell(num_measurements_night,1);
smap_composite.night.albedo                        = cell(num_measurements_night,1);
smap_composite.night.albedo_dca                    = cell(num_measurements_night,1);
smap_composite.night.roughness                     = cell(num_measurements_night,1);
smap_composite.night.roughness_dca                 = cell(num_measurements_night,1);
smap_composite.night.surface_flags                 = cell(num_measurements_night,1);
smap_composite.night.retrieval_quality_flags       = cell(num_measurements_night,1);
smap_composite.night.retrieval_quality_flags_dca   = cell(num_measurements_night,1);
smap_composite.night.retrieval_quality_flags_merge = cell(num_measurements_night,1);

% load data
one2eleven = @(x) x(1:11);
upscale    = @(x) kron(x,ones(9));
%
smap       = cell(num_rows,num_cols);
for i=1:num_rows
    for j=1:num_cols
        disp(['    Reading data for 36km EASE cell (',num2str(i),',',num2str(j),') of a ',num2str(num_rows), 'x', num2str(num_cols), ' array...'])
        
        % read out data for next 36km EASE cell
        dat = load([datapath,'/SMAP/',rstr,num2str(ease_rows(i)),cstr,num2str(ease_cols(j))]);
        
        % extract day measurements indices
        measurement_inds_day = NaN(num_measurements_day,1);
        for k=1:num_measurements_day
            % get indices corresponding to measurement date
            ind = find(cellfun(@(x) strcmp(one2eleven(datestr(x)),measurement_dates_day{k}),dat.smap.observation_time));
            if isempty(ind)
                % no observation at this time for this cell
                continue;
            elseif length(ind)==1
                if strcmp(dat.smap.ascending_or_descending{ind},'D')
                    measurement_inds_day(k) = ind;
                end
            elseif length(ind)==2
                first_is_ascending   = strcmp(dat.smap.ascending_or_descending{ind(1)},'A');
                second_is_ascending  = strcmp(dat.smap.ascending_or_descending{ind(2)},'A');
                first_is_descending  = strcmp(dat.smap.ascending_or_descending{ind(1)},'D');
                second_is_descending = strcmp(dat.smap.ascending_or_descending{ind(2)},'D');
                if first_is_ascending&&second_is_ascending
                    error('both measurements are ascending');
                elseif first_is_descending&&second_is_descending
                    error('both measurements are descending');
                else
                    if first_is_descending&&second_is_ascending
                        measurement_inds_day(k) = ind(1);
                    elseif first_is_ascending&&second_is_descending
                        measurement_inds_day(k) = ind(2);
                    end
                end
            end
        end
        
        % extract night measurements indices
        measurement_inds_night = NaN(num_measurements_night,1);
        for k=1:num_measurements_night
            % get indices corresponding to measurement date
            ind = find(cellfun(@(x) strcmp(one2eleven(datestr(x)),measurement_dates_night{k}),dat.smap.observation_time));
            if isempty(ind)
                % no observation at this time for this cell
                continue;
            elseif length(ind)==1
                if strcmp(dat.smap.ascending_or_descending{ind},'A')
                    measurement_inds_night(k) = ind;
                end
            elseif length(ind)==2
                first_is_ascending   = strcmp(dat.smap.ascending_or_descending{ind(1)},'A');
                second_is_ascending  = strcmp(dat.smap.ascending_or_descending{ind(2)},'A');
                first_is_descending  = strcmp(dat.smap.ascending_or_descending{ind(1)},'D');
                second_is_descending = strcmp(dat.smap.ascending_or_descending{ind(2)},'D');
                if first_is_ascending&&second_is_ascending
                    error('both measurements are ascending');
                elseif first_is_descending&&second_is_descending
                    error('both measurements are descending');
                else
                    if first_is_descending&&second_is_ascending
                        measurement_inds_night(k) = ind(2);
                    elseif first_is_ascending&&second_is_descending
                        measurement_inds_night(k) = ind(1);
                    end
                end
            end
        end
        
        % put day measurements into output data struct
        smap{i,j}.day.soil_moisture                 = cell(num_measurements_day,1);
        smap{i,j}.day.soil_moisture_dca             = cell(num_measurements_day,1);
        smap{i,j}.day.time                          = NaN(num_measurements_day,1);
        smap{i,j}.day.measurement_inds              = NaN(num_measurements_day,1);
        smap{i,j}.day.vegetation_opacity            = cell(num_measurements_day,1);
        smap{i,j}.day.vegetation_opacity_dca        = cell(num_measurements_day,1);
        smap{i,j}.day.vegetation_water_content      = cell(num_measurements_day,1);
        smap{i,j}.day.surface_temperature           = cell(num_measurements_day,1);
        smap{i,j}.day.albedo                        = cell(num_measurements_day,1);
        smap{i,j}.day.albedo_dca                    = cell(num_measurements_day,1);
        smap{i,j}.day.roughness                     = cell(num_measurements_day,1);
        smap{i,j}.day.roughness_dca                 = cell(num_measurements_day,1);
        smap{i,j}.day.vegetation_water_content      = cell(num_measurements_day,1);
        smap{i,j}.day.water_body_fraction           = cell(num_measurements_day,1);
        smap{i,j}.day.surface_flags                 = cell(num_measurements_day,1);
        smap{i,j}.day.retrieval_quality_flags       = cell(num_measurements_day,1);
        smap{i,j}.day.retrieval_quality_flags_dca   = cell(num_measurements_day,1);
        smap{i,j}.day.retrieval_quality_flags_merge = cell(num_measurements_day,1);
        for k=1:num_measurements_day
            if isnan(measurement_inds_day(k))
                smap{i,j}.day.soil_moisture{k}                 = NaN(num_pixels,num_pixels);
                smap{i,j}.day.soil_moisture_dca{k}             = NaN(num_pixels,num_pixels);
                smap{i,j}.day.time(k)                          = NaN;
                smap{i,j}.day.measurement_inds(k)              = NaN;
                smap{i,j}.day.vegetation_opacity{k}            = NaN(num_pixels,num_pixels);
                smap{i,j}.day.vegetation_opacity_dca{k}        = NaN(num_pixels,num_pixels);
                smap{i,j}.day.vegetation_water_content{k}      = NaN(num_pixels,num_pixels);
                smap{i,j}.day.surface_temperature{k}           = NaN(num_pixels,num_pixels);
                smap{i,j}.day.albedo{k}                        = NaN(num_pixels,num_pixels);
                smap{i,j}.day.albedo_dca{k}                    = NaN(num_pixels,num_pixels);
                smap{i,j}.day.roughness{k}                     = NaN(num_pixels,num_pixels);
                smap{i,j}.day.roughness_dca{k}                 = NaN(num_pixels,num_pixels);
                smap{i,j}.day.water_body_fraction{k}           = NaN(num_pixels,num_pixels);
                smap{i,j}.day.surface_flags{k}                 = NaN(num_pixels,num_pixels);
                smap{i,j}.day.retrieval_quality_flags{k}       = NaN(num_pixels,num_pixels);
                smap{i,j}.day.retrieval_quality_flags_dca{k}   = NaN(num_pixels,num_pixels);
                smap{i,j}.day.retrieval_quality_flags_merge{k} = NaN(num_pixels,num_pixels);
            else
                smap{i,j}.day.soil_moisture{k}                 = dat.smap.soil_moisture_sca{measurement_inds_day(k)};
                smap{i,j}.day.soil_moisture_dca{k}             = dat.smap.soil_moisture{measurement_inds_day(k)};
                smap{i,j}.day.time(k)                          = dat.smap.observation_time{measurement_inds_day(k)};
                smap{i,j}.day.measurement_inds(k)              = measurement_inds_day(k);
                smap{i,j}.day.vegetation_opacity{k}            = dat.smap.vegetation_opacity_sca{measurement_inds_day(k)};
                smap{i,j}.day.vegetation_opacity_dca{k}        = dat.smap.vegetation_opacity{measurement_inds_day(k)};
                smap{i,j}.day.vegetation_water_content{k}      = dat.smap.vegetation_water_content{measurement_inds_day(k)};
                smap{i,j}.day.surface_temperature{k}           = dat.smap.surface_temperature{measurement_inds_day(k)};
                smap{i,j}.day.albedo{k}                        = dat.smap.albedo_sca{measurement_inds_day(k)};
                smap{i,j}.day.albedo_dca{k}                    = dat.smap.albedo{measurement_inds_day(k)};
                smap{i,j}.day.roughness{k}                     = dat.smap.roughness_sca{measurement_inds_day(k)};
                smap{i,j}.day.roughness_dca{k}                 = dat.smap.roughness{measurement_inds_day(k)};
                smap{i,j}.day.retrieval_quality_flags{k}       = dat.smap.retrieval_quality_flags_sca{measurement_inds_day(k)};
                smap{i,j}.day.retrieval_quality_flags_dca{k}   = dat.smap.retrieval_quality_flags{measurement_inds_day(k)};
                rqf1                                           = dat.smap.retrieval_quality_flags_sca{measurement_inds_day(k)};
                rqf2                                           = dat.smap.retrieval_quality_flags{measurement_inds_day(k)};
                x1                                             = dec2bin(rqf1(:))-'0';
                x2                                             = dec2bin(rqf2(:))-'0';
                smap{i,j}.day.retrieval_quality_flags_merge{k} = reshape(bin2dec(num2str(x1 | x2)),[4,4]);
                smap{i,j}.day.water_body_fraction{k}           = dat.smap.water_body_fraction{measurement_inds_day(k)};
                smap{i,j}.day.surface_flags{k}                 = dat.smap.surface_flags{measurement_inds_day(k)};
            end
        end
        
        % put night measurements into output data struct
        smap{i,j}.night.soil_moisture                 = cell(num_measurements_night,1);
        smap{i,j}.night.soil_moisture_dca             = cell(num_measurements_night,1);
        smap{i,j}.night.time                          = NaN(num_measurements_night,1);
        smap{i,j}.night.measurement_inds              = NaN(num_measurements_night,1);
        smap{i,j}.night.vegetation_opacity            = cell(num_measurements_night,1);
        smap{i,j}.night.vegetation_opacity_dca        = cell(num_measurements_night,1);
        smap{i,j}.night.vegetation_water_content      = cell(num_measurements_night,1);
        smap{i,j}.night.surface_temperature           = cell(num_measurements_night,1);
        smap{i,j}.night.albedo                        = cell(num_measurements_night,1);
        smap{i,j}.night.albedo_dca                    = cell(num_measurements_night,1);
        smap{i,j}.night.roughness                     = cell(num_measurements_night,1);
        smap{i,j}.night.roughness_dca                 = cell(num_measurements_night,1);
        smap{i,j}.night.water_body_fraction           = cell(num_measurements_night,1);
        smap{i,j}.night.surface_flags                 = cell(num_measurements_night,1);
        smap{i,j}.night.retrieval_quality_flags       = cell(num_measurements_night,1);
        smap{i,j}.night.retrieval_quality_flags_dca   = cell(num_measurements_night,1);
        smap{i,j}.night.retrieval_quality_flags_merge = cell(num_measurements_night,1);
        for k=1:num_measurements_night
            if isnan(measurement_inds_night(k))
                smap{i,j}.night.soil_moisture{k}                 = NaN(num_pixels,num_pixels);
                smap{i,j}.night.soil_moisture_dca{k}             = NaN(num_pixels,num_pixels);
                smap{i,j}.night.time(k)                          = NaN;
                smap{i,j}.night.measurement_inds(k)              = NaN;
                smap{i,j}.night.vegetation_opacity{k}            = NaN(num_pixels,num_pixels);
                smap{i,j}.night.vegetation_opacity_dca{k}        = NaN(num_pixels,num_pixels);
                smap{i,j}.night.vegetation_water_content{k}      = NaN(num_pixels,num_pixels);
                smap{i,j}.night.surface_temperature{k}           = NaN(num_pixels,num_pixels);
                smap{i,j}.night.albedo{k}                        = NaN(num_pixels,num_pixels);
                smap{i,j}.night.albedo_dca{k}                    = NaN(num_pixels,num_pixels);
                smap{i,j}.night.roughness{k}                     = NaN(num_pixels,num_pixels);
                smap{i,j}.night.roughness_dca{k}                 = NaN(num_pixels,num_pixels);
                smap{i,j}.night.water_body_fraction{k}           = NaN(num_pixels,num_pixels);
                smap{i,j}.night.surface_flags{k}                 = NaN(num_pixels,num_pixels);
                smap{i,j}.night.retrieval_quality_flags{k}       = NaN(num_pixels,num_pixels);
                smap{i,j}.night.retrieval_quality_flags_dca{k}   = NaN(num_pixels,num_pixels);
                smap{i,j}.night.retrieval_quality_flags_merge{k} = NaN(num_pixels,num_pixels);
            else
                smap{i,j}.night.soil_moisture{k}                 = dat.smap.soil_moisture_sca{measurement_inds_night(k)};
                smap{i,j}.night.soil_moisture_dca{k}             = dat.smap.soil_moisture{measurement_inds_night(k)};
                smap{i,j}.night.time(k)                          = dat.smap.observation_time{measurement_inds_night(k)};
                smap{i,j}.night.measurement_inds(k)              = measurement_inds_night(k);
                smap{i,j}.night.vegetation_opacity{k}            = dat.smap.vegetation_opacity_sca{measurement_inds_night(k)};
                smap{i,j}.night.vegetation_opacity_dca{k}        = dat.smap.vegetation_opacity{measurement_inds_night(k)};
                smap{i,j}.night.vegetation_water_content{k}      = dat.smap.vegetation_water_content{measurement_inds_night(k)};
                smap{i,j}.night.surface_temperature{k}           = dat.smap.surface_temperature{measurement_inds_night(k)};
                smap{i,j}.night.albedo{k}                        = dat.smap.albedo_sca{measurement_inds_night(k)};
                smap{i,j}.night.albedo_dca{k}                    = dat.smap.albedo{measurement_inds_night(k)};
                smap{i,j}.night.roughness{k}                     = dat.smap.roughness_sca{measurement_inds_night(k)};
                smap{i,j}.night.roughness_dca{k}                 = dat.smap.roughness{measurement_inds_night(k)};
                smap{i,j}.night.water_body_fraction{k}           = dat.smap.water_body_fraction{measurement_inds_night(k)};
                smap{i,j}.night.surface_flags{k}                 = dat.smap.surface_flags{measurement_inds_night(k)};
                smap{i,j}.night.retrieval_quality_flags{k}       = dat.smap.retrieval_quality_flags_sca{measurement_inds_night(k)};
                smap{i,j}.night.retrieval_quality_flags_dca{k}   = dat.smap.retrieval_quality_flags{measurement_inds_night(k)};
                %
                rqf1                                             = dat.smap.retrieval_quality_flags_sca{measurement_inds_night(k)};
                rqf2                                             = dat.smap.retrieval_quality_flags{measurement_inds_night(k)};
                x1                                               = dec2bin(rqf1(:))-'0';
                x2                                               = dec2bin(rqf2(:))-'0';
                smap{i,j}.night.retrieval_quality_flags_merge{k} = reshape(bin2dec(num2str(x1 | x2)),[4,4]);
            end
        end
        
        % put day composite measurements into composite image
        row_inds = (num_pixels*(i-1)+1):(num_pixels*i); 
        col_inds = (num_pixels*(j-1)+1):(num_pixels*j);
        rem_inds = false(num_measurements_day,1);
        for k=1:num_measurements_day
            if (i==1)&&(j==1)
                smap_composite.day.soil_moisture{k}                 = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.soil_moisture_dca{k}             = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.time(k)                          = NaN;
                smap_composite.day.measurement_inds(k)              = NaN;
                smap_composite.day.vegetation_opacity{k}            = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.vegetation_opacity_dca{k}        = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.vegetation_water_content{k}      = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.surface_temperature{k}           = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.albedo{k}                        = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.albedo_dca{k}                    = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.roughness{k}                     = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.roughness_dca{k}                 = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.water_body_fraction{k}           = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.surface_flags{k}                 = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.retrieval_quality_flags{k}       = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.retrieval_quality_flags_dca{k}   = NaN(num_row_pixels,num_col_pixels);
                smap_composite.day.retrieval_quality_flags_merge{k} = NaN(num_row_pixels,num_col_pixels);
            end
            if not(isnan(measurement_inds_day(k)))
                smap_composite.day.soil_moisture{k}(row_inds,col_inds)                 = dat.smap.soil_moisture_sca{measurement_inds_day(k)};
                smap_composite.day.soil_moisture_dca{k}(row_inds,col_inds)             = dat.smap.soil_moisture{measurement_inds_day(k)};
                smap_composite.day.time(k)                                             = dat.smap.observation_time{measurement_inds_day(k)};
                smap_composite.day.measurement_inds(k)                                 = measurement_inds_day(k);
                smap_composite.day.vegetation_opacity{k}(row_inds,col_inds)            = dat.smap.vegetation_opacity_sca{measurement_inds_day(k)};
                smap_composite.day.vegetation_opacity_dca{k}(row_inds,col_inds)        = dat.smap.vegetation_opacity{measurement_inds_day(k)};
                smap_composite.day.vegetation_water_content{k}(row_inds,col_inds)      = dat.smap.vegetation_water_content{measurement_inds_day(k)};
                smap_composite.day.surface_temperature{k}(row_inds,col_inds)           = dat.smap.surface_temperature{measurement_inds_day(k)};
                smap_composite.day.albedo{k}(row_inds,col_inds)                        = dat.smap.albedo_sca{measurement_inds_day(k)};
                smap_composite.day.albedo_dca{k}(row_inds,col_inds)                    = dat.smap.albedo{measurement_inds_day(k)};
                smap_composite.day.roughness{k}(row_inds,col_inds)                     = dat.smap.roughness_sca{measurement_inds_day(k)};
                smap_composite.day.roughness_dca{k}(row_inds,col_inds)                 = dat.smap.roughness{measurement_inds_day(k)};
                smap_composite.day.water_body_fraction{k}(row_inds,col_inds)           = dat.smap.water_body_fraction{measurement_inds_day(k)};
                smap_composite.day.surface_flags{k}(row_inds,col_inds)                 = dat.smap.surface_flags{measurement_inds_day(k)};
                rqfsca                                                                 = dat.smap.retrieval_quality_flags_sca{measurement_inds_day(k)};
                rqfdca                                                                 = dat.smap.retrieval_quality_flags{measurement_inds_day(k)};
                smap_composite.day.retrieval_quality_flags{k}(row_inds,col_inds)       = rqfsca;
                smap_composite.day.retrieval_quality_flags_dca{k}(row_inds,col_inds)   = rqfdca;
                x1                                                                     = dec2bin(rqfsca(:))-'0';
                x2                                                                     = dec2bin(rqfdca(:))-'0';
                smap_composite.day.retrieval_quality_flags_merge{k}(row_inds,col_inds) = reshape(bin2dec(num2str(x1|x2)),[num_pixels,num_pixels]);
            else
                rem_inds(k) = true;
            end
        end
        
        % put night composite measurements into composite image
        row_inds = (num_pixels*(i-1)+1):(num_pixels*i); 
        col_inds = (num_pixels*(j-1)+1):(num_pixels*j);
        rem_inds = false(num_measurements_night,1);
        for k=1:num_measurements_night
            if (i==1)&&(j==1)
                smap_composite.night.soil_moisture{k}                 = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.soil_moisture_dca{k}             = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.time(k)                          = NaN;
                smap_composite.night.measurement_inds(k)              = NaN;
                smap_composite.night.vegetation_opacity{k}            = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.vegetation_opacity_dca{k}        = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.vegetation_water_content{k}      = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.surface_temperature{k}           = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.albedo{k}                        = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.albedo_dca{k}                    = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.roughness{k}                     = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.roughness_dca{k}                 = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.water_body_fraction{k}           = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.surface_flags{k}                 = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.retrieval_quality_flags{k}       = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.retrieval_quality_flags_dca{k}   = NaN(num_row_pixels,num_col_pixels);
                smap_composite.night.retrieval_quality_flags_merge{k} = NaN(num_row_pixels,num_col_pixels);
            end
            if not(isnan(measurement_inds_night(k)))
                smap_composite.night.soil_moisture{k}(row_inds,col_inds)                 = dat.smap.soil_moisture_sca{measurement_inds_night(k)};
                smap_composite.night.soil_moisture_dca{k}(row_inds,col_inds)             = dat.smap.soil_moisture{measurement_inds_night(k)};
                smap_composite.night.time(k)                                             = dat.smap.observation_time{measurement_inds_night(k)};
                smap_composite.night.measurement_inds(k)                                 = measurement_inds_night(k);
                smap_composite.night.vegetation_opacity{k}(row_inds,col_inds)            = dat.smap.vegetation_opacity_sca{measurement_inds_night(k)};
                smap_composite.night.vegetation_opacity_dca{k}(row_inds,col_inds)        = dat.smap.vegetation_opacity{measurement_inds_night(k)};
                smap_composite.night.vegetation_water_content{k}(row_inds,col_inds)      = dat.smap.vegetation_water_content{measurement_inds_night(k)};
                smap_composite.night.surface_temperature{k}(row_inds,col_inds)           = dat.smap.surface_temperature{measurement_inds_night(k)};
                smap_composite.night.albedo{k}(row_inds,col_inds)                        = dat.smap.albedo_sca{measurement_inds_night(k)};
                smap_composite.night.albedo_dca{k}(row_inds,col_inds)                    = dat.smap.albedo{measurement_inds_night(k)};
                smap_composite.night.roughness{k}(row_inds,col_inds)                     = dat.smap.roughness_sca{measurement_inds_night(k)};
                smap_composite.night.roughness_dca{k}(row_inds,col_inds)                 = dat.smap.roughness{measurement_inds_night(k)};
                smap_composite.night.water_body_fraction{k}(row_inds,col_inds)           = dat.smap.water_body_fraction{measurement_inds_night(k)};
                smap_composite.night.surface_flags{k}(row_inds,col_inds)                 = dat.smap.surface_flags{measurement_inds_night(k)};
                rqfsca                                                                   = dat.smap.retrieval_quality_flags_sca{measurement_inds_night(k)};
                rqfdca                                                                   = dat.smap.retrieval_quality_flags{measurement_inds_night(k)};
                smap_composite.night.retrieval_quality_flags{k}(row_inds,col_inds)       = rqfsca;
                smap_composite.night.retrieval_quality_flags_dca{k}(row_inds,col_inds)   = rqfdca;
                x1                                                                       = dec2bin(rqfsca(:))-'0';
                x2                                                                       = dec2bin(rqfdca(:))-'0';
                smap_composite.night.retrieval_quality_flags_merge{k}(row_inds,col_inds) = reshape(bin2dec(num2str(x1|x2)),[num_pixels,num_pixels]);
            else
                rem_inds(k) = true;
            end
        end
    end
end

% remove invalid measurements
switch lower(dataset_label)
    case {'smapvex16m','vex16m','16m','manitoba','ma','m'}
        num_msmts   = length(smap_composite.night.soil_moisture);
        rem_inds    = false(num_msmts,1);
        rem_inds(3) = true;
        %
        smap_composite.night.soil_moisture(rem_inds)                 = [];
        smap_composite.night.soil_moisture_dca(rem_inds)             = [];
        smap_composite.night.time(rem_inds)                          = [];
        smap_composite.night.measurement_inds(rem_inds)              = [];
        smap_composite.night.vegetation_opacity(rem_inds)            = [];
        smap_composite.night.vegetation_opacity_dca(rem_inds)        = [];
        smap_composite.night.vegetation_water_content(rem_inds)      = [];
        smap_composite.night.surface_temperature(rem_inds)           = [];
        smap_composite.night.albedo(rem_inds)                        = [];
        smap_composite.night.albedo_dca(rem_inds)                    = [];
        smap_composite.night.roughness(rem_inds)                     = [];
        smap_composite.night.roughness_dca(rem_inds)                 = [];
        smap_composite.night.water_body_fraction(rem_inds)           = [];
        smap_composite.night.surface_flags(rem_inds)                 = [];
        smap_composite.night.retrieval_quality_flags(rem_inds)       = [];
        smap_composite.night.retrieval_quality_flags_dca(rem_inds)   = [];
        smap_composite.night.retrieval_quality_flags_merge(rem_inds) = [];
end
disp('Done.')
end
%--------------------------------------------------------------------------

%==========================================================================
function shrm = compute_shrm(sm,hr_mean,unobserved,pm,params)
%COMPUTE_SHRM Compute scaled high-resolution mean
mean_nan = @(a) mean(mean( a(~isnan(a)) ));
num_blocks    = params.LR_SIDE_PIXELS;
num_pixels    = params.HR_SIDE_PIXELS;
num_subpixels = params.SMAP_RES/params.SENTINEL_RES;

% compute baselines for blockwise data
if isempty(hr_mean)
    shrm = NaN(12,12);
else
    SHRM_cell = cell(num_blocks,num_blocks); % cell for each block's SHRM (climatology) predictions
    hr_mean(unobserved)  = NaN; % set missing values to NaN
    hr_mean(logical(pm)) = NaN;
    
    for i=1:num_blocks
        for j=1:num_blocks
            % initialize cell to hold SHRM predictions for block (i,j)
            SHRM_cell{i,j} = NaN(num_subpixels*num_subpixels,1);

            % get HR index vector to index pixels in block (i,j)
            hr_ind_vec = zeros(num_pixels,num_pixels);
            hr_ind_vec(((i-1)*num_subpixels+1):(i*num_subpixels),((j-1)*num_subpixels+1):(j*num_subpixels)) = 1;
            hr_ind_vec = logical(hr_ind_vec(:));

            % get LR index vector for block (i,j)
            lr_ind_vec      = zeros(num_blocks,num_blocks);
            lr_ind_vec(i,j) = 1;
            lr_ind_vec      = logical(lr_ind_vec(:));

            % get the HR mean image for this block
            hr_mean_block_ij = hr_mean(hr_ind_vec);

            % get HR mean of block (i,j)
            mu_hr_ij = mean_nan(hr_mean_block_ij);
            if isnan(mu_hr_ij)
                continue;
            end
            
            % get lr mean of block (i,j)
            lr_mean_block_ij = sm(lr_ind_vec);
            if isnan(lr_mean_block_ij)
                continue;
            end
            
            % compute SHRM for block (i,j)
            shrm_block_ij = (lr_mean_block_ij/mu_hr_ij)*hr_mean_block_ij;
            SHRM_cell{i,j} = shrm_block_ij;
        end
    end
    shrm = get_shrm_matrix(SHRM_cell,num_blocks,params.HR_SIDE_PIXELS);
end
end
%--------------------------------------------------------------------------

%==========================================================================
function shrm_matrix = get_shrm_matrix(SHRM,num_blocks,num_pixels)
%GET_SHRM_MATRIX
num_subpixels = num_pixels / num_blocks;
shrm_matrix   = NaN(num_pixels,num_pixels);
for i=1:num_blocks
    for j=1:num_blocks
        hr_ind_vec = zeros(num_pixels,num_pixels);
        hr_ind_vec(((i-1)*num_subpixels+1):(i*num_subpixels),((j-1)*num_subpixels+1):(j*num_subpixels)) = 1;
        hr_ind_vec = logical(hr_ind_vec(:));
        shrm_matrix(hr_ind_vec) = SHRM{i,j}(:);
    end
end
end
%--------------------------------------------------------------------------

%==========================================================================
function missing_pixels = get_missing_pixel_matrix_smap(sm,vo,alb,rough,temp,rqf,sf,vwc,wbf,p)
%GET_MISSING_PIXEL_MATRIX_SMAP Return matrix indicating which pixels in the
%soil moisture image are 'missing', where values in matrix are 1 if pixel
%is 'present' and 0 if pixel is 'missing'.
sm_miss    = ~isnan(sm);
vo_miss    = ~isnan(vo);
alb_miss   = ~isnan(alb);
rough_miss = ~isnan(rough);
temp_miss  = ~isnan(temp);
miss       = sm_miss & vo_miss & alb_miss & rough_miss & temp_miss;
%
sm_in_range = (sm>p.MIN_REASONABLE_SM)&(sm<p.MAX_REASONABLE_SM); % acceptable range of soil moisture
switch lower(p.MIN_PIXEL_QUALITY)
    case {'none'} % use all pixels in the HR data, regardless of quality
        missing_pixels = miss & sm_in_range;
    case {'recommended', 'rec', ''} % only use recommended pixels
        % get matrix indicating recommended quality pixels
        [num_rows,num_cols] = size(rqf);
        rqfbin = fliplr(dec2bin(rqf,16)); % get binary representation of flags
        rec_quality = reshape(rqfbin(:,1)=='0', num_rows, num_cols);
        missing_pixels = miss & sm_in_range & rec_quality;
    case {'uncertain', 'un'} % we suppose uncertain quality is OK if it 
        % meets our custom conditions
        % - we remove pixels that are anomalous, have retrieval issues, or
        %   correspond to certain surface condition issues where we expect
        %   the soil moisture measurement may be unreliable.
        
        % set NaNs in RQF and SF to have all bits flipped
        rqf(isnan(rqf)) = 2^16 - 1;
        sf(isnan(sf))   = 2^16 - 1;
        
        % get retrieval quality flags
        [num_rows,num_cols] = size(rqf); % get number of rows and cols
        rqfbin = fliplr(dec2bin(rqf,16)); % get binary representation of flags
        rec_quality                 = reshape(rqfbin(:,1)=='0', num_rows, num_cols); % recommended quality
        retrieval_quality_uncertain = reshape(rqfbin(:,1)=='1', num_rows, num_cols); % uncertain quality
        retrieval_skipped           = reshape(rqfbin(:,2)=='1', num_rows, num_cols); 
        retrieval_unsuccessful      = reshape(rqfbin(:,3)=='1', num_rows, num_cols); 
        
        % determine which pixels have acceptable retrieval quality
        no_retrieval_failure = ~(retrieval_unsuccessful | retrieval_skipped);
        acceptable_quality   = (rec_quality | retrieval_quality_uncertain) & no_retrieval_failure;
        
        % get surface conditions
        sfbin = fliplr(dec2bin(sf,16)); % bit code representation
        water                 = reshape((sfbin(:,1)=='1')|(sfbin(:,2)=='1'),num_rows,num_cols); % 1 if water fraction is greater than 0.05
        coastal_area          = reshape(sfbin(:,3)=='1', num_rows, num_cols);   % 1 if distance to significant water body is less than 1 36km grid cell
        urban_area            = reshape(sfbin(:,4)=='1', num_rows, num_cols);   % 1 if urban area fraction is greater than 0.25
        precipitation         = reshape(sfbin(:,5)=='1', num_rows, num_cols);   % 1 if precipitation is greater than 1 mm/hr
        snow_or_ice           = reshape(sfbin(:,6)=='1', num_rows, num_cols);   % 1 if snow area fraction is greater than 0.05
        permanent_snow_or_ice = reshape(sfbin(:,7)=='1', num_rows, num_cols);   % 1 if permanent ice fraction is greater than 0.05
        frozen_ground         = reshape((sfbin(:,8)=='1')|(sfbin(:,9)=='1'), num_rows, num_cols); % 1 if frozen ground fraction greater than 0.05
        mountainous_terrain   = reshape(sfbin(:,10)=='1', num_rows, num_cols);  % 1 if stdev of slope in cell is greater than 3 percent
        %high_vwc             = reshape(sfbin(:,11)=='1', num_rows, num_cols);  % 1 if VWC is greater than 5 kg/m^2
        
        % get pixels with VWC over threshold
        vwc_below_thresh = vwc <= p.MAX_ACCEPTABLE_VWC_SMAP;
        
        % get pixels with WBF over threshold
        wbf_below_thresh = wbf <= p.MAX_ACCEPTABLE_WBF;
        
        % determine which pixels are acceptable based on surface conditions
        bad_surface_condition = water | snow_or_ice | permanent_snow_or_ice;
        if p.REMOVE_COASTAL
            bad_surface_condition = bad_surface_condition | coastal_area;
        end
        if p.REMOVE_URBAN
            bad_surface_condition = bad_surface_condition | urban_area;
        end
        if p.REMOVE_FROZEN
            bad_surface_condition = bad_surface_condition | frozen_ground;
        end
        if p.REMOVE_MOUNTAINOUS
            bad_surface_condition = bad_surface_condition | mountainous_terrain;
        end
        if p.REMOVE_PRECIP
            bad_surface_condition = bad_surface_condition | precipitation;
        end
        acceptable_surface_condition = ~bad_surface_condition;
        
        % compute missing pixel matrix
        missing_pixels = miss & sm_in_range & vwc_below_thresh & wbf_below_thresh & acceptable_quality & acceptable_surface_condition;
    otherwise
        error('Invalid minimum pixel quality specified.');
end
end
%--------------------------------------------------------------------------

%==========================================================================
function plot_results_against_pals(datapath,dataset_label,pals,pred)
%PLOT_EOF_PREDICTIONS_AND_BASELINES_AGAINST_COMPOSITE_IMAGES
mean_nan = @(a) mean(mean( a(~isnan(a)) ));

% define error metrics
RMSD   = @(x,y,m) sqrt( (1/length(x(m)))*sum((x(m) - y(m)).^2) );                                                        % root mean square difference / error
MD     = @(x,y,m) (1/length(x(m)))*sum(x(m)-y(m));                                                                       % mean difference or bias
ubRMSD = @(x,y,m) sqrt( RMSD(x,y,m)^2 - MD(x,y,m)^2 );                                                                   % unbiased RMSD
R      = @(x,y,m) sum((x(m)-mean(x(m))).*(y(m)-mean(y(m)))) / sqrt(sum((x(m)-mean(x(m))).^2)*sum((y(m)-mean(y(m))).^2)); % Pearson correlation
MAD    = @(x,y,m) 1/length(x(m))*sum( abs(x(m) - y(m)) );                                                                % mean absolute difference

switch lower(dataset_label)
    case {'smapvex16i','16i','vex16i','iowa','ia','i'}
        smap_measurement_dates_day = {'26-May-2016','28-May-2016',...      % dates are in UTC
            '31-May-2016','02-Jun-2016','03-Jun-2016','05-Jun-2016',...
            '03-Aug-2016','05-Aug-2016','06-Aug-2016','08-Aug-2016',...
            '11-Aug-2016','13-Aug-2016','14-Aug-2016','16-Aug-2016',...
            '19-Aug-2016'};
        measurement_dates_night = {'27-May-2016','28-May-2016',...
            '31-May-2016','02-Jun-2016','04-Jun-2016','05-Jun-2016',...
            '03-Aug-2016','05-Aug-2016','07-Aug-2016','08-Aug-2016',...
            '11-Aug-2016','13-Aug-2016','15-Aug-2016','16-Aug-2016',...
            '19-Aug-2016'};
        pals_measurement_dates  = {'08-June-2016','11-June-2016',...
            '14-June-2016','16-June-2016','19-June-2016','20-June-2016',...
            '14-July-2016','16-July-2016','18-July-2016','19-July-2016',...
            '21-July-2016','22-July-2016'};
        time_offset  = 5/24;
        m_pals       = cellfun(@(x) datestr(mean_nan(x.km1.time)-time_offset),pals,'UniformOutput',false);
        row_range    = 10:45;
        col_range    = 16:51;
        cmin         = 0.02;
        cmax         = 0.5726;
        inds_to_plot = [3,4,12];
    case {'smapvex16m','16m','vex16m','manitoba','ma','m'}
        measurement_dates_day = {'06-Jun-2016','08-Jun-2016',...
            '09-Jun-2016','11-Jun-2016','12-Jun-2016','14-Jun-2016',...
            '16-Jun-2016','19-Jun-2016','20-Jun-2016','22-Jun-2016',...
            '13-Jul-2016','14-Jul-2016','16-Jul-2016','18-Jul-2016',...
            '19-Jul-2016','21-Jul-2016','22-Jul-2016','24-Jul-2016'};
        measurement_dates_night = {'07-Jun-2016','08-Jun-2016',...
            '10-Jun-2016','12-Jun-2016','13-Jun-2016','15-Jun-2016',...
            '16-Jun-2016','18-Jun-2016','20-Jun-2016','21-Jun-2016',...
            '23-Jun-2016','12-Jul-2016','14-Jul-2016','15-Jul-2016',...
            '17-Jul-2016','18-Jul-2016','20-Jul-2016','22-Jul-2016',...
            '23-Jul-2016',};
        time_offset = 5/24;
        pals_times = zeros(length(pals),1);
        for k=1:length(pals_times)
            pals_times(k) = pals{k}.km1.time(find(not(isnan(pals{k}.km1.time)),1)) - time_offset;
        end
        m_pals       = mat2cell(pals_times',1,ones(1,length(pals_times')));
        row_range    = 25:60;
        col_range    = 10:45;
        cmin         = 0.02;
        cmax         = 0.7624;
        inds_to_plot = [6,12,13];
    case {'smapvex15a','15a','vex15a','arizona','az','a'}
        % get measurement times
        measurement_dates_day = {'31-Jul-2015','02-Aug-2015',...
            '03-Aug-2015','05-Aug-2015','08-Aug-2015','10-Aug-2015',...
            '11-Aug-2015','13-Aug-2015','18-Aug-2015','19-Aug-2015',...
            '21-Aug-2015'};
        measurement_dates_night = {'01-Aug-2015','02-Aug-2015',...
            '04-Aug-2015','07-Aug-2015','09-Aug-2015','10-Aug-2015',...
            '12-Aug-2015','15-Aug-2015','17-Aug-2015','18-Aug-2015',...
            '20-Aug-2015'};
        time_offset  = 7/24;
        cmin         = 0.02;
        cmax         = 0.5121;
        row_range    = 26:68;
        col_range    = 29:143;
        times_pals   = cellfun(@(x) datestr(mean_nan(x.km1.time)-time_offset),pals,'UniformOutput',false);
        m_pals       = cellfun(@(x) datenum(x),times_pals);
        inds_to_plot = [2,5,6];
    otherwise
        error('invalid dataset label');
end

% read ground truth data into pred
ff = @(x) cellfun(@(y) floor(datenum(y)),x);     % convert cell array to vector and truncate hour/min/seconds (just want to find measurement on same day)
[num_row_pixels,num_col_pixels] = size(pals{1}.km1.vsm);
num_validation_set_measurements = length(pred.day.lr);
num_figures                     = ceil(num_validation_set_measurements / 4);
%
pred.day.hr                     = cell(num_validation_set_measurements,1);
pred.day.hr_times               = cell(num_validation_set_measurements,1);
for k=1:num_validation_set_measurements
    switch lower(dataset_label)
        case {'smapvex16i','16i','vex16i','iowa','ia','i','smapvex16m','16m','vex16m','manitoba','ma','m'}
            ind   = find(floor(pred.day.measurement_times(k))==ff(m_pals));
        otherwise
            ind = find(floor(pred.day.measurement_times(k))==floor(m_pals));
    end
    if isempty(ind)
        pred.day.hr{k}       = NaN(num_row_pixels,num_col_pixels);
        pred.day.hr_times{k} = '';
    else
        pred.day.hr{k}       = pals{ind}.km1.vsm;
        pred.day.hr_times{k} = m_pals(ind);
    end
end

% compute errors against HR (PALS) data
day_or_night          = 'day';
errors.(day_or_night) = cell(num_validation_set_measurements,1);
min_LR                = 1;
min_RAW_SHRM          = 1;
min_MC_SHRM           = 1;
min_SR                = 1;
min_HR                = 1;
max_LR                = 0;
max_RAW_SHRM          = 0;
max_MC_SHRM           = 0;
max_SR                = 0;
max_HR                = 0;
for k=1:num_validation_set_measurements
    % read out data
    LR_k       = pred.day.lr{k};
    RAW_SHRM_k = pred.day.raw_shrm{k};
    MC_SHRM_k  = pred.day.mc_shrm{k};
    SR_k       = pred.day.sr{k};
    HR_k       = pred.day.hr{k};

    % get missing values in common
    miss_lr       = not(isnan(LR_k));
    miss_raw_shrm = not(isnan(RAW_SHRM_k));
    miss_mc_shrm  = not(isnan(MC_SHRM_k));
    miss_sr       = not(isnan(SR_k));
    miss_hr       = not(isnan(HR_k));
    M_k           = miss_lr & miss_raw_shrm & miss_mc_shrm & miss_sr & miss_hr;
    
    % remove non-common values
    LR_k(not(M_k))          = NaN;
    RAW_SHRM_k(not(M_k))    = NaN;
    MC_SHRM_k(not(M_k))     = NaN;
    SR_k(not(M_k))          = NaN;
    HR_k(not(M_k))          = NaN;

    % get min and max values
    min_LR       = min([min_LR,       min(min(LR_k))]);
    min_RAW_SHRM = min([min_RAW_SHRM, min(min(RAW_SHRM_k))]);
    min_MC_SHRM  = min([min_MC_SHRM,  min(min(MC_SHRM_k))]);
    min_SR       = min([min_SR,       min(min(SR_k))]);
    min_HR       = min([min_HR,       min(min(HR_k))]);
    %
    max_LR       = max([max_LR,       max(max(LR_k))]);
    max_RAW_SHRM = max([max_RAW_SHRM, max(max(RAW_SHRM_k))]);
    max_MC_SHRM  = max([max_MC_SHRM,  max(max(MC_SHRM_k))]);
    max_SR       = max([max_SR,       max(max(SR_k))]);
    max_HR       = max([max_HR,       max(max(HR_k))]);

    % compute errors
    errors.(day_or_night){k}.PALS.LR.RMSD             = RMSD(LR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.LR.MD               = MD(LR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.LR.ubRMSD           = ubRMSD(LR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.LR.MAD              = MAD(LR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.LR.R                = R(LR_k,HR_k,M_k);
    %
    errors.(day_or_night){k}.PALS.RAW_SHRM.RMSD       = RMSD(RAW_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.RAW_SHRM.MD         = MD(RAW_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.RAW_SHRM.ubRMSD     = ubRMSD(RAW_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.RAW_SHRM.MAD        = MAD(RAW_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.RAW_SHRM.R          = R(RAW_SHRM_k,HR_k,M_k);
    %
    errors.(day_or_night){k}.PALS.MC_SHRM.RMSD        = RMSD(MC_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.MC_SHRM.MD          = MD(MC_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.MC_SHRM.ubRMSD      = ubRMSD(MC_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.MC_SHRM.MAD         = MAD(MC_SHRM_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.MC_SHRM.R           = R(MC_SHRM_k,HR_k,M_k);
    %
    errors.(day_or_night){k}.PALS.SR.RMSD             = RMSD(SR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.SR.MD               = MD(SR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.SR.ubRMSD           = ubRMSD(SR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.SR.MAD              = MAD(SR_k,HR_k,M_k);
    errors.(day_or_night){k}.PALS.SR.R                = R(SR_k,HR_k,M_k);
    %
    errors.(day_or_night){k}.PALS.LR.num_pixels       = sum(sum(M_k));
    errors.(day_or_night){k}.PALS.RAW_SHRM.num_pixels = sum(sum(M_k));
    errors.(day_or_night){k}.PALS.MC_SHRM.num_pixels  = sum(sum(M_k));
    errors.(day_or_night){k}.PALS.SR.num_pixels       = sum(sum(M_k));
end

% get min and max of colorbar range
cmin = min([min_LR,min_RAW_SHRM,min_MC_SHRM,min_SR,min_HR]);
cmax = max([max_LR,max_RAW_SHRM,max_MC_SHRM,max_SR,max_HR]);

% plot composite test set images with error values displayed
load([datapath,'/SMAP_Color_SoilMoisture.mat']);
if inds_to_plot==-1
    num_figs = ceil(num_validation_set_measurements/5);
else
    num_figs = ceil(length(inds_to_plot)/5);
end
for n=1:num_figs
    figure
    for k=1:5
        if inds_to_plot==-1
            index = (n-1)*5 + k;
        else
            if k>length(inds_to_plot)
                break;
            else
                index = inds_to_plot(k);
            end
        end
        if index>num_validation_set_measurements
            break;
        end
        pals_missing = isnan(pred.day.hr{index});
        
        % LR
        subplot_ind = k;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        lr_im = pred.day.lr{index};
        lr_im(pals_missing) = NaN;
        imagesc(lr_im(row_range,col_range)), caxis manual, caxis([cmin,cmax])
        colormap(ax(subplot_ind),uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title({'SMAP SCA',datestr(pred.day.measurement_times(index)-time_offset)},'interpreter','latex')
        
        % RAW SHRM
        subplot_ind = k+5;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        shrm_im = pred.day.raw_shrm{index};
        shrm_im(pals_missing) = NaN;
        imagesc(shrm_im(row_range,col_range)), caxis manual, caxis([cmin,cmax])
        colormap(ax(subplot_ind),uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('Raw Data Climatology','interpreter','latex')
        
        % MC SHRM
        subplot_ind = k+10;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        shrm_im = pred.day.mc_shrm{index};
        shrm_im(pals_missing) = NaN;
        imagesc(shrm_im(row_range,col_range)), caxis manual, caxis([cmin,cmax])
        colormap(ax(subplot_ind),uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('MC Climatology','interpreter','latex')
        
        % SR
        subplot_ind = k+15;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        sr_im = pred.day.sr{index};
        sr_im(pals_missing) = NaN;
        imagesc(sr_im(row_range,col_range)), caxis manual, caxis([cmin,cmax])
        colormap(ax(subplot_ind),uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('EOF Superresolution','interpreter','latex')
        
        % HR
        subplot_ind = k+20;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        imagesc(pred.day.hr{index}(row_range,col_range)), caxis manual, caxis([cmin,cmax])
        colormap(ax(subplot_ind),uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('PALS','interpreter','latex')
        
        % display errors
        % LR
        r = 0.164 * mod(k-1,5);
        % 0.706, 0.692
        annotation('textbox',  [0.155+r, 0.707,  0.1, 0.1], 'String', "RMSD   = " + errors.(day_or_night){index}.PALS.LR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.65,  0.1, 0.1], 'String', "MD     = " + errors.(day_or_night){index}.PALS.LR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.696, 0.1, 0.1], 'String', "ubRMSD = " + errors.(day_or_night){index}.PALS.LR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.684,  0.1, 0.1], 'String', "R      = " + errors.(day_or_night){index}.PALS.LR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

        % RAW SHRM
        % 0.533, 0.519
        annotation('textbox',  [0.155+r, 0.534,  0.1, 0.1], 'String', "RMSD   = " + errors.(day_or_night){index}.PALS.RAW_SHRM.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.68,  0.1, 0.1], 'String', "MD     = " + errors.(day_or_night){index}.PALS.RAW_SHRM.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.522, 0.1, 0.1], 'String', "ubRMSD = " + errors.(day_or_night){index}.PALS.RAW_SHRM.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.510,  0.1, 0.1], 'String', "R      = " + errors.(day_or_night){index}.PALS.RAW_SHRM.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        
        % MC SHRM
        % 0.36, 0.346
        annotation('textbox',  [0.155+r, 0.362,  0.1, 0.1], 'String', "RMSD   = " + errors.(day_or_night){index}.PALS.MC_SHRM.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.68,  0.1, 0.1], 'String', "MD     = " + errors.(day_or_night){index}.PALS.MC_SHRM.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.350, 0.1, 0.1], 'String', "ubRMSD = " + errors.(day_or_night){index}.PALS.MC_SHRM.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.338,  0.1, 0.1], 'String', "R      = " + errors.(day_or_night){index}.PALS.MC_SHRM.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        
        % EOF SR
        % 0.187, 0.173
        annotation('textbox',  [0.155+r, 0.19,  0.1, 0.1], 'String',  "RMSD   = " + errors.(day_or_night){index}.PALS.SR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.68,  0.1, 0.1], 'String',  "MD     = " + errors.(day_or_night){index}.PALS.SR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.178,  0.1, 0.1], 'String', "ubRMSD = " + errors.(day_or_night){index}.PALS.SR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.166,  0.1, 0.1], 'String',  "R      = " + errors.(day_or_night){index}.PALS.SR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        % 0.188, 0.174, 0.160
    end
    
    % insert common colorbar
    h=colorbar; h.TickLabelInterpreter='latex'; colormap(uint8(SMAP_Color_Scale));
    caxis([cmin,cmax])
    set(h, 'Position', [0.9, 0.15, 0.02, 0.7])
    clear ax
end
end
%--------------------------------------------------------------------------

%==========================================================================
function filename_set = get_filename_set(dataset_label)
%GET_FILENAME_SET Gets the filenames for all files at a given SMAPVEX site,
%where each file represents the data for a specific 36km EASE grid cell.
switch lower(dataset_label)
    case {'smapvex15a','smapvex15','15','vex15','arizona','az'}
        ease_rows = 96:98;
        ease_cols = 185:188;
        %
        num_rows  = length(ease_rows);
        num_cols  = length(ease_cols);
        %
        R_str = 'R0';
        C_str = 'C';
    case {'smapvex16i','16i','i','iowa','vex16i','ia'}
        ease_rows = 66:67;
        ease_cols = 232:233;
        %
        num_rows  = length(ease_rows);
        num_cols  = length(ease_cols);
        %
        R_str = 'R0';
        C_str = 'C';
    case {'smapvex16m','16m','m','manitoba','vex16m','ma'}
        ease_rows = 48:49;
        ease_cols = 220:221;
        %
        num_rows  = length(ease_rows);
        num_cols  = length(ease_cols);
        %
        R_str     = 'R0';
        C_str     = 'C';
    otherwise
        error('invalid dataset label');
end

filename_set = cell(num_rows,num_cols);
for i=1:length(ease_rows)
    for j=1:length(ease_cols)
        filename_set{i,j} = [R_str,num2str(ease_rows(i)),C_str,num2str(ease_cols(j))];
    end
end
end
%--------------------------------------------------------------------------

%==========================================================================
function [raw_data_hr_means,mc_hr_means,eof,pm] = load_hr_means_and_eof_data(datapath,dataset_label,params)
%LOAD_HR_MEANS_AND_EOF_DATA
disp(['Loading EOF data and computing HR means for ',dataset_label])
filename_set        = get_filename_set(dataset_label);
[num_rows,num_cols] = size(filename_set);
num_subpixels       = params.SMAP_RES/params.SENTINEL_RES;
%
raw_data_hr_means   = cell(num_rows,num_cols);
mc_hr_means         = cell(num_rows,num_cols);
eof                 = cell(num_rows,num_cols);
pm                  = cell(num_rows,num_cols);
for i=1:num_rows
    for j=1:num_cols
        disp(['    Reading data for 36km EASE cell (',num2str(i),',',num2str(j),') of a ',num2str(num_rows), 'x', num2str(num_cols), ' array...'])
        
        % load the raw data for this EASE cell
        filepath_eof    = [datapath,'/EOF/',filename_set{i,j},'.mat'];
        eof_file_exists = (exist(filepath_eof,'file')==2);
        if eof_file_exists
            dat = load(filepath_eof);
            EOF = dat.EOF;
        else
            continue;
        end
        
        % prepare training data
        persistently_missing = EOF.persistently_missing;
        
        % get persistently missing pixels
        pm{i,j}.day   = persistently_missing.day;
        pm{i,j}.night = persistently_missing.night;
        
        % get raw data HR means
        raw_data_hr_means{i,j}.day.sm           = EOF.hr_mean.day;
        raw_data_hr_means{i,j}.day.unobserved   = isnan(raw_data_hr_means{i,j}.day.sm);
        %
        raw_data_hr_means{i,j}.night.sm         = EOF.hr_mean.night;
        raw_data_hr_means{i,j}.night.unobserved = isnan(raw_data_hr_means{i,j}.night.sm);
        
        % get MC HR means
        mv_array = EOF.params.day.mean_vector;
        for k=1:numel(mv_array)
            if isempty(mv_array{k})
                mv_array{k} = NaN(num_subpixels);
            end
        end
        mc_hr_means{i,j}.day.sm           = cell2mat(cellfun(@(x) reshape(x,[num_subpixels,num_subpixels]),mv_array,'UniformOutput',false));
        mc_hr_means{i,j}.day.unobserved   = isnan(mc_hr_means{i,j}.day.sm);
        %
        mv_array = EOF.params.night.mean_vector;
        for k=1:numel(mv_array)
            if isempty(mv_array{k})
                mv_array{k} = NaN(num_subpixels);
            end
        end
        mc_hr_means{i,j}.night.sm         = cell2mat(cellfun(@(x) reshape(x,[num_subpixels,num_subpixels]),mv_array,'UniformOutput',false));
        mc_hr_means{i,j}.night.unobserved = isnan(mc_hr_means{i,j}.night.sm);
        
        % get EOF data
        if isfield(EOF.params,'day')
            eof{i,j}.day.E                 = EOF.params.day.E;
            eof{i,j}.day.C                 = EOF.params.day.C;
            eof{i,j}.day.Lambda            = EOF.params.day.Lambda;
            eof{i,j}.day.Beta              = EOF.params.day.Beta;
            eof{i,j}.day.mean_vector       = EOF.params.day.mean_vector;
            eof{i,j}.day.test              = EOF.test.day;
            eof{i,j}.day.valid             = EOF.valid.day;
        end
        if isfield(EOF.params,'night')
            eof{i,j}.night.E               = EOF.params.night.E;
            eof{i,j}.night.C               = EOF.params.night.C;
            eof{i,j}.night.Lambda          = EOF.params.night.Lambda;
            eof{i,j}.night.Beta            = EOF.params.night.Beta;
            eof{i,j}.night.mean_vector     = EOF.params.night.mean_vector;
            eof{i,j}.night.test            = EOF.test.night;
            eof{i,j}.night.valid           = EOF.valid.night;
        end
    end
end
disp('    Done.')
end
%--------------------------------------------------------------------------

%==========================================================================
function output_struct = get_EOF_predictions_and_baselines(dataset_label,smap,raw_data_hr_means,mc_hr_means,eof,pm,num_components,R2_threshold,params)
%GET_EOF_PREDICTIONS_AND_BASELINES Compute the EOF-based high-resolution
%soil moisture predictions from each LR measurement, and also prepare some
%baseline methods to compare the predictions against.
% OUTPUT:
%   output_struct = struct containing the following fields:
%       output_struct.lr   = Nx1 cell array containing full LR images
%       output_struct.shrm = Nx1 cell array containing full SHRM images
%       output_struct.sr   = Nx1 cell array containing full SR images
mean_nan   = @(a) mean(mean( a(~isnan(a)) ));
num_pixels = params.HR_SIDE_PIXELS;
switch lower(dataset_label)
    case {'smapvex16i','vex16i','16i','iowa','ia'}
        ease_rows = 66:67;
        ease_cols = 232:233;
        measurement_dates_day = {'26-May-2016','28-May-2016',...           
            '31-May-2016','02-Jun-2016','03-Jun-2016','05-Jun-2016',...
            '03-Aug-2016','05-Aug-2016','06-Aug-2016','08-Aug-2016',...
            '11-Aug-2016','13-Aug-2016','14-Aug-2016','16-Aug-2016',...
            '19-Aug-2016'};
        measurement_dates_night = {'27-May-2016','28-May-2016',...
            '31-May-2016','02-Jun-2016','04-Jun-2016','05-Jun-2016',...
            '03-Aug-2016','05-Aug-2016','07-Aug-2016','08-Aug-2016',...
            '11-Aug-2016','13-Aug-2016','15-Aug-2016','16-Aug-2016',...
            '19-Aug-2016'};
    case {'smapvex16m','vex16m','16m','manitoba','ma'}
        ease_rows = 48:49;
        ease_cols = 220:221;
        measurement_dates_day = {'06-Jun-2016','08-Jun-2016',...
            '09-Jun-2016','11-Jun-2016','12-Jun-2016','14-Jun-2016',...
            '16-Jun-2016','19-Jun-2016','20-Jun-2016','22-Jun-2016',...
            '13-Jul-2016','14-Jul-2016','16-Jul-2016','18-Jul-2016',...
            '19-Jul-2016','21-Jul-2016','22-Jul-2016','24-Jul-2016'};
        measurement_dates_night = {'07-Jun-2016','08-Jun-2016',...
            '10-Jun-2016','12-Jun-2016','13-Jun-2016','15-Jun-2016',...
            '16-Jun-2016','18-Jun-2016','20-Jun-2016','21-Jun-2016',...
            '23-Jun-2016','12-Jul-2016','14-Jul-2016','15-Jul-2016',...
            '17-Jul-2016','18-Jul-2016','20-Jul-2016','22-Jul-2016',...
            '23-Jul-2016'};
    case {'smapvex15a','smapvex15','vex15','15','arizona','az'}
        ease_rows = 96:98;
        ease_cols = 185:188;
        measurement_dates_day = {'31-Jul-2015','02-Aug-2015',...
            '03-Aug-2015','05-Aug-2015','08-Aug-2015','10-Aug-2015',...
            '11-Aug-2015','13-Aug-2015','18-Aug-2015','19-Aug-2015',...
            '21-Aug-2015'};
        measurement_dates_night = {'01-Aug-2015','02-Aug-2015',...
            '04-Aug-2015','07-Aug-2015','09-Aug-2015','10-Aug-2015',...
            '12-Aug-2015','15-Aug-2015','17-Aug-2015','18-Aug-2015',...
            '20-Aug-2015'};
    otherwise
        error('invalid dataset label');
end
num_measurements_day   = length(measurement_dates_day);
num_measurements_night = length(measurement_dates_night);
num_rows               = length(ease_rows);
num_cols               = length(ease_cols);

% get average measurement time of DAY measurements over the SMAPVEX area
measurement_times_day = NaN(num_measurements_day,1);
assert(length(smap{1,1}.day.time)==num_measurements_day)
for k=1:num_measurements_day
    all_msmts_in_array = NaN(num_rows,num_cols);
    for i=1:num_rows
        for j=1:num_cols
            all_msmts_in_array(i,j) = smap{i,j}.day.time(k);
        end
    end
    measurement_times_day(k) = mean_nan(all_msmts_in_array);
end

% get average measurement time of NIGHT measurements over the SMAPVEX area
measurement_times_night = NaN(num_measurements_night,1);
assert(length(smap{1,1}.night.time)==num_measurements_night)
for k=1:num_measurements_night
    all_msmts_in_array = NaN(num_rows,num_cols);
    for i=1:num_rows
        for j=1:num_cols
            all_msmts_in_array(i,j) = smap{i,j}.night.time(k);
        end
    end
    measurement_times_night(k) = mean_nan(all_msmts_in_array);
end

% get full LR, RD-SHRM, MC-SHRM, and SR images for DAY data
% - get image for each cell and insert into appropriate spot in full image
upscale = @(x) kron(x,ones(9/params.SENTINEL_RES));
vec     = @(x) x(:);
%
lr       = cell(num_measurements_day,1);
vod      = cell(num_measurements_day,1);
vwc      = cell(num_measurements_day,1);
raw_shrm = cell(num_measurements_day,1);
mc_shrm  = cell(num_measurements_day,1);
sr       = cell(num_measurements_day,1);
for k=1:num_measurements_day
    % initialize arrays
    lr{k}       = NaN(num_rows*num_pixels,num_cols*num_pixels);
    vod{k}      = NaN(num_rows*num_pixels,num_cols*num_pixels);
    vwc{k}      = NaN(num_rows*num_pixels,num_cols*num_pixels);
    raw_shrm{k} = NaN(num_rows*num_pixels,num_cols*num_pixels);
    mc_shrm{k}  = NaN(num_rows*num_pixels,num_cols*num_pixels);
    sr{k}       = NaN(num_rows*num_pixels,num_cols*num_pixels);
    for i=1:num_rows
        for j=1:num_cols
            % get indexing information for this cell
            indexing_vec = false(num_rows*num_pixels,num_cols*num_pixels);
            indexing_vec( (num_pixels*(i-1)+1):(num_pixels*i), (num_pixels*(j-1)+1):(num_pixels*j) ) = true;
            
            % read out data
            sm             = smap{i,j}.day.soil_moisture{k};
            vo             = smap{i,j}.day.vegetation_opacity{k};
            albedo         = smap{i,j}.day.albedo{k};
            roughness      = smap{i,j}.day.roughness{k};
            temp           = smap{i,j}.day.surface_temperature{k};
            rqf            = smap{i,j}.day.retrieval_quality_flags{k};
            sf             = smap{i,j}.day.surface_flags{k};
            wbf            = smap{i,j}.day.water_body_fraction{k};
            vwc_k          = smap{i,j}.day.vegetation_water_content{k};
            
            % get missing pixels
            missing_pixels_smap_sca = get_missing_pixel_matrix_smap(sm,vo,albedo,roughness,temp,rqf,sf,vwc_k,wbf,params);
            if (sum(sum(missing_pixels_smap_sca))==0)
                continue;
            end
            
            % remove bad pixels from SCA data
            sm(not(missing_pixels_smap_sca))   = NaN;
            vo(not(missing_pixels_smap_sca))   = NaN;
            temp(not(missing_pixels_smap_sca)) = NaN;
            
            % insert the LR (SCA) image for this cell into the full LR image
            lr{k}(indexing_vec) = vec(upscale(sm));
            if sum(sum(isnan(lr{k}(indexing_vec))))==numel(lr{k}(indexing_vec))
                continue;
            end
            
            % insert the Raw-Data-SHRM for this cell
            if ~isempty(raw_data_hr_means{i,j})
                if ~isempty(raw_data_hr_means{i,j}.day)
                    shrm_cell                 = compute_shrm(sm,raw_data_hr_means{i,j}.day.sm,raw_data_hr_means{i,j}.day.unobserved,pm{i,j}.day,params);
                    raw_shrm{k}(indexing_vec) = shrm_cell(:);
                end
            end
            
            % insert the MC-SHRM for this cell
            if ~isempty(mc_hr_means{i,j})
                if ~isempty(mc_hr_means{i,j}.day)
                    shrm_cell                = compute_shrm(sm,mc_hr_means{i,j}.day.sm,mc_hr_means{i,j}.day.unobserved,pm{i,j}.day,params);
                    mc_shrm{k}(indexing_vec) = shrm_cell(:);
                end
            end
            
            % get the EOF superresolution prediction for this cell
            if isfield(eof{i,j},'day')
                sr{k}(indexing_vec)  = superresolve_LR_image_blockwise(sm,vo,temp,eof{i,j}.day,num_components,R2_threshold,params);
                vod{k}(indexing_vec) = upscale(vo);
                vwc{k}(indexing_vec) = upscale(vwc_k);
            end
        end
    end
end
output_struct.day.lr                = lr;
output_struct.day.vod               = vod;
output_struct.day.vwc               = vwc;
output_struct.day.raw_shrm          = raw_shrm;
output_struct.day.mc_shrm           = mc_shrm;
output_struct.day.sr                = sr;
output_struct.day.measurement_times = measurement_times_day;

% get full LR, SHRM, and SR images for NIGHT data
lr       = cell(num_measurements_night,1);
vod      = cell(num_measurements_night,1);
vwc      = cell(num_measurements_night,1);
raw_shrm = cell(num_measurements_night,1); 
mc_shrm  = cell(num_measurements_night,1);
sr       = cell(num_measurements_night,1);
for k=1:num_measurements_night
    % initialize arrays
    lr{k}       = NaN(num_rows*num_pixels,num_cols*num_pixels);
    vod{k}      = NaN(num_rows*num_pixels,num_cols*num_pixels);
    vwc{k}      = NaN(num_rows*num_pixels,num_cols*num_pixels);
    raw_shrm{k} = NaN(num_rows*num_pixels,num_cols*num_pixels);
    mc_shrm{k}  = NaN(num_rows*num_pixels,num_cols*num_pixels);
    sr{k}       = NaN(num_rows*num_pixels,num_cols*num_pixels);
    for i=1:num_rows
        for j=1:num_cols
            % get indexing information for this cell
            indexing_vec = false(num_rows*num_pixels,num_cols*num_pixels);
            indexing_vec( (num_pixels*(i-1)+1):(num_pixels*i), (num_pixels*(j-1)+1):(num_pixels*j) ) = true;
            
            % read out data
            sm             = smap{i,j}.night.soil_moisture{k};
            vo             = smap{i,j}.night.vegetation_opacity{k};
            albedo         = smap{i,j}.night.albedo{k};
            roughness      = smap{i,j}.night.roughness{k};
            temp           = smap{i,j}.night.surface_temperature{k};
            rqf            = smap{i,j}.night.retrieval_quality_flags{k};
            sf             = smap{i,j}.night.surface_flags{k};
            wbf            = smap{i,j}.night.water_body_fraction{k};
            vwc_k          = smap{i,j}.night.vegetation_water_content{k};
            
            % get missing pixels
            missing_pixels_smap_sca = get_missing_pixel_matrix_smap(sm,vo,albedo,roughness,temp,rqf,sf,vwc_k,wbf,params);
            if (sum(sum(missing_pixels_smap_sca))==0)
                continue;
            end
            
            % remove bad pixels from SCA data
            sm(not(missing_pixels_smap_sca))            = NaN;
            vo(not(missing_pixels_smap_sca))            = NaN;
            temp(not(missing_pixels_smap_sca))          = NaN;
            
            % insert the LR (SCA) image for this cell into the full LR image
            lr{k}(indexing_vec) = vec(upscale(sm));
            if sum(sum(isnan(lr{k}(indexing_vec))))==numel(lr{k}(indexing_vec))
                continue;
            end
            
            % insert the Raw-Data-SHRM for this cell
            if ~isempty(raw_data_hr_means{i,j})
                if ~isempty(raw_data_hr_means{i,j}.night)
                    shrm_cell                 = compute_shrm(sm,raw_data_hr_means{i,j}.night.sm,raw_data_hr_means{i,j}.night.unobserved,pm{i,j}.night,params);
                    raw_shrm{k}(indexing_vec) = shrm_cell(:);
                end
            end
            
            % insert the MC-SHRM for this cell
            if ~isempty(mc_hr_means{i,j})
                if ~isempty(mc_hr_means{i,j}.night)
                    shrm_cell                = compute_shrm(sm,mc_hr_means{i,j}.night.sm,mc_hr_means{i,j}.night.unobserved,pm{i,j}.night,params);
                    mc_shrm{k}(indexing_vec) = shrm_cell(:);
                end
            end
            
            % get the EOF superresolution prediction for this cell
            if isfield(eof{i,j},'night')
                sr{k}(indexing_vec)  = superresolve_LR_image_blockwise(sm,vo,temp,eof{i,j}.night,num_components,R2_threshold,params);
                vod{k}(indexing_vec) = upscale(vo);
                vwc{k}(indexing_vec) = upscale(vwc_k);
            end
        end
    end
end
output_struct.night.lr                = lr;
output_struct.night.vod               = vod;
output_struct.night.vwc               = vwc;
output_struct.night.raw_shrm          = raw_shrm;
output_struct.night.mc_shrm           = mc_shrm;
output_struct.night.sr                = sr;
output_struct.night.measurement_times = measurement_times_night;
end
%--------------------------------------------------------------------------

%==========================================================================
function valid_data = plot_results_against_smap_sentinel(EOF,raw_data_hr_means,mc_hr_means,pm,num_components,dataset_label,R2_threshold,datapath,params)
%COMPUTE_AND_PLOT_SMAP_SENTINEL_RESULTS_AND_ERRORS
switch lower(dataset_label)
    case {'smapvex16i','vex16i','16i','iowa','ia'}
        time_offset        = 5/24;
        cmin               = 0.02;
        cmax               = 0.5726;
        inds_to_plot.day   = [];      % -1 will plot all inds
        inds_to_plot.night = [4,2,5]; % -1 will plot all inds
    case {'smapvex16m','vex16m','16m','manitoba','ma'}
        time_offset        = 5/24;
        cmin               = 0.02;
        cmax               = 0.7624;
        inds_to_plot.day   = [7];
        inds_to_plot.night = [7,4];
    case {'smapvex15a','smapvex15','vex15','15','arizona','az'}
        time_offset        = 7/24;
        cmin               = 0.02;
        cmax               = 0.5121;
        inds_to_plot.day   = [2,3,5];
        inds_to_plot.night = [];
    otherwise
        error('invalid dataset label');
end

% Compute errors and plot validation set results
labels = {'day','night'};
for k=1:length(labels)
    day_or_night = labels{k};
    
    % Form composite test set images
    [valid_data.(day_or_night),measurement_times] = get_composite_test_set_images(EOF,'valid',day_or_night,...
        raw_data_hr_means,mc_hr_means,pm,num_components,R2_threshold,params);
    
    % Truncate unrealistic soil moisture values
    valid_data.(day_or_night) = truncate_unrealistic_values(valid_data.(day_or_night),params);
    
    % Compute errors on full set of composite test set images
    calculate_errors_and_plot_results(valid_data.(day_or_night),day_or_night,measurement_times,time_offset,...
        inds_to_plot.(day_or_night),cmin,cmax,datapath);
end
end
%--------------------------------------------------------------------------

%==========================================================================
function [data,measurement_times] = get_composite_test_set_images(EOF,test_or_valid,day_or_night,raw_data_hr_means,mc_hr_means,pm,num_components,R2_threshold,params)
%GET_COMPOSITE_TEST_SET_IMAGES
num_pixels = params.HR_SIDE_PIXELS;
[num_rows,num_cols] = size(EOF);
%
vec      = @(x) x(:);
upscale  = @(x) vec(kron(reshape(x,[4,4]),ones(9/params.SENTINEL_RES)));
%
data              = {};
measurement_times = [];
cell_index        = 1;
for i=1:num_rows
    for j=1:num_cols
        % read out the test set for cell (i,j)
        test             = EOF{i,j}.(day_or_night).(test_or_valid);
        num_measurements = size(test.Y,2);
        
        % get indexing vector for cell (i,j)
        index_vector = false(num_rows*num_pixels,num_cols*num_pixels);
        index_vector( ((i-1)*num_pixels+1):(i*num_pixels), ((j-1)*num_pixels+1):(j*num_pixels) ) = true;
        
        % extract data for cell (i,j)
        for k=1:num_measurements
            % read out data for test set sample k
            yk       = test.Y(:,k);
            mk       = test.M(:,k);
            lk       = test.L(:,k);
            vk       = test.V(:,k);
            tk       = test.T(:,k);
            
            % get LR image
            LR_k      = upscale(lk);
            VO_k      = upscale(vk);
            T_k       = upscale(tk);
            
            % get RAW-SHRM image
            RAW_SHRM_k = vec(compute_shrm(lk,raw_data_hr_means{i,j}.(day_or_night).sm,raw_data_hr_means{i,j}.(day_or_night).unobserved,pm{i,j}.(day_or_night),params));
            
            % get MC-SHRM image
            MC_SHRM_k = vec(compute_shrm(lk,mc_hr_means{i,j}.(day_or_night).sm,mc_hr_means{i,j}.(day_or_night).unobserved,pm{i,j}.(day_or_night),params));
            
            % get SR prediction
            SR_k = superresolve_LR_image_blockwise(lk,vk,tk,EOF{i,j}.(day_or_night),num_components,R2_threshold,params);
            
            % get missing pixel matrix
            m = mk & not(isnan(LR_k)) & not(isnan(MC_SHRM_k)) & not(isnan(RAW_SHRM_k)) & not(isnan(SR_k));
            
            % put data into composite image
            if (i==1)&&(j==1)
                data{cell_index}.LR                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                data{cell_index}.VO                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                data{cell_index}.T                      = NaN(num_pixels*num_rows,num_pixels*num_cols);
                data{cell_index}.RAW_SHRM               = NaN(num_pixels*num_rows,num_pixels*num_cols);
                data{cell_index}.MC_SHRM                = NaN(num_pixels*num_rows,num_pixels*num_cols);
                data{cell_index}.SR                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                data{cell_index}.HR                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                data{cell_index}.M                      = false(num_pixels*num_rows,num_pixels*num_cols);
                %
                data{cell_index}.LR(index_vector)       = LR_k(:);
                data{cell_index}.VO(index_vector)       = VO_k(:);
                data{cell_index}.T(index_vector)        = T_k(:);
                data{cell_index}.RAW_SHRM(index_vector) = RAW_SHRM_k(:);
                data{cell_index}.MC_SHRM(index_vector)  = MC_SHRM_k(:);
                data{cell_index}.SR(index_vector)       = SR_k(:);
                data{cell_index}.HR(index_vector)       = yk(:);
                data{cell_index}.M(index_vector)        = m(:);
                %
                data{cell_index}.measurement_time       = test.t_smap(k);
                data{cell_index}.sentinel_time          = test.t_sent(k);
                measurement_times(cell_index)           = test.t_smap(k);
                %
                cell_index = cell_index + 1;
            else
                % find corresponding index
                ind = find(  abs(test.t_smap(k) - measurement_times)< 2/24 );
                if isempty(ind)
                    data{cell_index}.LR                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                    data{cell_index}.VO                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                    data{cell_index}.T                      = NaN(num_pixels*num_rows,num_pixels*num_cols);
                    data{cell_index}.RAW_SHRM               = NaN(num_pixels*num_rows,num_pixels*num_cols);
                    data{cell_index}.MC_SHRM                = NaN(num_pixels*num_rows,num_pixels*num_cols);
                    data{cell_index}.SR                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                    data{cell_index}.HR                     = NaN(num_pixels*num_rows,num_pixels*num_cols);
                    data{cell_index}.M                      = false(num_pixels*num_rows,num_pixels*num_cols);
                    %
                    data{cell_index}.LR(index_vector)       = LR_k(:);
                    data{cell_index}.VO(index_vector)       = VO_k(:);
                    data{cell_index}.T(index_vector)        = T_k(:);
                    data{cell_index}.RAW_SHRM(index_vector) = RAW_SHRM_k(:);
                    data{cell_index}.MC_SHRM(index_vector)  = MC_SHRM_k(:);
                    data{cell_index}.SR(index_vector)       = SR_k(:);
                    data{cell_index}.HR(index_vector)       = yk(:);
                    data{cell_index}.M(index_vector)        = m(:);
                    %
                    data{cell_index}.measurement_time       = test.t_smap(k);
                    data{cell_index}.sentinel_time          = test.t_smap(k);
                    measurement_times(cell_index)           = test.t_smap(k);
                    %
                    cell_index = cell_index + 1;
                elseif length(ind)==1
                    data{ind}.LR(index_vector)       = LR_k(:);
                    data{ind}.VO(index_vector)       = VO_k(:);
                    data{ind}.T(index_vector)        = T_k(:);
                    data{ind}.RAW_SHRM(index_vector) = RAW_SHRM_k(:);
                    data{ind}.MC_SHRM(index_vector)  = MC_SHRM_k(:);
                    data{ind}.SR(index_vector)       = SR_k(:);
                    data{ind}.HR(index_vector)       = yk(:);
                    data{ind}.M(index_vector)        = m(:);
                else
                    error('something strange happened');
                end
            end
        end
    end
end

% order data by measurement time
[measurement_times,inds] = sort(measurement_times);
data                     = data(inds); 
end
%--------------------------------------------------------------------------

%==========================================================================
function data = truncate_unrealistic_values(data,params)
%TRUNCATE_UNREALISTIC_VALUES
fields     = {'RAW_SHRM','MC_SHRM','SR'};
%
num_obs    = length(data);
num_fields = length(fields);
for k=1:num_obs
    for j=1:num_fields
        data{k}.(fields{j})(data{k}.(fields{j})<params.MIN_REASONABLE_SM) = params.MIN_REASONABLE_SM;
        data{k}.(fields{j})(data{k}.(fields{j})>params.MAX_REASONABLE_SM) = params.MAX_REASONABLE_SM;
    end
end
end
%--------------------------------------------------------------------------

%==========================================================================
function calculate_errors_and_plot_results(data,day_or_night,measurement_times,time_offset,inds_to_plot,cmin,cmax,datapath)
%CALCULATE_ERRORS_AND_PLOT_RESULTS
switch lower(day_or_night)
    case {'d','day'}
        s = 'D';
    case {'n','a','night'}
        s = 'A';
    otherwise
        error('invalid option');
end

% define error metrics
RMSD   = @(x,y,m) sqrt( (1/length(x(m)))*sum((x(m) - y(m)).^2) );                                                        % root mean square difference / error
MD     = @(x,y,m) (1/length(x(m)))*sum(x(m)-y(m));                                                                       % mean difference or bias
ubRMSD = @(x,y,m) sqrt( RMSD(x,y,m)^2 - MD(x,y,m)^2 );                                                                   % unbiased RMSD
R      = @(x,y,m) sum((x(m)-mean(x(m))).*(y(m)-mean(y(m)))) / sqrt(sum((x(m)-mean(x(m))).^2)*sum((y(m)-mean(y(m))).^2)); % Pearson correlation
MAD    = @(x,y,m) 1/length(x(m))*sum( abs(x(m) - y(m)) );                                                                % mean absolute difference

% compute errors on full set of composite test set images
num_composite_measurements = length(measurement_times);
errors = cell(num_composite_measurements,1);
%
min_LR       = 1;
min_RAW_SHRM = 1;
min_MC_SHRM  = 1;
min_SR       = 1;
min_HR       = 1;
max_LR       = 0;
max_RAW_SHRM = 0;
max_MC_SHRM  = 0;
max_SR       = 0;
max_HR       = 0;
for k=1:num_composite_measurements
    % read out data
    LR_k       = data{k}.LR;
    RAW_SHRM_k = data{k}.RAW_SHRM;
    MC_SHRM_k  = data{k}.MC_SHRM;
    SR_k       = data{k}.SR;
    HR_k       = data{k}.HR;
    %
    M_k        = data{k}.M;
    
    % get min and max values
    min_LR       = min([min_LR,       min(min(LR_k))]);
    min_RAW_SHRM = min([min_RAW_SHRM, min(min(RAW_SHRM_k))]);
    min_MC_SHRM  = min([min_MC_SHRM,  min(min(MC_SHRM_k))]);
    min_SR       = min([min_SR,       min(min(SR_k))]);
    min_HR       = min([min_HR,       min(min(HR_k))]);
    %
    max_LR       = max([max_LR,       max(max(LR_k))]);
    max_RAW_SHRM = max([max_RAW_SHRM, max(max(RAW_SHRM_k))]);
    max_MC_SHRM  = max([max_MC_SHRM,  max(max(MC_SHRM_k))]);
    max_SR       = max([max_SR,       max(max(SR_k))]);
    max_HR       = max([max_HR,       max(max(HR_k))]);
    
    % compute errors
    errors{k}.Sent.LR.RMSD             = RMSD(LR_k,HR_k,M_k);
    errors{k}.Sent.LR.MD               = MD(LR_k,HR_k,M_k);
    errors{k}.Sent.LR.ubRMSD           = ubRMSD(LR_k,HR_k,M_k);
    errors{k}.Sent.LR.R                = R(LR_k,HR_k,M_k);
    errors{k}.Sent.LR.MAD              = MAD(LR_k,HR_k,M_k);
    %
    errors{k}.Sent.RAW_SHRM.RMSD       = RMSD(RAW_SHRM_k,HR_k,M_k);
    errors{k}.Sent.RAW_SHRM.MD         = MD(RAW_SHRM_k,HR_k,M_k);
    errors{k}.Sent.RAW_SHRM.ubRMSD     = ubRMSD(RAW_SHRM_k,HR_k,M_k);
    errors{k}.Sent.RAW_SHRM.R          = R(RAW_SHRM_k,HR_k,M_k);
    errors{k}.Sent.RAW_SHRM.MAD        = MAD(RAW_SHRM_k,HR_k,M_k);
    %
    errors{k}.Sent.MC_SHRM.RMSD        = RMSD(MC_SHRM_k,HR_k,M_k);
    errors{k}.Sent.MC_SHRM.MD          = MD(MC_SHRM_k,HR_k,M_k);
    errors{k}.Sent.MC_SHRM.ubRMSD      = ubRMSD(MC_SHRM_k,HR_k,M_k);
    errors{k}.Sent.MC_SHRM.R           = R(MC_SHRM_k,HR_k,M_k);
    errors{k}.Sent.MC_SHRM.MAD         = MAD(MC_SHRM_k,HR_k,M_k);
    %
    errors{k}.Sent.SR.RMSD             = RMSD(SR_k,HR_k,M_k);
    errors{k}.Sent.SR.MD               = MD(SR_k,HR_k,M_k);
    errors{k}.Sent.SR.ubRMSD           = ubRMSD(SR_k,HR_k,M_k);
    errors{k}.Sent.SR.R                = R(SR_k,HR_k,M_k);
    errors{k}.Sent.SR.MAD              = MAD(SR_k,HR_k,M_k);
    %
    errors{k}.Sent.LR.num_pixels       = sum(sum(M_k));
    errors{k}.Sent.RAW_SHRM.num_pixels = sum(sum(M_k));
    errors{k}.Sent.MC_SHRM.num_pixels  = sum(sum(M_k));
    errors{k}.Sent.SR.num_pixels       = sum(sum(M_k));
end

% get min and max of colorbar range
cmin_data = min([min_LR,min_RAW_SHRM,min_MC_SHRM,min_SR,min_HR]);
cmax_data = max([max_LR,max_RAW_SHRM,max_MC_SHRM,max_SR,max_HR]);

if cmax_data > cmax
    warning(['Data has max value=',num2str(cmax_data),...
        ' but plotting with cmax=',num2str(cmax)]);
end
if cmin_data < cmin
    warning(['Data has min value=',num2str(cmin_data),...
        ' but plotting with cmin=',num2str(cmin)]);
end

% plot composite test set images with error values displayed
load([datapath,'/SMAP_Color_SoilMoisture.mat']); % get colormap
if inds_to_plot==-1
    num_figs = ceil(num_composite_measurements/5);
else
    num_figs = ceil(length(inds_to_plot)/5);
end
for n=1:num_figs
    figure
    for k=1:5
        if inds_to_plot==-1
            index = (n-1)*5 + k;
        else
            if k>length(inds_to_plot)
                break;
            else
                index = inds_to_plot(k);
            end
        end
        
        if index>num_composite_measurements
            break;
        end

        % LR (SCA)
        subplot_ind = k;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        imagesc(data{index}.LR), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title({['SMAP SCA (',s,')'],datestr(data{index}.measurement_time-time_offset)},'interpreter','latex')

        % RAW-SHRM
        subplot_ind = k+5;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        imagesc(data{index}.RAW_SHRM), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('Raw Data Climatology','interpreter','latex')

        % MC-SHRM
        subplot_ind = k+10;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        imagesc(data{index}.MC_SHRM), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('MC Climatology','interpreter','latex')

        % SR
        subplot_ind = k+15;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        imagesc(data{index}.SR), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('EOF Superresolution','interpreter','latex')

        % HR
        subplot_ind = k+20;
        ax(subplot_ind) = subplot(5,5,subplot_ind);
        imagesc(data{index}.HR), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('SMAP-Sentinel','interpreter','latex')

        % display errors
        % LR SCA
        r = 0.164 * mod(k-1,5);
        annotation('textbox',  [0.155+r, 0.707,  0.1, 0.1], 'String',  "RMSD   = " + errors{index}.Sent.LR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.65,  0.1, 0.1], 'String',   "MD     = " + errors{index}.Sent.LR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.696, 0.1, 0.1], 'String',   "ubRMSD = " + errors{index}.Sent.LR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.684,  0.1, 0.1], 'String', "R      = " + errors{index}.Sent.LR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        % 0.728, 0.714

        % Raw Data SHRM
        annotation('textbox',  [0.155+r, 0.534,  0.1, 0.1], 'String', "RMSD   = " + errors{index}.Sent.RAW_SHRM.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.68,  0.1, 0.1], 'String',  "MD     = " + errors{index}.Sent.RAW_SHRM.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.522, 0.1, 0.1], 'String',  "ubRMSD = " + errors{index}.Sent.RAW_SHRM.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.510,  0.1, 0.1], 'String',  "R      = " + errors{index}.Sent.RAW_SHRM.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %0.585, 0.571

        % MC SHRM
        annotation('textbox',  [0.155+r, 0.362,  0.1, 0.1], 'String', "RMSD   = " + errors{index}.Sent.MC_SHRM.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.68,  0.1, 0.1], 'String',  "MD     = " + errors{index}.Sent.MC_SHRM.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.350, 0.1, 0.1], 'String',   "ubRMSD = " + errors{index}.Sent.MC_SHRM.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.338,  0.1, 0.1], 'String',  "R      = " + errors{index}.Sent.MC_SHRM.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        % 0.443, 0.429

        % EOF SR
        annotation('textbox',  [0.155+r, 0.19,  0.1, 0.1], 'String', "RMSD   = " + errors{index}.Sent.SR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.155+r, 0.68,  0.1, 0.1], 'String',  "MD     = " + errors{index}.Sent.SR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.178,  0.1, 0.1], 'String', "ubRMSD = " + errors{index}.Sent.SR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.155+r, 0.166,  0.1, 0.1], 'String',  "R      = " + errors{index}.Sent.SR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        % 0.3, 0.286
    end

    % insert common colorbar
    h=colorbar; h.TickLabelInterpreter='latex'; colormap(uint8(SMAP_Color_Scale));
    caxis([cmin,cmax])
    set(h, 'Position', [0.9, 0.15, 0.02, 0.7])
    clear ax
end
end
%--------------------------------------------------------------------------

%==========================================================================
function [] = perform_EOF_interpolation_experiment(datapath,vex,EOF,data,R2_thresh,params)
%PERFORM_EOF_INTERPOLATION_EXPERIMENT
missing_percent        = [0,0.1,0.2,0.3,0.4,0.5,0.6];
num_missingness_levels = length(missing_percent);

% get time offset
switch lower(vex)
    case {'ia','i','smapvex16i'}
        vex_label        = 'SMAPVEX16-Iowa';
        time_offset      = 5;
        measurement_inds = [4,3,4,2,6,5];
        day_or_night     = {'night','night','day','night','day','night'};
        grid_size_x      = 72;
        grid_size_y      = 72;
    case {'ma','m','smapvex16m'}
        vex_label        = 'SMAPVEX16-Manitoba';
        time_offset      = 5;
        measurement_inds = [6,7,5,4,2,7];
        day_or_night     = {'day','night','day','night','night','day'};
        grid_size_x      = 72;
        grid_size_y      = 72;
    case {'az','a','smapvex15a'}
        vex_label        = 'SMAPVEX15-Arizona';
        time_offset      = 7;
        measurement_inds = [4,2,1,3,6,5];
        day_or_night     = {'day','day','day','day','day','day'};
        grid_size_x      = 144;
        grid_size_y      = 108;
end
num_measurements = length(measurement_inds);
cmin = 0.02;
cmax = 0.7624;

% get common colorbar
min_HR = 1;
max_HR = 0;
for k=1:num_measurements
    HR_k   = data.(day_or_night{k}){measurement_inds(k)}.HR;
    min_HR = min([min_HR, min(min(HR_k))]);
    max_HR = max([max_HR, max(max(HR_k))]);
end

% randomly delete pixels at different levels
M_original = cell(num_measurements,1);
M          = cell(num_measurements,1);
for k=1:num_measurements
    M_original{k}       = not(isnan(data.(day_or_night{k}){measurement_inds(k)}.HR));
    num_original_pixels = sum(sum(M_original{k}));
    M{k}                = cell(length(missing_percent),1);
    for j=1:length(missing_percent)
        num_pixels_to_keep = floor(num_original_pixels*(1-missing_percent(j)));
        valid_inds         = find(M_original{k});
        inds_to_keep       = randperm(num_original_pixels,num_pixels_to_keep);
        %
        im                           = false(size(M_original{k}));
        im(valid_inds(inds_to_keep)) = true;
        M{k}{j}                      = im;
    end
end

% plot images at different levels of missingness
load([datapath,'/SMAP_Color_SoilMoisture.mat']); % get colormap
figure
missing_images = cell(num_missingness_levels,1);
for j=1:num_missingness_levels
    missing_images{j} = cell(num_measurements,1);
    for k=1:num_measurements
        HR_k = data.(day_or_night{k}){measurement_inds(k)}.HR;
        missing_images{j}{k}               = HR_k;
        missing_images{j}{k}(not(M{k}{j})) = NaN;
        %
        ind = (j-1)*num_measurements + k;
        subplot(num_missingness_levels,num_measurements,ind)
        imagesc(missing_images{j}{k}), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        if ind<=num_measurements
            title({['SMAP-Sentinel ',num2str(k)],...
                datestr(data.(day_or_night{k}){measurement_inds(k)}.measurement_time - time_offset/24),...
                [num2str(missing_percent(j)*100),'\% Missing']},'interpreter','latex')
        else
            title([num2str(missing_percent(j)*100),'\% Missing'],'interpreter','latex')
        end
    end
end

% insert common colorbar
h=colorbar; h.TickLabelInterpreter='latex'; colormap(uint8(SMAP_Color_Scale));
caxis([cmin,cmax])
set(h, 'Position', [0.93, 0.15, 0.02, 0.7])
clear ax

% compute different interpolations
% - interpolate missing pixels with corresponding SMAP pixel values
LR_interp = cell(num_measurements,1);
for k=1:num_measurements
    base_LR_im   = data.(day_or_night{k}){measurement_inds(k)}.LR .* M{k}{1};
    perm_missing = not(M{k}{1});
    LR_interp{k} = cell(num_missingness_levels-1,1);
    for j=2:num_missingness_levels
        fill_pixels                  = not(M{k}{1} & M{k}{j}) & not(perm_missing);
        LR_interp{k}{j}              = missing_images{1}{k};
        LR_interp{k}{j}(fill_pixels) = base_LR_im(fill_pixels);
    end
end

% - interpolate missing pixels using bicubic interpolation
BC_interp = cell(num_measurements,1);
for k=1:num_measurements
    perm_missing = not(M{k}{1});
    BC_interp{k} = cell(num_missingness_levels-1,1);
    for j=2:num_missingness_levels
        fill_pixels     = not(M{k}{1} & M{k}{j}) & not(perm_missing);
        BC_interp{k}{j} = missing_images{j}{k};
        v               = missing_images{j}{k};
        
        % remove missing samples
        [x,y]           = meshgrid(1:grid_size_x,1:grid_size_y);
        x(not(M{k}{j})) = [];
        y(not(M{k}{j})) = [];
        v(not(M{k}{j})) = [];
        
        % specify sampling points
        [xq,yq] = meshgrid(1:grid_size_x,1:grid_size_y);
        xq = xq(fill_pixels);
        yq = yq(fill_pixels);
        
        % perform cubic interpolation
        vq = griddata(x,y,v,xq,yq,'cubic');
        
        % insert interpolated values
        vq(vq<0.02) = 0.02;
        vq(vq>1)    = 1;
        BC_interp{k}{j}(fill_pixels) = vq;
    end
end

% - interpolate pixels using EOFs
num_subpixels = 9;
num_pixels    = 4*num_subpixels;
EOF_interp    = cell(num_measurements,1);
for k=1:num_measurements
    perm_missing  = not(M{k}{1});
    EOF_interp{k} = cell(num_missingness_levels-1,1);
    for m=2:num_missingness_levels
        fill_pixels      = not(M{k}{1} & M{k}{m}) & not(perm_missing);
        EOF_interp{k}{m} = missing_images{m}{k};
        for r=1:size(EOF,1)
            for c=1:size(EOF,2)
                for i=1:size(EOF{1}.(day_or_night{k}).E,1)
                    for j=1:size(EOF{1}.(day_or_night{k}).E,2)
                        % get HR index vector to index pixels in block (i,j)
                        hr_ind_vec = zeros(size(perm_missing));
                        hr_ind_vec( ((r-1)*num_pixels+(i-1)*num_subpixels+1):((r-1)*num_pixels+i*num_subpixels),((c-1)*num_pixels+(j-1)*num_subpixels+1):((c-1)*num_pixels+j*num_subpixels)) = 1;
                        hr_ind_vec = logical(hr_ind_vec(:));
                        
                        % check if there are enough pixels to fill-in here
                        if sum(sum(fill_pixels(hr_ind_vec)))==0
                            continue;
                        end
                        fill_in_pixels = fill_pixels(hr_ind_vec);
                        
                        % get the EOFs for this 9km EASE cell
                        E = EOF{r,c}.(day_or_night{k}).E{i,j};
                        if isempty(E)
                            continue;
                        end
                        
                        % get the target image
                        hr      = missing_images{m}{k}(hr_ind_vec);
                        obs_pix = not(isnan(hr) | isnan(E(:,1)));
                        
                        % compute missing pixels and return gap-filled
                        % image
                        hr_gap_filled = gap_fill_with_eofs(EOF{r,c}.(day_or_night{k}),...
                            i,j,R2_thresh,hr,obs_pix,missing_percent(m),...
                            fill_in_pixels,params);
                        
                        % use pixels from reconstruction to fill-in missing
                        % pixels (EOF interpolation)
                        EOF_interp{k}{m}(hr_ind_vec) = hr_gap_filled(:);
                    end
                end
            end
        end
    end
end

% - fill in missing pixels using SR prediction
SR_interp = cell(num_measurements,1);
for k=1:num_measurements
    base_SR_im   = data.(day_or_night{k}){measurement_inds(k)}.SR .* M{k}{1};
    perm_missing = not(M{k}{1});
    SR_interp{k} = cell(num_missingness_levels-1,1);
    for j=2:num_missingness_levels
        fill_pixels                  = not(M{k}{1} & M{k}{j}) & not(perm_missing);
        SR_interp{k}{j}              = missing_images{1}{k};
        pred                         = base_SR_im(fill_pixels);
        pred(pred<0.02)              = 0.02;
        pred(pred>1)                 = 1;
        SR_interp{k}{j}(fill_pixels) = pred;
    end
end

% define error metrics
RMSD   = @(x,y,m) sqrt( (1/length(x(m)))*sum((x(m) - y(m)).^2) );                                                        % root mean square difference / error
MD     = @(x,y,m) (1/length(x(m)))*sum(x(m)-y(m));                                                                       % mean difference or bias
ubRMSD = @(x,y,m) sqrt( RMSD(x,y,m)^2 - MD(x,y,m)^2 );                                                                   % unbiased RMSD
MAD    = @(x,y,m) 1/length(x(m))*sum( abs(x(m) - y(m)) );                                                                % mean absolute difference
R      = @(x,y,m) sum((x(m)-mean(x(m))).*(y(m)-mean(y(m)))) / sqrt(sum((x(m)-mean(x(m))).^2)*sum((y(m)-mean(y(m))).^2)); % Pearson correlation

% compute errors on full set of composite test set images
errors         = cell(num_measurements,1);
cmin_matrix    = ones(num_measurements,num_missingness_levels);
cmax_matrix    = ones(num_measurements,num_missingness_levels);
%
min_LR_INTERP  = 1;
min_BC_INTERP  = 1;
min_SR_INTERP  = 1;
min_EOF_INTERP = 1;
min_HR         = 1;
max_LR_INTERP  = 0;
max_BC_INTERP  = 0;
max_SR_INTERP  = 0;
max_EOF_INTERP = 0;
max_HR         = 0;
for k=1:num_measurements
    perm_missing = not(M{k}{1});
    errors{k}    = cell(num_missingness_levels,1);
    for j=2:num_missingness_levels
        % read out images
        LR_k  = LR_interp{k}{j};
        BC_k  = BC_interp{k}{j};
        SR_k  = SR_interp{k}{j};
        EOF_k = EOF_interp{k}{j};
        HR_k  = data.(day_or_night{k}){measurement_inds(k)}.HR;
        
        % get pattern of filled-in (originally missing) pixels
        M_k = not(M{k}{1} & M{k}{j}) & not(perm_missing);
        M_k = M_k & not(isnan(LR_k.*BC_k.*SR_k.*EOF_k.*HR_k));
        
        % get min and max values
        min_LR_INTERP  = min([min_LR_INTERP,  min(min(LR_k))]);
        min_BC_INTERP  = min([min_BC_INTERP,  min(min(BC_k))]);
        min_SR_INTERP  = min([min_SR_INTERP,  min(min(SR_k))]);
        min_EOF_INTERP = min([min_EOF_INTERP, min(min(EOF_k))]);
        min_HR         = min([min_HR,         min(min(HR_k))]);
        %
        max_LR_INTERP  = max([max_LR_INTERP,  max(max(LR_k))]);
        max_BC_INTERP  = max([max_BC_INTERP,  max(max(BC_k))]);
        max_SR_INTERP  = max([max_SR_INTERP,  max(max(SR_k))]);
        max_EOF_INTERP = max([max_EOF_INTERP, max(max(EOF_k))]);
        max_HR         = max([max_HR,         max(max(HR_k))]);
        
        % compute errors
        errors{k}{j}.LR.RMSD        = RMSD(LR_k,HR_k,M_k);
        errors{k}{j}.LR.MD          = MD(LR_k,HR_k,M_k);
        errors{k}{j}.LR.ubRMSD      = ubRMSD(LR_k,HR_k,M_k);
        errors{k}{j}.LR.MAD         = MAD(LR_k,HR_k,M_k);
        errors{k}{j}.LR.R           = R(LR_k,HR_k,M_k);
        %
        errors{k}{j}.BC.RMSD        = RMSD(BC_k,HR_k,M_k);
        errors{k}{j}.BC.MD          = MD(BC_k,HR_k,M_k);
        errors{k}{j}.BC.ubRMSD      = ubRMSD(BC_k,HR_k,M_k);
        errors{k}{j}.BC.MAD         = MAD(BC_k,HR_k,M_k);
        errors{k}{j}.BC.R           = R(BC_k,HR_k,M_k);
        %
        errors{k}{j}.SR.RMSD        = RMSD(SR_k,HR_k,M_k);
        errors{k}{j}.SR.MD          = MD(SR_k,HR_k,M_k);
        errors{k}{j}.SR.ubRMSD      = ubRMSD(SR_k,HR_k,M_k);
        errors{k}{j}.SR.MAD         = MAD(SR_k,HR_k,M_k);
        errors{k}{j}.SR.R           = R(SR_k,HR_k,M_k);
        %
        errors{k}{j}.EOF.RMSD       = RMSD(EOF_k,HR_k,M_k);
        errors{k}{j}.EOF.MD         = MD(EOF_k,HR_k,M_k);
        errors{k}{j}.EOF.ubRMSD     = ubRMSD(EOF_k,HR_k,M_k);
        errors{k}{j}.EOF.MAD        = MAD(EOF_k,HR_k,M_k);
        errors{k}{j}.EOF.R          = R(EOF_k,HR_k,M_k);
        %
        errors{k}{j}.LR.num_pixels  = sum(sum(M_k));
        errors{k}{j}.BC.num_pixels  = sum(sum(M_k));
        errors{k}{j}.SR.num_pixels  = sum(sum(M_k));
        errors{k}{j}.EOF.num_pixels = sum(sum(M_k));
        
        % get min and max of colorbar range
        cmin_matrix(k,j) = min([min_LR_INTERP,min_BC_INTERP,min_SR_INTERP,min_EOF_INTERP,min_HR]); 
        cmax_matrix(k,j) = max([max_LR_INTERP,max_BC_INTERP,max_SR_INTERP,max_EOF_INTERP,max_HR]);
    end
end

% display interpolation results
% plot composite test set images with error values displayed
load([datapath,'/SMAP_Color_SoilMoisture.mat']); % get colormap
for j=2:num_missingness_levels
    figure
    for k=1:num_measurements
        % HR image with Missing Pixels
        subplot_ind = k;
        ax(subplot_ind) = subplot(6,num_measurements,subplot_ind);
        imagesc(missing_images{j}{k}), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title({vex_label,'SMAP-Sentinel',...
                datestr(data.(day_or_night{k}){measurement_inds(k)}.measurement_time - time_offset/24),...
                [num2str(missing_percent(j)*100),'\% Missing']},'interpreter','latex')
        
        % LR Interp
        subplot_ind = num_measurements + k;
        ax(subplot_ind) = subplot(6,num_measurements,subplot_ind);
        imagesc(LR_interp{k}{j}), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('SMAP Interpolation','interpreter','latex')
        
        % BC Interp
        subplot_ind = 2*num_measurements + k;
        ax(subplot_ind) = subplot(6,num_measurements,subplot_ind);
        imagesc(BC_interp{k}{j}), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('Bicubic Interpolation','interpreter','latex')
        
        % SR Interp
        subplot_ind = 3*num_measurements + k;
        ax(subplot_ind) = subplot(6,num_measurements,subplot_ind);
        imagesc(SR_interp{k}{j}), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('SR Interpolation','interpreter','latex')
        
        % EOF Interp
        subplot_ind = 4*num_measurements + k;
        ax(subplot_ind) = subplot(6,num_measurements,subplot_ind);
        imagesc(EOF_interp{k}{j}), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('EOF Interpolation','interpreter','latex')
        
        % GT HR Image
        subplot_ind = 5*num_measurements + k;
        ax(subplot_ind) = subplot(6,num_measurements,subplot_ind);
        imagesc(data.(day_or_night{k}){measurement_inds(k)}.HR), caxis manual, caxis([cmin,cmax])
        colormap(uint8(SMAP_Color_Scale));
        axis square, set(gca,'xtick',[],'ytick',[])
        title('Ground Truth','interpreter','latex')
        
        % display errors
        % LR Interp
        r = 0.135 * mod(k-1,6);
        annotation('textbox',  [0.15+r, 0.585,  0.1, 0.1], 'String',  "RMSD   = " + errors{k}{j}.LR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.15+r, 0.65,  0.1, 0.1], 'String',   "MD     = " + errors{k}{j}.LR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.15+r, 0.571, 0.1, 0.1], 'String',   "ubRMSD = " + errors{k}{j}.LR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox',  [0.15+r, 0.6825,  0.1, 0.1], 'String', "R      = " + errors{k}{j}.LR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        
        % BC Interp
        annotation('textbox',  [0.15+r, 0.443,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{j}.BC.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.15+r, 0.68,  0.1, 0.1], 'String',  "MD     = " + errors{k}{j}.BC.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.15+r, 0.429, 0.1, 0.1], 'String',  "ubRMSD = " + errors{k}{j}.BC.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox',  [0.15+r, 0.51,  0.1, 0.1], 'String',  "R      = " + errors{k}{j}.BC.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        
        % SR Interp
        annotation('textbox',  [0.15+r, 0.3,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{j}.SR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.15+r, 0.68,  0.1, 0.1], 'String',  "MD     = " + errors{k}{j}.SR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.15+r, 0.285, 0.1, 0.1], 'String',   "ubRMSD = " + errors{k}{j}.SR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox',  [0.15+r, 0.51,  0.1, 0.1], 'String',  "R      = " + errors{k}{j}.SR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        
        % EOF Interp
        annotation('textbox',  [0.15+r, 0.157,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{j}.EOF.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox', [0.15+r, 0.68,  0.1, 0.1], 'String',  "MD     = " + errors{k}{j}.EOF.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        annotation('textbox',  [0.15+r, 0.143,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{j}.EOF.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
        %annotation('textbox',  [0.15+r, 0.337,  0.1, 0.1], 'String',  "R      = " + errors{k}{j}.EOF.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
    end

    % insert common colorbar
    h=colorbar; h.TickLabelInterpreter='latex'; colormap(uint8(SMAP_Color_Scale));
    caxis([cmin,cmax])
    set(h, 'Position', [0.93, 0.15, 0.02, 0.7])
    clear ax
end

% plot single example in standalone figure
figure
%
k = 2; % measurement index
m = 2; % missing percentage index
%
% LR Interp
ax(2) = subplot(5,5,1);
imagesc(LR_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('SMAP Interpolation','interpreter','latex')

% BC Interp
ax(3) = subplot(5,5,2);
imagesc(BC_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('Bicubic Interpolation','interpreter','latex')

% SR Interp
ax(4) = subplot(5,5,3);
imagesc(SR_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('SR Interpolation','interpreter','latex')

% EOF Interp
ax(5) = subplot(5,5,4);
imagesc(EOF_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('EOF Interpolation','interpreter','latex')

% Ground Truth HR Image
% ax(6) = subplot(5,5,5);
% imagesc(data.(day_or_night{k}){measurement_inds(k)}.HR), caxis manual, caxis([cmin,cmax])
% colormap(uint8(SMAP_Color_Scale));
% axis square, set(gca,'xtick',[],'ytick',[])
% title('Ground Truth','interpreter','latex')
% title({vex_label,datestr(data.(day_or_night{k}){measurement_inds(k)}.measurement_time - time_offset/24),},'interpreter','latex')

% Missing Pixel Image
ax(6) = subplot(5,5,5);
imagesc(missing_images{m}{k}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title({'SMAP-Sentinel',[num2str(missing_percent(m)*100),'\% Missing']},'interpreter','latex')


% display errors
% LR Interp
annotation('textbox',  [0.160, 0.7000,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.LR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.160, 0.6850,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.LR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.160, 0.6850,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.LR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.160, 0.6700,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.LR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% BC Interp
annotation('textbox',  [0.320, 0.7000,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.BC.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.320, 0.6850,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.BC.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.320, 0.6850,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.BC.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.320, 0.6700,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.BC.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% SR Interp
annotation('textbox',  [0.485, 0.7000,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.SR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.485, 0.6850,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.SR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.485, 0.6850,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.SR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.485, 0.6700,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.SR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% EOF Interp
annotation('textbox',  [0.650, 0.7000,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.EOF.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.650, 0.6850,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.EOF.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.650, 0.6850,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.EOF.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.650, 0.6700,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.EOF.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)


k = 4; % measurement index
m = 4; % missingness index
% LR Interp
ax(8) = subplot(5,5,11);
imagesc(LR_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('SMAP Interpolation','interpreter','latex')

% BC Interp
ax(9) = subplot(5,5,12);
imagesc(BC_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('Bicubic Interpolation','interpreter','latex')

% SR Interp
ax(10) = subplot(5,5,13);
imagesc(SR_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('SR Interpolation','interpreter','latex')

% EOF Interp
ax(11) = subplot(5,5,14);
imagesc(EOF_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('EOF Interpolation','interpreter','latex')

% Ground Truth HR Image
% ax(12) = subplot(5,5,15);
% imagesc(data.(day_or_night{k}){measurement_inds(k)}.HR), caxis manual, caxis([cmin,cmax])
% colormap(uint8(SMAP_Color_Scale));
% axis square, set(gca,'xtick',[],'ytick',[])
% title('Ground Truth','interpreter','latex')
% title(datestr(data.(day_or_night{k}){measurement_inds(k)}.measurement_time - time_offset/24),'interpreter','latex')

% Missing Pixel Image
ax(12) = subplot(5,5,15);
imagesc(missing_images{m}{k}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title({[num2str(missing_percent(m)*100),'\% Missing']},'interpreter','latex')


% display errors
% LR Interp
annotation('textbox',  [0.160, 0.3600,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.LR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.160, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.LR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.160, 0.3450,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.LR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.160, 0.3300,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.LR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% BC Interp
annotation('textbox',  [0.320, 0.3600,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.BC.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.320, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.BC.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.320, 0.3450,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.BC.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.320, 0.3300,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.BC.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% SR Interp
annotation('textbox',  [0.485, 0.3600,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.SR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.485, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.SR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.485, 0.3450,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.SR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.485, 0.3300,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.SR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% EOF Interp
annotation('textbox',  [0.650, 0.3600,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.EOF.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.650, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.EOF.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.650, 0.3450,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.EOF.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.650, 0.3300,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.EOF.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)


k = 3; % measurement index
m = 6; % missingness index
% LR Interp
ax(13) = subplot(5,5,21);
imagesc(LR_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('SMAP Interpolation','interpreter','latex')

% BC Interp
ax(14) = subplot(5,5,22);
imagesc(BC_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('Bicubic Interpolation','interpreter','latex')

% SR Interp
ax(15) = subplot(5,5,23);
imagesc(SR_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('SR Interpolation','interpreter','latex')

% EOF Interp
ax(16) = subplot(5,5,24);
imagesc(EOF_interp{k}{m}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title('EOF Interpolation','interpreter','latex')

% Ground Truth HR Image
% ax(17) = subplot(5,5,25);
% imagesc(data.(day_or_night{k}){measurement_inds(k)}.HR), caxis manual, caxis([cmin,cmax])
% colormap(uint8(SMAP_Color_Scale));
% axis square, set(gca,'xtick',[],'ytick',[])
% title('Ground Truth','interpreter','latex')
% title(datestr(data.(day_or_night{k}){measurement_inds(k)}.measurement_time - time_offset/24),'interpreter','latex')

% Missing Pixel Image
ax(17) = subplot(5,5,25);
imagesc(missing_images{m}{k}), caxis manual, caxis([cmin,cmax])
colormap(uint8(SMAP_Color_Scale));
axis square, set(gca,'xtick',[],'ytick',[])
title({[num2str(missing_percent(m)*100),'\% Missing']},'interpreter','latex')

% display errors
% LR Interp
annotation('textbox',  [0.210, 0.0300,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.LR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.210, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.LR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.210, 0.0150,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.LR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.210, 0.0000,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.LR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% BC Interp
annotation('textbox',  [0.380, 0.0300,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.BC.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.380, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.BC.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.380, 0.0150,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.BC.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.380, 0.0000,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.BC.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% SR Interp
annotation('textbox',  [0.540, 0.0300,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.SR.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.540, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.SR.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.540, 0.0150,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.SR.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.540, 0.0000,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.SR.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% EOF Interp
annotation('textbox',  [0.700, 0.0300,  0.1, 0.1], 'String', "RMSD   = " + errors{k}{m}.EOF.RMSD,   'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
% annotation('textbox',[0.700, 0.3450,  0.1, 0.1], 'String', "MD     = " + errors{k}{m}.EOF.MD,     'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.700, 0.0150,  0.1, 0.1], 'String', "ubRMSD = " + errors{k}{m}.EOF.ubRMSD, 'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)
annotation('textbox',  [0.700, 0.0000,  0.1, 0.1], 'String', "R      = " + errors{k}{m}.EOF.R,      'LineStyle', 'none', 'interpreter', 'latex', 'fontsize', 8)

% insert common colorbar
h=colorbar; h.TickLabelInterpreter='latex'; colormap(uint8(SMAP_Color_Scale));
caxis([cmin,cmax])
set(h, 'Position', [0.9, 0.15, 0.02, 0.7])
clear ax

% compute table of aggregate errors for each missing value percentage
err.LR.RMSD    = zeros(num_missingness_levels,1);
err.LR.MD      = zeros(num_missingness_levels,1);
err.LR.ubRMSD  = zeros(num_missingness_levels,1);
err.LR.MAD     = zeros(num_missingness_levels,1);
err.LR.R       = zeros(num_missingness_levels,1);
%
err.BC.RMSD    = zeros(num_missingness_levels,1);
err.BC.MD      = zeros(num_missingness_levels,1);
err.BC.ubRMSD  = zeros(num_missingness_levels,1);
err.BC.MAD     = zeros(num_missingness_levels,1);
err.BC.R       = zeros(num_missingness_levels,1);
%
err.SR.RMSD    = zeros(num_missingness_levels,1);
err.SR.MD      = zeros(num_missingness_levels,1);
err.SR.ubRMSD  = zeros(num_missingness_levels,1);
err.SR.MAD     = zeros(num_missingness_levels,1);
err.SR.R       = zeros(num_missingness_levels,1);
%
err.EOF.RMSD   = zeros(num_missingness_levels,1);
err.EOF.MD     = zeros(num_missingness_levels,1);
err.EOF.ubRMSD = zeros(num_missingness_levels,1);
err.EOF.MAD     = zeros(num_missingness_levels,1);
err.EOF.R      = zeros(num_missingness_levels,1);
for j=2:num_missingness_levels
    for k=1:num_measurements
        err.LR.RMSD(j)   = err.LR.RMSD(j)   + errors{k}{j}.LR.RMSD;
        err.LR.MD(j)     = err.LR.MD(j)     + errors{k}{j}.LR.MD;
        err.LR.ubRMSD(j) = err.LR.ubRMSD(j) + errors{k}{j}.LR.ubRMSD;
        err.LR.MAD(j)    = err.LR.MAD(j)    + errors{k}{j}.LR.MAD;
        err.LR.R(j)      = err.LR.R(j)      + errors{k}{j}.LR.R;
        %
        err.BC.RMSD(j)   = err.BC.RMSD(j)   + errors{k}{j}.BC.RMSD;
        err.BC.MD(j)     = err.BC.MD(j)     + errors{k}{j}.BC.MD;
        err.BC.ubRMSD(j) = err.BC.ubRMSD(j) + errors{k}{j}.BC.ubRMSD;
        err.BC.MAD(j)    = err.BC.MAD(j)    + errors{k}{j}.BC.MAD;
        err.BC.R(j)      = err.BC.R(j)      + errors{k}{j}.BC.R;
        %
        err.SR.RMSD(j)   = err.SR.RMSD(j)   + errors{k}{j}.SR.RMSD;
        err.SR.MD(j)     = err.SR.MD(j)     + errors{k}{j}.SR.MD;
        err.SR.ubRMSD(j) = err.SR.ubRMSD(j) + errors{k}{j}.SR.ubRMSD;
        err.SR.MAD(j)    = err.SR.MAD(j)    + errors{k}{j}.SR.MAD;
        err.SR.R(j)      = err.SR.R(j)      + errors{k}{j}.SR.R;
        %
        err.EOF.RMSD(j)   = err.EOF.RMSD(j)   + errors{k}{j}.EOF.RMSD;
        err.EOF.MD(j)     = err.EOF.MD(j)     + errors{k}{j}.EOF.MD;
        err.EOF.ubRMSD(j) = err.EOF.ubRMSD(j) + errors{k}{j}.EOF.ubRMSD;
        err.EOF.MAD(j)    = err.EOF.MAD(j)    + errors{k}{j}.EOF.MAD;
        err.EOF.R(j)      = err.EOF.R(j)      + errors{k}{j}.EOF.R;
    end
    err.LR.RMSD(j)   = err.LR.RMSD(j)/num_measurements;
    err.LR.MD(j)     = err.LR.MD(j)/num_measurements;
    err.LR.ubRMSD(j) = err.LR.ubRMSD(j)/num_measurements;
    err.LR.MAD(j)    = err.LR.MAD(j)/num_measurements;
    err.LR.R(j)      = err.LR.R(j)/num_measurements;
    %
    err.BC.RMSD(j)   = err.BC.RMSD(j)/num_measurements;
    err.BC.MD(j)     = err.BC.MD(j)/num_measurements;
    err.BC.ubRMSD(j) = err.BC.ubRMSD(j)/num_measurements;
    err.BC.MAD(j)    = err.BC.MAD(j)/num_measurements;
    err.BC.R(j)      = err.BC.R(j)/num_measurements;
    %
    err.SR.RMSD(j)   = err.SR.RMSD(j)/num_measurements;
    err.SR.MD(j)     = err.SR.MD(j)/num_measurements;
    err.SR.ubRMSD(j) = err.SR.ubRMSD(j)/num_measurements;
    err.SR.MAD(j)    = err.SR.MAD(j)/num_measurements;
    err.SR.R(j)      = err.SR.R(j)/num_measurements;
    %
    err.EOF.RMSD(j)   = err.EOF.RMSD(j)/num_measurements;
    err.EOF.MD(j)     = err.EOF.MD(j)/num_measurements;
    err.EOF.ubRMSD(j) = err.EOF.ubRMSD(j)/num_measurements;
    err.EOF.MAD(j)    = err.EOF.MAD(j)/num_measurements;
    err.EOF.R(j)      = err.EOF.R(j)/num_measurements;
end

% print table of aggregate errors for each missing value percentage
interp_method = {'LR','BC','SR','EOF'};
disp(newline)
disp('Interpolation Errors:')
for j=2:num_missingness_levels
    disp(['  Missing Level ',num2str(100*missing_percent(j)),'%:'])
    for k=1:length(interp_method)
        disp(['    ',interp_method{k},' Interp Errors:'])
        disp(['      Average RMSD   = ',num2str(err.(interp_method{k}).RMSD(j))])
        disp(['      Average MD     = ',num2str(err.(interp_method{k}).MD(j))])
        disp(['      Average ubRMSD = ',num2str(err.(interp_method{k}).ubRMSD(j))])
        disp(['      Average MAD    = ',num2str(err.(interp_method{k}).MAD(j))])
        disp(['      Average R      = ',num2str(err.(interp_method{k}).R(j))])
    end
end

end
%--------------------------------------------------------------------------
