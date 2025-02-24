%plot_global_errorts
% Plots Figures 6 and 7 from "Low Rank Gap-Filling and Downscaling for SMAP 
% Soil Moisture Datasets" by Beale, Bras, and Romberg. Ecohydrology 2025.
sentinel_res  = 1;
max_vwc       = 3;
err_type      = 'avg_R'; % options: 'avg_RMSD', 'avg_ubRMSD', 'avg_R', 'avg_MAD'

%% get paths and parameters
paths  = get_paths(sentinel_res,max_vwc);
params = get_params(sentinel_res,max_vwc);

%% load global results struct
load_struct    = load([paths.GLOBAL_RESULTS_DIR,'global_results']);
global_results = load_struct.global_results;
clear load_struct

%% plot set of images showing global errors for different methods
plot_global_error_images(global_results,[err_type,'_36km'],params);

%% plot histograms of global errors for different methods
plot_global_error_histograms(global_results,sentinel_res,params);

%% Ancillary Functions
%==========================================================================
function paths = get_paths(sentinel_res,max_vwc)
%GET_PATHS
param_str                = [num2str(sentinel_res),'km',num2str(max_vwc),'VWC/'];
paths.GLOBAL_RESULTS_DIR = ['./data/Global_Error_Data/',param_str];
end
%--------------------------------------------------------------------------

%==========================================================================
function params = get_params(sentinel_res,max_vwc)
%GET_PARAMS
params.SENTINEL_RES                = sentinel_res; % Sentinel pixel resolution [km]
params.MAX_ACCEPTABLE_VWC_SENTINEL = max_vwc;      % remove SMAP-Sentinel pixels with VWC > this threshold [kg/m^2]
end
%--------------------------------------------------------------------------

%==========================================================================
function [] = plot_global_error_images(GR,ERR_TYPE,params)
%PLOT_GLOBAL_ERROR_IMAGES
day_or_night = {'day','night'};
switch lower(ERR_TYPE)
    case {'avg_rmsd_36km','rmsd','avg_rmsd','avg_rmsd_9km'}
        err_str = 'RMSD';
    case {'avg_ubrmsd_36km','ubrmsd','avg_ubrmsd','avg_ubrmsd_9km'}
        err_str = 'ubRMSD';
    case {'avg_r_36km','r','avg_r','avg_r_9km'}
        err_str = 'R';
    case {'avg_mad_36km','mad','avg_mad','avg_mad_9km'}
        err_str = 'MAD';
    case {'avg_md_36km','md','avg_md','avg_md_9km'}
        err_str = 'MD';
end
resvwc_str = [' (',num2str(params.SENTINEL_RES),'km, ',...
    num2str(params.MAX_ACCEPTABLE_VWC_SENTINEL),'\,$kg/m^2$ Max. VWC)'];

% get common colorbar across all baselines
cmin = Inf;
cmax = -Inf;
baselines = {'LR','BC','RAW_SHRM','MC_SHRM'};
for k=1:length(day_or_night)
    for j=1:length(baselines)
        non_nan_inds = not(isnan(GR.images.(day_or_night{k}).(baselines{j}).(ERR_TYPE)));
        min_err      = min(GR.images.(day_or_night{k}).(baselines{j}).(ERR_TYPE)(non_nan_inds));
        max_err      = max(GR.images.(day_or_night{k}).(baselines{j}).(ERR_TYPE)(non_nan_inds));
        cmin         = min([cmin,min_err]);
        cmax         = max([cmax,max_err]);
    end
end
for k=1:length(day_or_night)
    non_nan_inds = not(isnan(GR.images.(day_or_night{k}).SR.(ERR_TYPE)));
    min_err      = min(GR.images.(day_or_night{k}).SR.(ERR_TYPE)(non_nan_inds));
    max_err      = max(GR.images.(day_or_night{k}).SR.(ERR_TYPE)(non_nan_inds));
    cmin         = min([cmin,min_err]);
    cmax         = max([cmax,max_err]);
end

% get lat/lon arrays
[latvec,lonvec] = ease_inverse(36,1:406,1:964);
[lon,lat]       = meshgrid(lonvec,latvec);

% plot images
for k=1:length(day_or_night)
    res = GR.images.(day_or_night{k}).SR.(ERR_TYPE);
        
    % SR error image displayed on world map
    figure
    h=worldmap('World');
    load coastlines;
    hold on
    geoshow(lat,lon,double(res),'DisplayType','texturemap');
    set(findall(h,'Tag','MLabel'),'visible','off')
    plotm(coastlat,coastlon,'LineWidth',1);
    title(['EOF ',upper(day_or_night{k}(1)),day_or_night{k}(2:end),...
        ' Method Average ',err_str,resvwc_str],...
        'interpreter','latex')
    h=colorbar; h.TickLabelInterpreter='latex';
    caxis manual, caxis([cmin,cmax]);
    
    figure
    h=worldmap('World'); load coastlines; hold on;
    geoshow(lat,lon,double(GR.images.(day_or_night{k}).LR.(ERR_TYPE)),'DisplayType','texturemap');
    set(findall(h,'Tag','MLabel'),'visible','off')
    plotm(coastlat,coastlon,'LineWidth',1);
    title(['SMAP ',upper(day_or_night{k}(1)),day_or_night{k}(2:end),' Data Average ',err_str,resvwc_str],'interpreter','latex')
    h=colorbar; h.TickLabelInterpreter='latex';
    caxis manual, caxis([cmin,cmax]);
    %
    figure
    h=worldmap('World'); load coastlines; hold on;
    geoshow(lat,lon,double(GR.images.(day_or_night{k}).RAW_SHRM.(ERR_TYPE)),'DisplayType','texturemap');
    set(findall(h,'Tag','MLabel'),'visible','off')
    plotm(coastlat,coastlon,'LineWidth',1);
    title(['Raw Data Climatology ',upper(day_or_night{k}(1)),day_or_night{k}(2:end),' Data Average ',err_str,resvwc_str],'interpreter','latex')
    h=colorbar; h.TickLabelInterpreter='latex';
    caxis manual, caxis([cmin,cmax]);
end

% plot global average error distribution image
res_day   = GR.images.day.SR.(ERR_TYPE);
res_night = GR.images.night.SR.(ERR_TYPE);
res_avg   = mean(cat(3,res_day,res_night),3,'omitnan');

figure
h=worldmap('World'); h.TickLabelInterpreter='latex';
load coastlines; hold on;
geoshow(lat,lon,double(res_avg),'DisplayType','texturemap');
set(findall(h,'Tag','MLabel'),'visible','off')
plotm(coastlat,coastlon,'LineWidth',1);
mlabel off; plabel off; gridm off;
title(['Global Average ',err_str,' Distribution for EOF Method (',num2str(params.SENTINEL_RES),'km Resolution)'],'interpreter','latex')
h=colorbar; h.TickLabelInterpreter='latex';
caxis manual; caxis([-1,1]);
end
%--------------------------------------------------------------------------

%==========================================================================
function [] = plot_global_error_histograms(GR,sentinel_res,params)
%PLOT_GLOBAL_ERROR_HISTOGRAMS
err_types = {'avg_RMSD_36km','avg_ubRMSD_36km','avg_MAD_36km','avg_R_36km'};
for j=1:length(err_types)
    switch lower(err_types{j})
        case {'avg_rmsd_36km'}
            err_str = 'RMSD';
            xmin    = 0;
            xmax    = 0.15;
            ymin    = 0;
            if sentinel_res==1
                ymax    = 5500;
            elseif sentinel_res==3
                ymax    = 3100;
            end
            edges   = linspace(0,0.15,80);
        case {'avg_ubrmsd_36km'}
            err_str = 'ubRMSD';
            xmin    = 0;
            xmax    = 0.14;
            ymin    = 0;
            if sentinel_res==1
                ymax    = 5000;
            elseif sentinel_res==3
                ymax    = 3100;
            end
            edges   = linspace(0,0.14,80);
        case {'avg_r_36km'}
            err_str = 'R';
            xmin    = -1;
            xmax    = 1;
            ymin    = 0;
            if sentinel_res==1
                ymax    = 5200;
            elseif sentinel_res==3
                ymax    = 2800;
            end
            edges   = linspace(-1,1,120);
        case {'avg_mad_36km'}
            err_str = 'MAD';
            xmin    = 0;
            xmax    = 0.15;
            ymin    = 0;
            if sentinel_res==1
                ymax    = 7100;
            elseif sentinel_res==3
                ymax    = 4000;
            end
            edges   = linspace(0,0.15,80);
    end
    
    figure
    % LR
    subplot(5,2,1)
    lr_data_day   = GR.images.day.LR.(err_types{j});
    lr_data_night = GR.images.night.LR.(err_types{j});
    lr_data = [lr_data_day(:);lr_data_night(:)];
    histogram(lr_data,edges)
    set(gca,'TickLabelInterpreter','latex')
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    hold on
    mean_value   = mean(lr_data(:),'omitnan');
    median_value = median(lr_data(:),'omitnan');
    plot([mean_value,mean_value],[ymin,ymax],'-r','LineWidth',2)
    plot([median_value,median_value],[ymin,ymax],'-b','LineWidth',2)
    title(['Average ',err_str,' Histogram for SMAP'],'interpreter','latex')
    xlabel('','interpreter','latex')
    ylabel('Number of Instances','interpreter','latex')
    legend({'',['Mean = ',num2str(mean_value)],['Median = ',num2str(median_value)]},'interpreter','latex')
    
    % RAW_SHRM
    subplot(5,2,3)
    rs_data_day   = GR.images.day.RAW_SHRM.(err_types{j});
    rs_data_night = GR.images.night.RAW_SHRM.(err_types{j});
    rs_data = [rs_data_day(:); rs_data_night(:)];
    histogram(rs_data,edges)
    set(gca,'TickLabelInterpreter','latex')
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    hold on
    mean_value   = mean(rs_data(:),'omitnan');
    median_value = median(rs_data(:),'omitnan');
    plot([mean_value,mean_value],[ymin,ymax],'-r','LineWidth',2)
    plot([median_value,median_value],[ymin,ymax],'-b','LineWidth',2)
    title(['Average ',err_str,' Histogram for Raw Data Climatology'],'interpreter','latex')
    xlabel('','interpreter','latex')
    ylabel('Number of Instances','interpreter','latex')
    legend({'',['Mean = ',num2str(mean_value)],['Median = ',num2str(median_value)]},'interpreter','latex')

    % SR
    subplot(5,2,5)
    sr_data_day   = GR.images.day.SR.(err_types{j});
        sr_data_night = GR.images.night.SR.(err_types{j});
    sr_data = [sr_data_day(:);sr_data_night(:)];
    histogram(sr_data,edges)
    set(gca,'TickLabelInterpreter','latex')
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    hold on
    mean_value   = mean(sr_data(:),'omitnan');
    median_value = median(sr_data(:),'omitnan');
    plot([mean_value,mean_value],[ymin,ymax],'-r','LineWidth',2)
    plot([median_value,median_value],[ymin,ymax],'-b','LineWidth',2)
    title(['Average ',err_str,' Histogram for EOF Method'],'interpreter','latex')
    xlabel('','interpreter','latex')
    ylabel('Number of Instances','interpreter','latex')
    legend({'',['Mean = ',num2str(mean_value)],['Median = ',num2str(median_value)]},'interpreter','latex')
end

end
%--------------------------------------------------------------------------
