function sr_output = superresolve_LR_image_blockwise(sm,vo,temp,eof,num_components,R2_threshold,params)
%SUPERRESOLVE_LR_IMAGES_BLOCKWISE Superresolve a low-resolution image using
%low-rank EOF superresolution method, applying method to each 9km x 9km
%EASE cell.
% - for each low-resolution measurement:
%   - for each 9km x 9km block in the area:
%       1. compute coefficients of EOFs using fitted parametric models
%       2. compute weighted linear combination of the EOFs using 
%          coefficients to obtain the superresolved image
num_blocks    = params.LR_SIDE_PIXELS;
num_pixels    = params.HR_SIDE_PIXELS;
num_subpixels = params.SMAP_RES/params.SENTINEL_RES;
%
sr_output = NaN(num_pixels,num_pixels);
for i=1:num_blocks
    for j=1:num_blocks
        % make sure we were able to fit curves for this block
        if isempty(eof.Beta{i,j})||isempty(eof.C{i,j})||isempty(eof.E{i,j})
            continue;
        end
        
        % get indexing vectors
        lr_ind_vec      = zeros(num_blocks,num_blocks);
        lr_ind_vec(i,j) = 1;
        lr_ind_vec      = logical(lr_ind_vec(:));
        %
        hr_ind_vec      = zeros(num_pixels,num_pixels);
        hr_ind_vec(((i-1)*num_subpixels+1):(i*num_subpixels),((j-1)*num_subpixels+1):(j*num_subpixels)) = 1;
        hr_ind_vec      = logical(hr_ind_vec(:));

        % get observed SMAP measurement for pixel (i,j)
        u = sm(lr_ind_vec);
        v = vo(lr_ind_vec);
        t = temp(lr_ind_vec);
        if isnan(u)||isnan(v)||isnan(t)
            continue;
        end
        
        % standardize variables
        u_mean = eof.Beta{i,j}.means.uk;
        v_mean = eof.Beta{i,j}.means.vk;
        t_mean = eof.Beta{i,j}.means.tk;
        %
        u_std  = eof.Beta{i,j}.std.uk;
        v_std  = eof.Beta{i,j}.std.vk;
        t_std  = eof.Beta{i,j}.std.tk;
        %
        us     = (u - u_mean)/u_std; % standardized mean soil moisture
        vs     = (v - v_mean)/v_std; % standardized mean vegetation opacity
        ts     = (t - t_mean)/t_std; % standardized mean surface temp
        
        % get regression matrix
        P = [1,us,vs,ts];
        
        % compute the EOF coefficients
        c = (P*eof.Beta{i,j}.curve_params)'; % predicted EOF coefficients from learned coefficient model

        % compute SR prediction
        SR = eof.mean_vector{i,j};
        if num_components==-1
            % use only EOFs whose coefficients have R^2 values > threshold
            valid_EOFs = (eof.Beta{i,j}.R2 >= R2_threshold);
            if not(sum(valid_EOFs)==0)
                SR = SR + eof.E{i,j}(:,valid_EOFs)*c(valid_EOFs); % compute the superresolved image
            end
        else
            max_components = min(length(c),num_components);               % max number of EOFs
            SR = SR + eof.E{i,j}(:,1:max_components)*c(1:max_components); % compute the superresolved image
        end
        
        % put SR prediction into output image
        sr_output(hr_ind_vec) = SR;
    end
end
sr_output = sr_output(:);
end
