function hr_gap_filled = gap_fill_with_eofs(EOF,i,j,R2_thresh,hr,obs_pix,...
    missing_percent,fill_in_pixels,params)
%GAP_FILL_WITH_EOFS Gap-fill missing pixels in a high-resolution image
%using learned EOFs.

% Identify the valid EOFs using R^2 threshold
valid_EOFs = EOF.Beta{i,j}.R2 >= R2_thresh;
E          = EOF.E{i,j}(:,valid_EOFs);

% Subtract the mean vector from the observed HR image (containing missing
% pixels)
z_hr = hr - EOF.mean_vector{i,j};

% Decompose the observed pixels over the equivalently masked EOFs
if (missing_percent==0.9)
    % Tikhonov regularize coefficients for high-missingness cases
    AtA = E(obs_pix,:)'*E(obs_pix,:);
    if cond(AtA)>30
        ss     = svd(AtA);
        delta  = ss(1)/2;
        coeffs = (AtA + delta*eye(size(AtA))) \ (E(obs_pix,:)' * z_hr(obs_pix));
    else
        coeffs = pinv(E(obs_pix,:)) * z_hr(obs_pix);
    end
else
    coeffs = pinv(E(obs_pix,:)) * z_hr(obs_pix);
end

% Compute prediction
z_hat = E * coeffs;
pred  = z_hat + EOF.mean_vector{i,j};

% Use pixels from reconstruction to fill-in missing pixels (EOF interpolation)
hr(fill_in_pixels)              = pred(fill_in_pixels);
hr(hr<params.MIN_REASONABLE_SM) = params.MIN_REASONABLE_SM;
hr(hr>params.MAX_REASONABLE_SM) = params.MAX_REASONABLE_SM;
hr_gap_filled                   = hr(:);
end
