function wvspec = bulk_statistics_from_dspec(wvspec, dirvec, dspec)
%% bulk = BULK_STATISTICS_FROM_DSPEC(wvspec, dirvec, dspec)
%
%   inputs
%       - dspec:
%       - theta: the direction vector associated with the directional
%                spectrum. In degrees, clockwise from true north, from
%                zero towards 360).
%
%   outputs
%       - bulk: 
%
%
% Turn dspec into bulk statistics. Check if the computer has
% parallel computing.
%
% Direction guess (guess variable should be an option for
% offshore propagating waves, because mean direction is not
% good for it).
%
% Adapted from bulk_statistics_from_dspec_onshore_offshore_V3.m
% by Pat Collins.
%
% Olavo Badaro Marques.

% Computed statistiscal quantities are:
%   - preak frequency:
%   - signifcant wave height:
%
%   - dir_mean_f:
%   - spread_mean_f:
%   - dir_mean:
%   - spread_mean:
%
%   - dir_peak:
%   - spread_peak:


%%

% Compute direction in radians and turn into column
dirtheta = deg2rad(dirvec);
dirtheta = dirtheta(:);

% % %
% % guess = deg2rad(guess);


%% Compute sta

% Find index of the peak in the power spectrum
[~, ind] = max(wvspec.S_f, [], 'omitnan');

% Frequency of the peak in the spectrum
f_peak = wvspec.f(ind);
% % T_peak = 1./f_peak;

% Mean period
T_mean = trapz(wvspec.f, wvspec.S_f) ./ trapz(wvspec.f, wvspec.f .* wvspec.S_f);

% Significant wave height (from the full Spectrum)
h_sig = 4.*sqrt(trapz(wvspec.f, wvspec.S_f));


%%

%
Ntime = length(wvspec.dtime);
Nfrequency = length(wvspec.f);

%
vec_like_dtime = NaN(size(wvspec.dtime));
%
dir_peak = vec_like_dtime;
spread_peak = vec_like_dtime;
%
theta0 = NaN(Ntime, Nfrequency);
sigma0 = NaN(Ntime, Nfrequency);
%
ind_f_pk = vec_like_dtime;

%
disp(['Total number of timestamps is = ' num2str(length(wvspec.dtime))])

% Loop over time (Skip the last one because I need to fix the code
% that computes spectrum)
tic
for timestep = 1:(length(wvspec.dtime)-1)

    % the index of the peak frequency
    [~, ind_f_pk(timestep)] = max(wvspec.S_f(:, timestep));

    % Loop over frequency
    for ii = 1 : length(wvspec.f)
    
        % Create a circular version of the directional spread 
        D = dspec.D_f_theta(ii, :, timestep);
        D = D(:);
        %
        [~, ind_max] = max(D);
        guess = dirtheta(ind_max);
        
        % use fzero to solve for theta0, per timestep, per frequency.
        theta0(timestep, ii) = fzero(@(theta0) trapz(dirtheta, sin(dirtheta-theta0).*D), guess);
        sigma0(timestep, ii) = rad2deg(sqrt(4.*trapz(dirtheta, (sin((dirtheta-theta0(timestep, ii))./2).^2).*D)));
    end

    % Print progress message
    fractional_progress = round(100*timestep/length(wvspec.dtime));
    %
    if fractional_progress~=0 && mod(fractional_progress, 10)==0
        disp(['Done with timestamp ' num2str(timestep) ' out of ' num2str(length(wvspec.dtime))])
        toc
    end
end


%
theta0 = mod(theta0, 2*pi);

%%
% this weight is for the directional averaging because it is a circular
% variable, the normal energy weighting has to be modified. 
% Lines have been added here to remove the impact when the spread
% function is uniform
weight = (wvspec.f(2) - wvspec.f(1)) .* wvspec.S_f ./ trapz(wvspec.f, wvspec.S_f);
weight(isnan(theta0')) = 0;

% THIS NEEDS TO BE BROKEN DOWN FOR CLARITY!!!!
dir_mean = mod(atan2d(sum(weight.*sind(rad2deg(theta0(:,:)'))), sum(weight.*cosd(rad2deg(theta0(:,:)')))), 360);
spread_mean = trapz(wvspec.f, wvspec.S_f.*sigma0') ./ trapz(wvspec.f, wvspec.S_f);

% these lines pick out the direction and spread for the peak frequency
for ii_dt = 1 : (length(wvspec.dtime) - 1)
    %
    dir_peak(ii_dt) = rad2deg(theta0(ii_dt, ind_f_pk(ii_dt)));
    spread_peak(ii_dt) = sigma0(ii_dt, ind_f_pk(ii_dt));
end


%% Assign bulk statistics to output variable

% Bulk statistic from elevation spectrum
% % IMLM.S_f = dspec.S_f;
wvspec.f_peak = f_peak;
% % IMLM.T_peak = T_peak;
wvspec.T_mean = T_mean;
wvspec.h_sig = h_sig;

% Bulk statistic from directional wave spectrum
dspec.dir_mean_f = rad2deg(theta0)';
dspec.spread_mean_f = sigma0';
dspec.dir_mean = dir_mean(:);
dspec.spread_mean = spread_mean(:);
dspec.dir_peak = dir_peak(:);
dspec.spread_peak = spread_peak(:);


% % 
% % bulk.IMLM = IMLM;