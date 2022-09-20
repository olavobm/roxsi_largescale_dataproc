function bulk_statistics_from_dspec()
%%
%
%
%
%
%
%
%
%
% Adapted from bulk_statistics_from_dspec_onshore_offshore_V3.m
% by Pat Collins.
%
% Olavo Badaro Marques.


    % preallocate for speed
    dir_peak = nan(size(dspec.dtime));
    spread_peak = nan(size(dspec.dtime));
    theta0 = nan(length(dspec.dtime), length(dspec.f));
    sigma0 = nan(length(dspec.dtime), length(dspec.f));
    ind_f_pk = nan(size(dspec.dtime));
    
    %reshape theta and convert to radians
    theta = deg2rad(theta);
    theta = theta(:);
    guess = deg2rad(guess);

    %generate theta, the variable used for integration etc, note here that
    %we will need to loop around from 0 to 360 (the vale from 0 to 360 so
    %that we can to the integration correctly 


    %%

[~, ind] = max(dspec.S_f_IMLM,[],'omitnan');
f_peak = dspec.f(ind);
T_peak = 1./f_peak;
T_mean = trapz(dspec.f, dspec.S_f_IMLM) ./ trapz(dspec.f, dspec.f .* dspec.S_f_IMLM);
h_sig = 4.*sqrt(trapz(dspec.f, dspec.S_f_IMLM));

dir_peak = nan(size(dspec.dtime));
spread_peak = nan(size(dspec.dtime));
theta0 = nan(length(dspec.dtime), length(dspec.f));
sigma0 = nan(length(dspec.dtime), length(dspec.f));
ind_f_pk = nan(size(dspec.dtime));

for timestep = 1:length(dspec.dtime)

% the index of the peak frequency
[~, ind_f_pk(timestep)] = max(dspec.S_f_IMLM(:,timestep));

for ii = 1 : length(dspec.f)

%create a circular version of the directional spread 
D = dspec.D_f_theta_IMLM(ii, :, timestep);
D = D(:);

% use fzero to solve for theta0, per timestep, per frequency.
theta0(timestep,ii) = fzero(@(theta0) trapz(theta, sin(theta-theta0).*D) , guess);
sigma0(timestep,ii) = rad2deg(sqrt(4.*trapz(theta, (sin((theta-theta0(timestep,ii))./2).^2).*D)));

end
end

theta0 = mod(theta0, 2*pi);

% this weight is for the directional averaging because it is a circular
% variable, the normal energy weighting has to be modified. 
% Lines have been added here to remove the impact when the spread
% function is uniform
weight = (dspec.f(2)-dspec.f(1)).* dspec.S_f_IMLM ./ trapz(dspec.f, dspec.S_f_IMLM);
    weight(isnan(theta0')) = 0;
dir_mean = mod(atan2d(sum(weight.*sind(rad2deg(theta0(:,:)'))),sum(weight.*cosd(rad2deg(theta0(:,:)')))),360);
spread_mean = trapz(dspec.f, dspec.S_f_IMLM.*sigma0') ./ trapz(dspec.f, dspec.S_f_IMLM);

% these lines pick out the direction and spread for the peak frequency
for ii_dt = 1 : length(dspec.dtime)

dir_peak(ii_dt) = rad2deg(theta0(ii_dt,ind_f_pk(ii_dt)));
spread_peak(ii_dt) = sigma0(ii_dt, ind_f_pk(ii_dt));

end

IMLM.S_f = dspec.S_f_IMLM;
IMLM.f_peak = f_peak;
IMLM.T_peak = T_peak;
IMLM.T_mean = T_mean;
IMLM.h_sig = h_sig;

IMLM.dir_mean_f = rad2deg(theta0)';
IMLM.spread_mean_f = sigma0';
IMLM.dir_mean = dir_mean(:);
IMLM.spread_mean = spread_mean(:);
IMLM.dir_peak = dir_peak(:);
IMLM.spread_peak = spread_peak(:);

bulk.IMLM = IMLM;