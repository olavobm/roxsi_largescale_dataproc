function bulkstats = bulkstats_from_dirspec(frequency, theta, See_theta, D_theta)
%% bulkstats = BULKSTATS_FROM_DIRSPEC(frequency, theta, See_theta, D_theta)
%
%   inputs
%       - theta: in degrees.
%       - See_theta:
%       - D_theta:
%
%   outputs
%       - bulkstats: structure with fields containing bulk statistics.
%
%
% BULKSTATS_FROM_DIRSPEC.m
%
%
%
%


%% Compute elevation spectra by integrating the directional spectra
%
% (HOW DO I INTEGRATE THIS PROPERLY???)

%
% % See = trapz(theta, See_theta, ???);    % get right dimensions!!!

%%

[meandir_f, meandirspread_f] = wavedirection_from_dirspectra(theta, See_theta, D_theta);


%%

%
df = frequency(2) - frequency(1);

% This weight is for the directional averaging because it is a circular
% variable, the normal energy weighting has to be modified.
weight = (df.* See) ./ trapz(frequency, See);

% Lines have been added here to remove the impact when the spread
% function is uniform (????)
weight(isnan(theta0')) = 0;
theta0(isnan(theta0)) = NaN;  % ??? only the EMEM had this pat's code
    
    
%%
dir_mean = mod(atan2d(sum(weight.*sind(rad2deg(theta0(:,:)'))),sum(weight.*cosd(rad2deg(theta0(:,:)')))),360);
    sigma0(isnan(sigma0)) = nan;    % ??? only the EMEM had this pat's code
spread_mean = trapz(dspec.f, dspec.S_f.*sigma0') ./ trapz(dspec.f, dspec.S_f);



%%
return

%%



EMEM.dir_mean_f = rad2deg(theta0)';
EMEM.spread_mean_f = sigma0';
EMEM.dir_mean = dir_mean(:);
EMEM.spread_mean = spread_mean(:);
EMEM.dir_peak = dir_peak(:);
EMEM.spread_peak = spread_peak(:);


%%

%
dir_peak = NaN(size(dspec.dtime));
spread_peak = NaN(size(dspec.dtime));
%
theta0 = NaN(length(dspec.dtime), length(dspec.f));
sigma0 = NaN(length(dspec.dtime), length(dspec.f));
%
ind_f_pk = nan(size(dspec.dtime));



%%



%%
% these lines pick out the direction and spread for the peak frequency
for ii_dt = 1 : length(dspec.dtime)
    dir_peak(ii_dt) = rad2deg(theta0(ii_dt,ind_f_pk(ii_dt)));
    spread_peak(ii_dt) = sigma0(ii_dt, ind_f_pk(ii_dt));
end


%% Organize output

EMEM.dir_mean_f = rad2deg(theta0)';
EMEM.spread_mean_f = sigma0';
EMEM.dir_mean = dir_mean(:);
EMEM.spread_mean = spread_mean(:);
EMEM.dir_peak = dir_peak(:);
EMEM.spread_peak = spread_peak(:);

bulk.EMEM = EMEM;


        
%%
return


%%

dir_mean_f = rad2deg(theta0)';
spread_mean_f = sigma0';

dir_mean = dir_mean(:);
spread_mean = spread_mean(:);

dir_peak = dir_peak(:);
spread_peak = spread_peak(:);




    %generate theta, the variable used for integration etc, note here that
    %we will need to loop around from 0 to 360 (the vale from 0 to 360 so
    %that we can to the integration correctly 