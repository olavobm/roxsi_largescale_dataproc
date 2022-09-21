function bulkoutput = bulkstats_from_wave_dirspectrum(freq, dirvec, Dir_f_fcn, wvdirspec, wvspec)
%% bulkoutput = BULKSTATS_FROM_WAVE_DIRSPECTRUM(freq, dirvec, Dir_f_fcn, wvdirspec, wvspec)
%
%   inputs
%       - freq: in Hz.
%       - dirvec:
%       - Dir_f_fcn:
%       - wvdirspec:
%       - wvspec:
%
%   outputs
%       - bulkoutput: 
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
dirinradians = deg2rad(dirvec);
dirinradians = dirinradians(:);

% % % May want to include this later
% % guess = deg2rad(guess);

% % % I should be able to compute this from
% % % integrating the directional spectrum,
% % % and all methods give the same result, right?!
% % % If so, then I don't need to input wvspec.
% % % 
% % waveSpec = wvspec.S_f


%%

% % %
% % if size(wvdirspec, 3)~=1
% %     wvdirspec
% % end

%
Nfrequency = length(freq);
Ntime = size(wvdirspec, 3);


%%

%
prealloc_aux = NaN(1, Ntime);
%
dir_peak = prealloc_aux;
spread_peak = prealloc_aux;
%
ind_f_pk = prealloc_aux;
%
theta0 = NaN(Nfrequency, Ntime);
sigma0 = NaN(Nfrequency, Ntime);

%
disp(['-------------- Total number of time stamps is = ' num2str(Ntime) ' --------------'])


%% Compute frequency spectra of the mean direction and mean direction spread

% Loop over time
%
% (Skip the last one because I need to fix the code that computes spectrum)
tic
for i1 = 1:(Ntime-1)

    % Loop over frequency
    for i2 = 1 : Nfrequency
    
        % Create a circular version of the directional spread 
        D = Dir_f_fcn(i2, :, i1);
        D = D(:);
        %
        [~, ind_max] = max(D);
        guess = dirinradians(ind_max);
        
        % use fzero to solve for theta0, per timestep, per frequency.
        theta0(i2, i1) = fzero(@(theta0find) trapz(dirinradians, sin(dirinradians-theta0find).*D), guess);
        sigma0(i2, i1) = rad2deg(sqrt(4.*trapz(dirinradians, (sin((dirinradians-theta0(i2, i1))./2).^2).*D)));
    end

    % Print progress message
    fractional_progress = round(1000*i1/Ntime)/10;
    %
    if fractional_progress~=0 && mod(fractional_progress, 10)==0
        disp(['Done with timestamp ' num2str(i1) ' out of ' num2str(Ntime)])
        toc
    end
end


%
theta0 = mod(theta0, 2*pi);


%% Compute mean direction, mean directional spread, peak direction
% and directional spread of the peak direction
%
% this weight is for the directional averaging because it is a circular
% variable, the normal energy weighting has to be modified.
%
% Lines have been added here to remove the impact when the spread
% function is uniform (???????????)

%
weight = (freq(2) - freq(1)) .* wvspec ./ trapz(freq, wvspec);
weight(isnan(theta0)) = 0;

% Sum over frequencies
dir_mean = mod(atan2d(sum(weight.*sin(theta0)), ...
                      sum(weight.*cos(theta0))), ...
               360);
%
spread_mean = trapz(freq, wvspec.*sigma0) ./ trapz(freq, wvspec);

% these lines pick out the direction and spread for the peak frequency
for ii_dt = 1:(Ntime - 1)

    % the index of the peak frequency
    [~, ind_f_pk(ii_dt)] = max(wvspec(:, ii_dt));

    %
    dir_peak(ii_dt) = rad2deg(theta0(ind_f_pk(ii_dt), ii_dt));
    spread_peak(ii_dt) = sigma0(ind_f_pk(ii_dt), ii_dt);
end


%% Assign bulk statistics to output variable

%
bulkoutput.dir_mean_f = rad2deg(theta0);
bulkoutput.spread_mean_f = sigma0;

%
bulkoutput.dir_mean = dir_mean;
bulkoutput.spread_mean = spread_mean;
%
bulkoutput.dir_peak = dir_peak;
bulkoutput.spread_peak = spread_peak;

