function [See, frequency, dtimespec, k, bottomdepthavg, Spp] = presdata_to_See(dtime, pressure, habobs, bottomdepth, windowfft, windowavgfft, dtimelims)
%% [See, frequency, dtimespec, k, bottomdepthavg, Spp] = PRESDATA_TO_SEE(dtime, pressure, habobs, bottomdepth, windowfft, windowavgfft, dtimelims)
%
%   inputs
%       - dtime:
%       - pressure:
%       - habobs:
%       - bottomdepth
%       - windowfft:
%       - windowavgfft
%       - dtimelims (optional):
%
%   outputs
%       - See:
%       - frequency:
%       - dtimespec:
%       - k:
%       - bottomdepthavg:
%       - Spp:
%
%
% PRESDATA_TO_SEE.m is a high-level function to compute spectrum of
% sea surface (vertical) elevation (See) from a timeseries of pressure,
% given at a constant height above the bottom habobs, where the
% timeseries of bottom depth is given by bottomdepth.
%
% Spectra are calculated using the parameters windowfft and
% windowavgfft. 
%
% Olavo Badaro Marques.


%% First compute pressure spectra

%
[Spp, dtimespec, frequency, dof, avgpres] = ...
            spectra_scalar_reg(dtime, pressure, ...
                               windowfft, windowavgfft, dtimelims);

Spp = Spp.';


%% Compute average bottom depth on dtimespec grid

% using inefficient code...

%
bottomdepthavg = NaN(length(dtimespec), 1);

%
for i = 1:length(dtimespec)
    %
    linlims_aux = (dtime >= (dtimespec(i) - seconds(windowavgfft/2))) & ...
                  (dtime < (dtimespec(i) + seconds(windowavgfft/2)));
    %
    bottomdepthavg(i) = mean(bottomdepth(linlims_aux), 'omitnan');

end


%% Compute wavenumber

%
freq_bounds = [(1/5000), (1/0.5)];

%
k = wave_freqtok(frequency, bottomdepthavg, [], freq_bounds);


%% Compute transfer function and elevation spectra
    
%
Tspec_p2eta = wave_spec_transferfcn('ptoeta', habobs, k, bottomdepthavg);


%%

See = Spp .* Tspec_p2eta;
