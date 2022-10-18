function [Hsig, freq_mean, freq_peak] = bulkstats_from_wave_spectrum(freq, wvspec, freqband)
%%  [Hsig, freq_mean, freq_peak] = BULKSTATS_FROM_WAVE_SPECTRUM(freq, wvspec, freqband)
%
%   inputs
%       - freq: in Hz.
%       - wvspec: wave spectrum. If a matrix, 
%       - freqband: 1x2 vector with frequency band (in Hz)
%                   where mean frequency and significant wave height
%                   will be computed.
%
%   outputs
%       - Hsig: 
%       - freq_mean:
%       - freq_peak:
%
%
% Adapted from bulk_statistics_from_dspec_onshore_offshore_V3.m
% by Pat Collins.
%
% Olavo Badaro Marques.


%%

% If vector, make sure it's a column vector
% (because I want time to be in the column dimension)
if isvector(wvspec)
    wvspec = wvspec(:);
end

%
if ~exist('freqband', 'var')
    freqband = [0, +Inf];
end


%%

%
linfreqband = (freq >= freqband(1)) & (freq <= freqband(2));

% Subset input within frequency band
freq = freq(linfreqband);
wvspec = wvspec(linfreqband, :);


%% Significant wave height

%
Hsig = 4.*sqrt(trapz(freq, wvspec));


%%

% Find index of the peak for each column in the wave height spectrum
[~, ind] = max(wvspec, [], 1, 'includenan');

% Frequency of the peak in the spectrum/spectra
freq_peak = freq(ind);

% Mean period and frequency
T_mean = trapz(freq, wvspec) ./ trapz(freq, freq.* wvspec);
freq_mean = 1./T_mean;



