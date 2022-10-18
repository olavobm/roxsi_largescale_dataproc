function [freq_peak, freq_mean, Hsig] = bulkstats_from_wave_spectrum(freq, wvspec, freqband)
%%  [freq_peak, freq_mean, Hsig] = BULKSTATS_FROM_WAVE_SPECTRUM(freq, wvspec, freqband)
%
%   inputs
%       - freq: in Hz.
%       - wvspec: wave spectrum. If a matrix, 
%       - freqband: 1x2 vector with frequency band (in Hz)
%                   where mean frequency and significant wave height
%                   will be computed.
%
%   outputs
%       - freq_peak:
%       - freq_mean:
%       - Hsig: 
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


%%

% Find index of the peak for each column in the wave height spectrum
[~, ind] = max(wvspec, [], 1, 'includenan');

% Frequency of the peak in the spectrum/spectra
freq_peak = freq(ind);

% Mean period and frequency
T_mean = trapz(freq, wvspec) ./ trapz(freq, freq.* wvspec);
freq_mean = 1./T_mean;


%% Significant wave height

%
Hsig = 4.*sqrt(trapz(freq, wvspec));
