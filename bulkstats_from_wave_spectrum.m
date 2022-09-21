function [freq_peak, freq_mean, Hsig] = bulkstats_from_wave_spectrum(freq, wvspec)
%%  [freq_peak, freq_mean, Hsig] = BULKSTATS_FROM_WAVE_SPECTRUM(freq, wvspec)
%
%   inputs
%       - freq: in Hz.
%       - wvspec: wave spectrum. If a matrix, 
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
if isvector(wvspec)
    wvspec = wvspec(:);
end


%%

% Find index of the peak for each column in the wave height spectrum
[~, ind] = max(wvspec, [], 1, 'includenan');

% Frequency of the peak in the spectrum
freq_peak = freq(ind);

% Mean period and frequency
T_mean = trapz(freq, wvspec) ./ trapz(freq, freq.* wvspec);
freq_mean = 1./T_mean;


%% Significant wave height

%
Hsig = 4.*sqrt(trapz(freq, wvspec));
