function [meanfreq, peakfreq] = wave_bulk_frequency(frequency, See, freqlims, fdim)
%%  [meanfreq, peakfreq] = WAVE_BULK_FREQUENCY(frequency, See, freqlims, fdim)
%
%   inputs
%       - frequency: frequency vector.
%       - See: surface vertical elevation spectrum(a). It can be a vector
%              or a 2D matrix, but not a higher dimensional array.
%       - freqlims: frequency limits for computing bulk quantities.
%       - fdim: the dimension of frequency in the See array.
%
%   outputs
%       - meanfreq: mean frequency.
%       - peakfreq: peak frequency.
%
%
% WAVE_BULK_FREQUENCY.m computes mean frequency and peak frequency from
% surface elevation spectrum(a). This function assumes that dimensions
% of frequency and See are consistent (e.g. frequency in Hz and See in
% meters squared per Hz).


%%

if ~exist('fdim', 'var')
    fdim = 1;
end


%%

if ~exist('freqlims', 'var') || isempty(freqlims)
    freqlims = [-Inf, +Inf];
end

%% 

%
frequency = frequency(:);

%
if isvector(See)
    if isrow(See)
        ltranspose = true;
    else
        ltranspose = false;
    end
%
else
    %
    if fdim==2
        ltranspose = true;
    else
        ltranspose = false;
    end
end

%
if ltranspose
    See = See.';
end


%% Compute mean frequency

%
linfreqlims = (frequency >= freqlims(1)) & (frequency <= freqlims(2));

%
freq_array = repmat(frequency(linfreqlims), 1, size(See, 2));

% Zero'th and first moments of the surface elevation spectrum(a)
m0 = trapz(frequency(linfreqlims), See(linfreqlims, :), 1);
m1 = trapz(frequency(linfreqlims), See(linfreqlims, :) .* freq_array, 1);

% % %
% % Hsig = 4.*sqrt(m0);

%
meanfreq = m1./m0;


%% Compute peak frequency

%
[~, ind_peakfreq] = max(See, [], 1);

%
peakfreq = frequency(ind_peakfreq);


%% Transpose back if necessary

%
if ~isvector(See)
    if ltranspose
        meanfreq = meanfreq.';
% %         peakfreq = peakfreq.';
    end
end


