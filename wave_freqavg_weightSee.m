function [xavg] = wave_freqavg_weightSee(frequency, See, xspec, freqlims)
%% [xavg] = WAVE_FREQAVG_WEIGHTSEE(frequency, See, xspec, freqlims)
%
%   inputs
%       - frequency:
%       - See:
%       - xspec:
%       - freqlims:
%
%   outputs
%       - xavg:
%
%
% WAVE_FREQAVG_WEIGHTSEE.m
%
%
%


%%

linfreqlims = (frequency >= freqlims(1)) & (frequency <= freqlims(2));


%%

%
m0 = trapz(frequency(linfreqlims), See(linfreqlims));
% % m1 = trapz(f(linfreqlims), See(linfreqlims) .* f(linfreqlims));

%
xavg = 1./m0 .* trapz(frequency(linfreqlims), xspec(linfreqlims).*See(linfreqlims));
