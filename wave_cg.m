function [cg, kH] = wave_cg(k, H, g)
%% [cg, kH] = WAVE_CG(k, H, g)
%
%   inputs
%       - k: wavenumber (in radians per meter).
%       - H: bottom depth (in meters, positive).
%       - g (optional): acceleration of gravity.
%
%   outputs
%       - cg: group velocity (in m/s).
%
%
% WAVE_CG.m computes group velocity of surface gravity waves
% using linear wave theory.
%
%
% Olavo Badaro Marques, 18/01/2023.


%%

if ~exist('g', 'var')
    g = 9.8;
end


%%

%
kH = k .* H;

%
cp = wave_cp(k, H, g);

%
cg = cp .* 0.5 .* (1 + (2*kH./sinh(2*kH)));

