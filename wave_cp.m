function cp = wave_cp(k, H, g)
%% cg = WAVE_CP(k, H, g)
%
%   inputs
%       - k: wavenumber (in radians per meter).
%       - H: bottom depth (in meters, positive).
%       - g (optional): acceleration of gravity.
%
%   outputs
%       - cp: phase speed (in m/s).
%
%
% WAVE_CG.m computes phase speed of surface gravity waves
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
cp = sqrt(g * tanh(k.*H) ./ k);