function [f, a1, a2, b1, b2, Ezz] = wave_dirmoments_xyz(x, y, z, dt, nfft, noverlap, window)
%% [f, a1, a2, b1, b2, Ezz] = WAVE_DIRMOMENTS_XYZ(x, y, z, dt, nfft, noverlap, window)
%
%   inputs
%       - x, y, z:
%       - dt:
%       - nfft:
%       - overlap:
%       - window:
%
%   outputs
%       - f: frequency vector (in Hz).
%       - a1, a2, b1, b2:
%       - Ezz:
%       - Hsig:
%       - T_mean, dir_mean, spread_mean:
%       - T_peak, dir_peak, spread_peak:
%       - a1_bar, b1_bar:
%       - DOF: degrees of freedom.
%
%
% WAVE_DIRMOMENTS_UVP.m takes displacement data and calculates
% the wave spectra using the fourier method. This approach is coded in a
% similar way to the sofar script, their documentation is useful here.
%
% Based on PC Apr 2022. Adapted by Olavo Marques Oct/2022.

%%

% % if ~exist('freq_lims', 'var')
% %     freq_lims = [0, Inf];
% % end


%%

%
x = detrend(x);
y = detrend(y);
z = detrend(z);


%%
 
%
[Ezz, f] = cpsd(z, z, window, noverlap, nfft, 1./dt);
Exx = cpsd(x, x, window, noverlap, nfft, 1./dt);
Eyy = cpsd(y, y, window, noverlap, nfft, 1./dt);
Exz = cpsd(x, z, window, noverlap, nfft, 1./dt);
Eyz = cpsd(y, z, window, noverlap, nfft, 1./dt);
Exy = cpsd(x, y, window, noverlap, nfft, 1./dt);


%% From wave buoy displacements

%
a1 = imag(Exz) ./ (sqrt(Ezz.*(Exx + Eyy)));
b1 = imag(Eyz) ./ (sqrt(Ezz.*(Exx + Eyy)));
a2 = (Exx - Eyy)./(Exx + Eyy);
b2 = (2.*real(Exy)) ./ (Exx + Eyy);




%%
return

%%

%
numwin = floor(length(x)./nfft) + floor(length(x)./nfft) - 1;
% DOF = 8/3 * 2 * numwin;    % what is this 8/3 ??? numwin looks like the
%                            % correct number of chunks over which fft was
%                            % computed for
DOF = 2 * numwin;
