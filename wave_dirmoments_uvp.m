function [f, a1, a2, b1, b2, Ezz] = wave_dirmoments_uvp(u, v, p, dt, nfft, noverlap, window)
%% [f, a1, a2, b1, b2, Ezz] = WAVE_DIRMOMENTS_UVP(u, v, p, zhabUV, zhabP, dt, nfft, noverlap, window)
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
x = detrend(u);
y = detrend(v);
z = detrend(p);


%%
 
%
[Ezz, f] = cpsd(z, z, window, noverlap, nfft, 1./dt);
Exx = cpsd(x, x, window, noverlap, nfft, 1./dt);
Eyy = cpsd(y, y, window, noverlap, nfft, 1./dt);
Exz = cpsd(x, z, window, noverlap, nfft, 1./dt);
Eyz = cpsd(y, z, window, noverlap, nfft, 1./dt);
Exy = cpsd(x, y, window, noverlap, nfft, 1./dt);


% % %% For Spotter
% % 
% % %
% % a1 = imag(Exz) ./ (sqrt(Ezz.*(Exx + Eyy)));
% % b1 = imag(Eyz) ./ (sqrt(Ezz.*(Exx + Eyy)));
% % a2 = (Exx - Eyy)./(Exx + Eyy);
% % b2 = (2.*real(Exy)) ./ (Exx + Eyy);


%% From pressure, ADCP (u, v)

%
a1 = real(Exz) ./ (sqrt(Ezz.*(Exx + Eyy)));
b1 = real(Eyz) ./ (sqrt(Ezz.*(Exx + Eyy)));
a2 = (Exx - Eyy)./(Exx + Eyy);
b2 = (2.*real(Exy)) ./ (Exx + Eyy);

%%

return

%%

linfreqlims = (f >= freq_lims(1)) & (f <= freq_lims(2));

%% Compute bulk statistics

%
m0 = trapz(f(linfreqlims), Ezz(linfreqlims));
m1 = trapz(f(linfreqlims), Ezz(linfreqlims) .* f(linfreqlims));
%
a1_bar = 1./m0 .* trapz(f(linfreqlims), a1(linfreqlims).*Ezz(linfreqlims));
b1_bar = 1./m0 .* trapz(f(linfreqlims), b1(linfreqlims).*Ezz(linfreqlims));

% % % Same as above, but using the full spectrum
% % m0 = trapz(f, Ezz);
% % m1 = trapz(f, Ezz .* f);
% % a1_bar = 1./m0 .* trapz(f, a1.*Ezz);
% % b1_bar = 1./m0 .* trapz(f, b1.*Ezz);

% % %
% % Hsig = 4.*sqrt(m0);
% % 
% % %
% % T_mean = m0./m1;

%
dir_mean = mod(90 - atan2d(b1_bar, a1_bar), 360);
%
spread_mean = rad2deg(sqrt(2.*(1 - sqrt(a1_bar.^2 + b1_bar.^2))));


%%
return


%% Compute bulk statistics -- peak quantities

%
[~, f_peak_ind] = max(Ezz);

%
f_peak = f(f_peak_ind);
T_peak = 1./f_peak;

%
dir_peak = mod(90 - atan2d(b1(f_peak_ind), a1(f_peak_ind)),360);
spread_peak = rad2deg(sqrt(2.*(1 - sqrt(a1(f_peak_ind).^2 + b1(f_peak_ind).^2))));


%%

%
numwin = floor(length(x)./nfft) + floor(length(x)./nfft) - 1;
% DOF = 8/3 * 2 * numwin;    % what is this 8/3 ??? numwin looks like the
%                            % correct number of chunks over which fft was
%                            % computed for
DOF = 2 * numwin;
