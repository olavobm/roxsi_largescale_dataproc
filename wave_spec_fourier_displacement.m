function [a1, a2, b1, b2, ...
          Ezz, Hsig, T_mean, ...
          dir_mean, spread_mean, f_peak, T_peak, dir_peak, spread_peak, ...
          f, DoF] = wave_spec_fourier_displacement(x, y, z, dt, nfft, noverlap, window)
%% function [a1, a2, b1, b2, ...
%           Ezz, Hsig, T_mean, ...
%           dir_mean, spread_mean, f_peak, ...
%           T_peak, dir_peak, spread_peak, f, DoF] = ...
%                          wave_spec_fourier_displacement(x, y, z, dt, nfft, noverlap, window)
%
%
% WAVE_SPEC_FOURIER_DISPLACEMENT.m takes displacement data and calculates
% the wave spectra using the fourier method. This approach is coded in a
% similar way to the sofar script, their documentation is useful here.
%
% PC Apr 2022

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

%%
%
a1 = imag(Exz) ./ (sqrt(Ezz.*(Exx + Eyy)));
b1 = imag(Eyz) ./ (sqrt(Ezz.*(Exx + Eyy)));
a2 = (Exx - Eyy)./(Exx + Eyy);
b2 = (2.*real(Exy)) ./ (Exx + Eyy);


%%
%
m0 = trapz(f, Ezz);
m1 = trapz(f, Ezz .* f);
a1_bar = 1./m0 .* trapz(f, a1.*Ezz);
b1_bar = 1./m0 .* trapz(f, b1.*Ezz);

%
Hsig = 4.*sqrt(m0);

%
T_mean = m0./m1;
dir_mean = mod(90 - atan2d(b1_bar,a1_bar),360);

spread_mean = rad2deg(sqrt(2.*(1 - sqrt(a1_bar.^2 + b1_bar.^2))));
[~, f_peak_ind] = max(Ezz);

%
f_peak = f(f_peak_ind);
T_peak = 1./f_peak;

%
dir_peak = mod(90 - atan2d(b1(f_peak_ind),a1(f_peak_ind)),360);
spread_peak = rad2deg(sqrt(2.*(1 - sqrt(a1(f_peak_ind).^2 + b1(f_peak_ind).^2))));

%%

%
numwin = floor(length(x)./nfft) + floor(length(x)./nfft) - 1;

DoF = 8/3 * 2 * numwin;

