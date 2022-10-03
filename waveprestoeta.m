function waveeta = waveprestoeta(wavep, hi, H, dtsampling, freqcutoff)
%% waveeta = WAVEPRESTOETA(wavep, hi, H, dtsampling, freqcutoff)
%
%   inputs
%       - wavep: wave pressure anomaly (in dbar).
%       - hi: height above the bottom (in meters) of the
%             instrument/pressure data.
%       - H: bottom depth (in meters).
%       - dtsampling: sampling period (in seconds).
%       - freqcutoff: high-frequency cut off (in Hertz).
%
%   outputs
%       - waveeta: surface displacement (in meters) relative
%                  to the mean position.
%
%
% WAVEPRESTOETA.m function to 
%
% NOTE THAT THIS IS VERY INEFFICIENT IF NUMBER OF DATA
% IS GREATER THAN ?????.
%
%
% Olavo Badaro Marques.



%%

%
rho0 = 1030;
g = 9.8;

%%

%
Npts = length(wavep);

% Frequency resolution
df = 1/(Npts*dtsampling);


%%

% Frequency vector, in Hz (does NOT start at 0)
frequency = (1:floor(Npts/2)) .* df;

%
labovecutoff = (frequency > freqcutoff);


%% Compute Fourier coefficients

pressure_fft_coefs = fft(wavep);


%% Function handle to compute wavenumbers

% (PS: frequency in Hz and k in radians per meter)
disp_rel = @(k, freq) g*k*tanh(k*H) - (2*pi*freq)^2;

%
kvec = NaN(size(frequency));

%
for i = 1:length(frequency)
    %
    disp_rel_atfreq0 = @(k) disp_rel(k, frequency(i));

    %
    kvec(i) = fzero(disp_rel_atfreq0, [(2*pi/(5000)), (2*pi/(1))]);

end


%% Create the wavenumber vector that matches the fft coefficients

%
kvec_NaNabovecutoff = kvec;
kvec_NaNabovecutoff(labovecutoff) = NaN;

% Indice to take into account both even or odd Npts
ind_topedge = floor((Npts-1)/2);

% NaN for mean, then positive frequencies, then negative frequencies
kvec_forfft = [NaN, kvec_NaNabovecutoff(1:ind_topedge), ...
                    fliplr(kvec_NaNabovecutoff)];


%% Compute surface displacement Fourier coefficients

%
pres2ssh = cosh(kvec_forfft .* H) .* cosh(kvec_forfft .* hi);
pres2ssh(isnan(pres2ssh)) = 0;

%
eta_fft_coefs = pres2ssh(:) .* pressure_fft_coefs(:);


%% Do inverse fft to compute timeseries of surface elevation

waveeta = real(ifft(eta_fft_coefs));






