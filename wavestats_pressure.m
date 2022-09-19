function statsout = wavestats_pressure(timevec, pressuredata, windowfft, windowavg, timespeccenter)
% statsout = WAVESTATS_PRESSURE(timevec, pressuredata, windowfft, windowavg, timespeccenter)
%
%   inputs
%       - timevec:
%       - pressuredata: timegridded pressure data (vector)????
%       - windowfft: in seconds
%       - windowavg: in seconds
%       - timespeccenter:
%
%   outputs
%       -
%
%
%
%
%

%%

%
dt = timevec(2) - timevec(1);    % datenum or datetime????

%
fs = 1/dt;     % sampling frequency
fnyq = fs/2;   % Nyquist frequency


%%

%
nfft = windowfft*fs;   % 6 minute chunks

%
df = fs/(nfft-1);   % This is the frequency resolution
nnyq = nfft/2 + 1;

% Positive frequency vector
fm = (0:(nnyq-1)) * df;


%% Frequency bands for infragravity (IG)
% and sea and swell (SS) bands (in Hz)

%
freq_lims_IG = [0.005, 0.03];
freq_lims_SS = [0.045, 0.3];

%
lin_IGband = (fm > freq_lims_IG(1)) & (fm < freq_lims_IG(2));
lin_SSband = (fm > freq_lims_SS(1)) & (fm < freq_lims_SS(2));


%%

% output variables
mean_depth = nan*ones(nt,1);
varP = nan*ones(nt,1);
varPdetrend = nan*ones(nt,1);
Spp = nan*ones(nt,nnyq);
meanfreqSS = nan*ones(nt,1);
peakfreqSS = nan*ones(nt,1);


%%




%%

ii = find((SOLOD.time_dnum>=timeL2_dnum(i)-half_t_interval) &  (SOLOD.time_dnum<timeL2_dnum(i)+half_t_interval));

    Pinterval = SOLOD.Pwater(ii);          % time series of pressure in interval
    nint = length(Pinterval);
    dummy_time_sec = [1:nint]'*dt;        % dummy time variable
    mean_depth(i) = mean(Pinterval) + doffz;      % mean depth in interval

    Pdemean = Pinterval - mean_depth(i);  % demeaned pressure as m
    varP(i) = (Pdemean'*Pdemean)/nint;       % full variance of P
    [m,b] = lsfit(dummy_time_sec,Pdemean);
    Pdetrend = Pdemean - (m*dummy_time_sec + b);  % linear detrend
    varPdetrend(i) = (Pdetrend'*Pdetrend)/nint;       % variance of detrend P

    clf
    subplot(3,1,1)
    plot(dummy_time_sec,Pinterval),grid
    xlabel('time of the hour (sec)')
    ylabel('water pressure (m)')
    title(sprintf('ID=%s: %s  depth=%.2f (m)',SOLOD.ID,datestr(timeL2_dnum(i),0),mean_depth(i)))

    if mean_depth(i)>0.3,  % only if deep enough, calculate wave stats
        
      Ap = calculate_fft(Pdetrend,nfft,fs);   %The fourier coefficients
      Spp0 = mean( Ap .* conj(Ap) ) / (nfft * df);  % need to normalize by freq resolution in hz
      Spp0 = 2*Spp0(1:nnyq);  % keep only to the nyquist frequency
      tmp = sum(Spp0)*df;
      Spp00 = Spp0 * varPdetrend(i)/tmp;   %normalize to keep detrended variance
                                           %    [sum(Spp00)  sum(Spp0)]
      Spp(i,:) = Spp0;
    

      doffz = 0.02;
      grav = 9.81;
      k = get_wavenum(2*pi*fm, mean_depth(i));
      correction = cosh( k*mean_depth(i))./cosh(k*doffz);
      See0 = Spp0 .* correction.^2;
      See(i,:) = See0;
    
      HsigSS(i) = 4*sqrt( sum(See0(i_ss))*df);
      HsigIG(i) = 4*sqrt( sum(See0(i_ig))*df);

      [mx,imx] = max(Spp0(i_ss));
      peakfreqSS(i) = fm(i_ss(imx));
      meanfreqSS(i) = sum( fm(i_ss).*See0(i_ss))/sum(See0(i_ss));
      DOF = 2*(nint/nfft)*2 ;  % last *2 is for the 50% overlap with the Hanning windowed data


