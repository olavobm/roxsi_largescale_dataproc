function statsout = wavestats_pressure(timevec, pressuredata, windowfft, windowavg, timespeclims, dt_step)
% statsout = WAVESTATS_PRESSURE(timevec, pressuredata, windowfft, windowavg, timespeclims, dt_step)
%
%   inputs
%       - timevec: time grid of the pressure data (the time grid CAN NOT
%                  be defined .....????)
%       - pressuredata: timegridded pressure data (vector)????
%       - windowfft: in seconds
%       - windowavg: in seconds
%       - timespeclims:
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

% Number of points 
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

%
if ~exist('dt_step', 'var')
    %
    dt_step = windowavg/86400;    % window in days/datenum
    %
    lindependentintervals = true;
else
    lindependentintervals = false;    % not necessarily
end

%
if ~exist('timespeclims', 'var') || isempty(timespeclims)
    
    %
    tstart_dnum = timevec(1) + dt_step;
    tend_dnum = timevec(end) - dt_step;
    %
    timespeclims = [tstart_dnum, tend_dnum];

end

%
timespecintervals = timespeclims(1): dt_step : timespeclims(2);

%
nt = length(timespecintervals);
half_t_interval = dt_step/2;


%%

% output variables
mean_depth = nan*ones(nt,1);
varP = nan*ones(nt,1);
varPdetrend = nan*ones(nt,1);
Spp = nan*ones(nt,nnyq);

%
meanfreqSS = nan*ones(nt,1);
peakfreqSS = nan*ones(nt,1);


%% Find the number of points per averaging window

%
npts_ininterval = fs*3600*(timespecintervals(2) - timespecintervals(1));
%
npts_ininterval = npts_ininterval - 1;


%%

if lindependentintervals

    %
    time_bounds = [(timespecintervals - half_t_interval); ...
                   (timespecintervals + half_t_interval)];

    %
    ind_start_data = find(~isnan(pressuredata), 1, 'first');
    ind_end_data = find(~isnan(pressuredata), 1, 'last');


    %
    ind_first_wholeinterval = find(time_bounds(1, :) >= timevec(ind_start_data), 1, 'first');

    %
    ind_data_start_firstinterval = find(timevec==time_bounds(1, ind_first_wholeinterval));

    %
    nt_intervals = floor(length(ind_data_start_firstinterval : ind_end_data) / npts_ininterval);

    %
    ind_data_wholeintervals = ind_data_start_firstinterval + ...
                                   (0 : 1 : (nt_intervals*npts_ininterval));
    %
    pressure_perinterval = reshape(pressuredata(ind_data_wholeintervals), npts_ininterval, nt_intervals);

    %
    

    %
    % ------------------------------------------------



%     %
%     ind_start_grid = find(time_bounds(1, :) >= timevec(ind_start_data), 1, 'first');
%     ind_end_grid = find(time_bounds(2, :) <= timevec(ind_end_data), 1, 'last');
    %
    ind_start_grid = find(time_bounds(1, 2:end-1) == timevec(ind_start_data));
    ind_end_grid = find(time_bounds(2, 2:end-1) == timevec(ind_end_data));
    %
    ind_in_grid = ind_start_grid:ind_end_grid;

    % not right!!! inds_in_grid is not in the data!!! THIS MUST/SHOULD BE AN INTEGER!!!
    n_intervalsindata = length(ind_in_grid) / npts_ininterval;

    %
    reshape(pressuredata(ind_in_grid), npts_ininterval, n_intervalsindata);
    %
    

else
    error('Complicated option not implemented')

end

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


