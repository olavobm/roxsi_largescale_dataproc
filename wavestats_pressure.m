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
% By the way, the advantage of using datetime instead of datenum is
% that time differences can be obtained exactly, instead of some
% including rounding differences that appear with datenum.
%
%

%% Check inputs

%
if ~isdatetime(timevec)
    error('Time grid is NOT in datetime format.')
end

%
if isempty(timespeclims.TimeZone)
    timespeclims.TimeZone = timevec.TimeZone;
end


%%

% Time difference, in seconds, of the data time grid
% (the code assumes it is a equally time grid)
dt = (timevec(2) - timevec(1));
dt = seconds(dt);

%
fs = 1/dt;     % sampling frequency
fnyq = fs/2;   % Nyquist frequency


%%

% Number of points per chunk
nfft = windowfft*fs;

% Only allow for even nfft (you could want otherwise
% is sampling period is 1.5s, but that seems unlikely)
if mod(nfft, 2)~=0
    error('nfft is odd. It must be even.')
end

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


%%

%
if ~exist('dt_step', 'var')
    %
    dt_step = windowavg/(24*3600);    % window in days/datenum
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

% % % output variables
% % mean_depth = nan*ones(nt,1);
% % varP = nan*ones(nt,1);
% % varPdetrend = nan*ones(nt,1);
% % Spp = nan*ones(nt,nnyq);
% % 
% % %
% % meanfreqSS = nan*ones(nt,1);
% % peakfreqSS = nan*ones(nt,1);


%% Find the number of points per averaging window

%
npts_ininterval = fs*windowavg;


%%

%
statsout.dtime = timespecintervals;
statsout.frequency = fm;
%
statsout.avgpres = NaN(1, length(timespecintervals));
statsout.presSpec = NaN(length(fm), length(timespecintervals));


%%

if lindependentintervals

    %
    time_bounds = [(timespecintervals - half_t_interval); ...
                   (timespecintervals + half_t_interval)];

    %
    ind_start_data = find(~isnan(pressuredata), 1, 'first');
%     ind_end_data = find(~isnan(pressuredata), 1, 'last');
    ind_end_data = find(~isnan(pressuredata), 1, 'last');


    %
    ind_first_wholeinterval = find(time_bounds(1, :) >= timevec(ind_start_data), 1, 'first');
    %
    ind_last_wholeinterval = find(time_bounds(2, :) < timevec(ind_end_data), 1, 'last');

    %
    ind_data_start_firstinterval = find(timevec==time_bounds(1, ind_first_wholeinterval));
    ind_data_end_lastinterval = find(timevec==time_bounds(2, ind_last_wholeinterval));

    % This MUST/SHOULD be an integer
    nt_intervals = length(ind_data_start_firstinterval : (ind_data_end_lastinterval-1)) / npts_ininterval;

    %
    ind_data_wholeintervals = ind_data_start_firstinterval + ...
                                   (0 : 1 : ((nt_intervals*npts_ininterval - 1)));
    % -1 because the round/0'th index is at the beginning and not the end of a whole interval


    %
    pressure_perinterval = reshape(pressuredata(ind_data_wholeintervals), npts_ininterval, nt_intervals);


    % ------------------------------------
% %     % Fill partial intervals at edges if
% %     % necessary
% %     if ind_first_wholeinterval ~= 1
% % 
% %     end
% %     %
% %     if ind_last_wholeinterval ~= size(time_bounds, 2)
% % 
% %     end
    % ------------------------------------

    % Update once I include partial edges
    ind_first_interval = ind_first_wholeinterval;
    ind_last_interval = ind_last_wholeinterval;
    %
    ind_all_intervals = ind_first_interval:ind_last_interval;

    % Now loop over intervals and compute spectra
    for i2 = 1:length(ind_all_intervals)

        % Reshape each interval such that each
        % column is a chunk where FFT will be
        % computed -- NO OVERLAP OR WINDOWING
        pres_chunks = reshape(pressure_perinterval(:, i2), nfft, (npts_ininterval/nfft));

        % Compute mean pressure in each chunk
        mean_pres_chunks = mean(pres_chunks, 1);

        % Remove the mean
        pres_anomaly_chunks = pres_chunks - repmat(mean_pres_chunks, nfft, 1);

        % Detrend each column with a 1st degree polynomial (i.e. a line)
        pres_detrended_chunks = detrend(pres_anomaly_chunks, 1);

% %         % Compute variance of data and detrended data
% %         varPdetrend(i) = sum((pres_anomaly_chunks .* pres_anomaly_chunks), 1)./nfft;
% %         varPdetrend(i) = sum((pres_detrended_chunks .* pres_detrended_chunks), 1)./nfft;


        % Compute variance of the full interval
        var_data_aux = sum((pres_detrended_chunks(:).^2))./nfft;


        %% Compute power spectra pressure

        % Compute fft of each column
        fourier_coefs = fft(pres_detrended_chunks, [], 1);
        
        % Select only coefficients associated with the mean
        % (0 frequency) and positive frequencies (the first
        % coefficient should be zero because the data was detrended)
        fourier_coefs = fourier_coefs(1:nnyq, :);

        % Normalize by frequency resolution
        pres_all_Spec = (fourier_coefs .* conj(fourier_coefs)) ./ (nfft * df);

        % Factor of two to account for negative frequencies
        % don't include the last one because there is no
        % negative frequency counterpart for even nfft)
        pres_all_Spec(2:end-1) = 2*pres_all_Spec(2:end-1);

        % Average across all chunks in the interval
        pres_avg_Spec = mean(pres_all_Spec, 2);

        % Normalize spectrum to maintain variance of
        % the data before it was detrended (doesn't
        % include first coefficient, though it's zero
        % after detrending).
        var_Spec_tmp = sum(pres_avg_Spec(2:end))*df;

        % Renormalization factor
        factor_renormalize = var_data_aux/var_Spec_tmp;

        % Renormalize average spectrum
        pres_avg_Spec_renormalized = pres_avg_Spec;
        pres_avg_Spec_renormalized(2:end) = pres_avg_Spec(2:end) * factor_renormalize;

% %         % Check Parseval's theorem
% %         [var_data_aux,   (sum(pres_avg_Spec_renormalized(2:end)).*df)]

        % Assign to output
        statsout.presSpec(:, ind_all_intervals(i2)) = pres_avg_Spec_renormalized;
        

%         FALK'S CODE DOES NOT USE THE NORMALIZED SPECTRUM!!! IT IS
%         COMPUTED, BUT IT IS NOT PASSED TO OUTPUT OR TO COMPUTATION
%         OF SEA SURFACE HEIGHT SPECTRUM


        %% Convert pressure spectra to sea-surface height spectra

        % What is the wavenumber k? It's straightforward to compute
        % either shallow-water or deep-water k. But what about the
        % general k? What's the numerical approach to do it?

        % % % ----------------------------------------
        % % % Go from pressure to SSH spectrum:
        % % doffz = 0.02;
        % % grav = 9.81;
        % % 
        % % k = get_wavenum(2*pi*fm, mean_depth(i));
        % % correction = cosh( k*mean_depth(i))./cosh(k*doffz);
        % % See0 = Spp0 .* correction.^2;
        % % See(i,:) = See0;
        % % % ----------------------------------------

% % % Find wavenumber from frequency -- need to
% % % loop over frequency and mean depth.
% % % Since frequency is known before the loop,
% % % should create a function handle before the
% % % loop, then only loop over mean depth
% % 
% % bla = @(x) 10*x*tanh(x*10) - (2*pi/10)^2
% % kloin = fzero(bla, [2*pi/4000, 2*pi/0.001]);

% %     %
% %     if mean_depth(i)>0.3,  % only if deep enough, calculate wave stats
% %     end
        %
    end


else
    error('Complicated option not implemented')

end

% Degrees of freedom
% DOF = 2*(nint/nfft)*2 ;  % last *2 is for the 50% overlap with the Hanning windowed data
statsout.DOF = NaN;


%% Compute bulk statistical quantities

%
inds_IGband = find((fm > freq_lims_IG(1)) & (fm < freq_lims_IG(2)));
inds_SSband = find((fm > freq_lims_SS(1)) & (fm < freq_lims_SS(2)));

%
statsout.freqlims_IG = freq_lims_IG;
statsout.PsigIG = 4*sqrt(sum(statsout.presSpec(inds_IGband, :)) * df);

%
statsout.freqlims_SS = freq_lims_SS;
statsout.PsigSS = 4*sqrt(sum(statsout.presSpec(inds_SSband, :)) * df);

%
[~, ind_peak] = max(statsout.presSpec(inds_SSband, :), [], 1);

%
statsout.peakfreqSS = statsout.frequency(inds_SSband(ind_peak));
statsout.meanfreqSS = (statsout.frequency(inds_SSband) * ...
                       statsout.presSpec(inds_SSband, :)) ./ ...
                           sum(statsout.presSpec(inds_SSband, :));






        % % HsigSS(i) = 4*sqrt( sum(See0(i_ss))*df);
        % % HsigIG(i) = 4*sqrt( sum(See0(i_ig))*df);
        % % 
        % % [mx,imx] = max(Spp0(i_ss));
        % % peakfreqSS(i) = fm(i_ss(imx));
        % % meanfreqSS(i) = sum( fm(i_ss).*See0(i_ss))/sum(See0(i_ss));



        % % DOF = 2*(nint/nfft)*2 ;  % last *2 is for the 50% overlap with the Hanning windowed data


%%







% % Pinterval = SOLOD.Pwater(ii);          % time series of pressure in interval
% % nint = length(Pinterval);
% % dummy_time_sec = [1:nint]'*dt;        % dummy time variable
% % mean_depth(i) = mean(Pinterval) + doffz;      % mean depth in interval
% % 
% % Pdemean = Pinterval - mean_depth(i);  % demeaned pressure as m
% % varP(i) = (Pdemean'*Pdemean)/nint;       % full variance of P
% % [m,b] = lsfit(dummy_time_sec,Pdemean);
% % Pdetrend = Pdemean - (m*dummy_time_sec + b);  % linear detrend
% % varPdetrend(i) = (Pdetrend'*Pdetrend)/nint;       % variance of detrend P
% % 
% % 
% % 
% % if mean_depth(i)>0.3,  % only if deep enough, calculate wave stats
% %     
% %   Ap = calculate_fft(Pdetrend,nfft,fs);   %The fourier coefficients
% %   Spp0 = mean( Ap .* conj(Ap) ) / (nfft * df);  % need to normalize by freq resolution in hz
% %   Spp0 = 2*Spp0(1:nnyq);  % keep only to the nyquist frequency
% %   tmp = sum(Spp0)*df;
% %   Spp00 = Spp0 * varPdetrend(i)/tmp;   %normalize to keep detrended variance
% %                                        %    [sum(Spp00)  sum(Spp0)]
% %   Spp(i,:) = Spp0;



%%
%%

% % ii = find((SOLOD.time_dnum>=timeL2_dnum(i)-half_t_interval) &  (SOLOD.time_dnum<timeL2_dnum(i)+half_t_interval));
% % 
% %     Pinterval = SOLOD.Pwater(ii);          % time series of pressure in interval
% %     nint = length(Pinterval);
% %     dummy_time_sec = [1:nint]'*dt;        % dummy time variable
% %     mean_depth(i) = mean(Pinterval) + doffz;      % mean depth in interval
% % 
% %     Pdemean = Pinterval - mean_depth(i);  % demeaned pressure as m
% %     varP(i) = (Pdemean'*Pdemean)/nint;       % full variance of P
% %     [m,b] = lsfit(dummy_time_sec,Pdemean);
% %     Pdetrend = Pdemean - (m*dummy_time_sec + b);  % linear detrend
% %     varPdetrend(i) = (Pdetrend'*Pdetrend)/nint;       % variance of detrend P
% % 
% %     clf
% %     subplot(3,1,1)
% %     plot(dummy_time_sec,Pinterval),grid
% %     xlabel('time of the hour (sec)')
% %     ylabel('water pressure (m)')
% %     title(sprintf('ID=%s: %s  depth=%.2f (m)',SOLOD.ID,datestr(timeL2_dnum(i),0),mean_depth(i)))
% % 
% %     if mean_depth(i)>0.3,  % only if deep enough, calculate wave stats
% %         
% %       Ap = calculate_fft(Pdetrend,nfft,fs);   %The fourier coefficients
% %       Spp0 = mean( Ap .* conj(Ap) ) / (nfft * df);  % need to normalize by freq resolution in hz
% %       Spp0 = 2*Spp0(1:nnyq);  % keep only to the nyquist frequency
% %       tmp = sum(Spp0)*df;
% %       Spp00 = Spp0 * varPdetrend(i)/tmp;   %normalize to keep detrended variance
% %                                            %    [sum(Spp00)  sum(Spp0)]
% %       Spp(i,:) = Spp0;
% %     
% % 
% %       doffz = 0.02;
% %       grav = 9.81;
% %       k = get_wavenum(2*pi*fm, mean_depth(i));
% %       correction = cosh( k*mean_depth(i))./cosh(k*doffz);
% %       See0 = Spp0 .* correction.^2;
% %       See(i,:) = See0;
% %     
% %       HsigSS(i) = 4*sqrt( sum(See0(i_ss))*df);
% %       HsigIG(i) = 4*sqrt( sum(See0(i_ig))*df);
% % 
% %       [mx,imx] = max(Spp0(i_ss));
% %       peakfreqSS(i) = fm(i_ss(imx));
% %       meanfreqSS(i) = sum( fm(i_ss).*See0(i_ss))/sum(See0(i_ss));
% %       DOF = 2*(nint/nfft)*2 ;  % last *2 is for the 50% overlap with the Hanning windowed data
% % 
% % 
