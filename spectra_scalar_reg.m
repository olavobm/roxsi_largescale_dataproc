function [Sxx, timespec, freqvec, dof, avgx] = spectra_scalar_reg(xtime, xdata, windowfft, windowavg, timespeclims, dt_step)
%% [Sxx, freqspec, dof] = SPECTRA_SCALAR_REG(xtime, xdata, windowfft, windowavg, timespeclims, dt_step)
%
%   inputs
%       - xtime:
%       - xdata:
%       - windowfft: in seconds
%       - windowavg: in seconds
%       - timespec (optional, but recommended):
%
%   outputs
%       - Sxx:
%       - timespec
%       - freqspec
%       - dof
%       - avgx
%
%
% SPECTRA_SCALAR_REG.m
%
%

%% Check inputs

%
if ~isdatetime(xtime)
    error('Time vector of the data is NOT in datetime format.')
end
% % %
% % if ~isdatetime(timespec)
% %     error('Time grid is NOT in datetime format.')
% % end

% % %
% % if isempty(timespec.TimeZone)
% %     timespec.TimeZone = xtime.TimeZone;
% % else
% %     %
% %     if ~isequal(xtime.TimeZone, timespec.TimeZone)
% %         error('Time zones don''t match.')
% %     end
% % end

%% check data is time-gridded!!!

%%

%
if ~exist('fracoverlap', 'var')
% %     fracoverlap = 0;
    fracoverlap = 0.5;    % EITHER 0 OR 0.5. NO OTHER OPTION IMPLEMENTED
end

%
if ~exist('dt_step', 'var')
    %
    dt_step = windowavg/(24*3600);    % window in units of days/datenum
end


%
if ~exist('timespeclims', 'var')
    
    %
    tstart_dnum = xtime(1) + dt_step;
    tend_dnum = xtime(end) - dt_step;
    %
    timespeclims = [tstart_dnum, tend_dnum];

end


%%
% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------------------------------------------------------------------


%%

%
timespec = timespeclims(1): dt_step : timespeclims(2);
%
ntspec = length(timespec);
half_t_interval = dt_step/2;

%
time_bounds = [(timespec - half_t_interval); ...
               (timespec + half_t_interval)];


%%

% Time difference, in seconds, of the data time grid
% (the code assumes it is a equally time grid)
dt = (xtime(2) - xtime(1));
dt = seconds(dt);

%
fs = 1/dt;     % sampling frequency in Hertz
fnyq = fs/2;   % Nyquist frequency


%% The number of points per averaging window

%
npts_ininterval = fs*windowavg;


%%

% Number of points per chunk
nfft = windowfft*fs;

% Only allow for even nfft (you could want otherwise
% is sampling period is 1.5s, but that seems unlikely)
if mod(nfft, 2)~=0
    error('nfft is odd. It must be even.')
end

%
df = fs/(nfft-1);   % This is the frequency resolution (in Hz)
nnyq = nfft/2 + 1;

% Positive frequency vector (in Hz)
fm = (0:(nnyq-1)) * df;



%%
% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------------------------------------------------------------------


%%

%
ind_start_data = find(~isnan(xdata), 1, 'first');
ind_end_data = find(~isnan(xdata), 1, 'last');


%%

%
ind_first_wholeinterval = find(time_bounds(1, :) >= xtime(ind_start_data), 1, 'first');
%
ind_last_wholeinterval = find(time_bounds(2, :) < xtime(ind_end_data), 1, 'last');

%
ind_data_start_firstinterval = find(xtime==time_bounds(1, ind_first_wholeinterval));
ind_data_end_lastinterval = find(xtime==time_bounds(2, ind_last_wholeinterval));


%%

% Number of whole intervals (i.e. whole intervals
% means that it doesn't include data at the edges)
nt_whole_intervals = length(ind_data_start_firstinterval : (ind_data_end_lastinterval-1)) / npts_ininterval;

%
ind_data_wholeintervals = ind_data_start_firstinterval + ...
                               (0 : 1 : ((nt_whole_intervals*npts_ininterval - 1)));
% -1 because the round/0'th index is at the beginning and not the end of a whole interval

%
xdata_perinterval = reshape(xdata(ind_data_wholeintervals), npts_ininterval, nt_whole_intervals);


%% IF there is more data at the edges, add them now

%
if ind_first_wholeinterval ~= 1
    %
    inds_edge_1 = ind_start_data:(ind_data_start_firstinterval - 1);

    %
    xdata_perinterval = [NaN(npts_ininterval, 1), xdata_perinterval];
    %
    inds_fill_aux = (npts_ininterval - (length(inds_edge_1) - 1)) : npts_ininterval;
    xdata_perinterval(inds_fill_aux, 1) = xdata(inds_edge_1);
    %
    ind_first_interval = ind_first_wholeinterval - 1;

else

    ind_first_interval = ind_first_wholeinterval;

end


%
if ind_last_wholeinterval ~= size(time_bounds, 2)
    %
    inds_edge_2 = (ind_data_wholeintervals(end) + 1):ind_end_data;

    %
    xdata_perinterval = [xdata_perinterval, NaN(npts_ininterval, 1)];
    %
    inds_fill_aux = 1:length(inds_edge_2);
    xdata_perinterval(inds_fill_aux, end) = xdata(inds_edge_2);
    %
    ind_last_interval = ind_last_wholeinterval + 1;
%
else

    ind_last_interval = ind_last_wholeinterval;


end


%%

%
ind_all_intervals = ind_first_interval:ind_last_interval;


%% Create indices to deal with overlap

%
indstart_firstchunk = 1;
indend_lastchunk = npts_ininterval;

%
if (fracoverlap==0.5)
    %
    indstart_firstchunk = [indstart_firstchunk, (nfft/2)];
    %
    indend_lastchunk = [indend_lastchunk, (npts_ininterval - (1 + (nfft*fracoverlap)))];
end

% Compute number of individual chunks for which fft will be computed
indnumberchunks = (indend_lastchunk - indstart_firstchunk + 1)./nfft;
%
indcols = NaN(length(indnumberchunks), 2);
indcols(1, 1) = 1;
%
inds_cumsum = cumsum(indnumberchunks);
%
indcols(2:end, 1) = inds_cumsum(1:end-1)+1;
indcols(:, 2) = inds_cumsum;


keyboard

for i = 1:length(indstart_firstchunk)
    %
    ind_edge_chunk_aux = indstart_firstchunk(i) : nfft : indend_lastchunk(i);

    keyboard
    if ind_edge_chunk_aux(end)==npts_ininterval
        
    end
% %     %
% %     indnumberchunks = length(indstart_firstchunk);
% %     %
% %     indcols = []
end

% Compute upper bound of DOF (as if the chunks were fully independent)

keyboard

%%
% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------------------------------------------------------------------

%%

%
freqvec = fm;

%
Sxx = NaN(length(fm), ntspec);
dof = NaN(1, ntspec);


%%

%
avgx = mean(xdata_perinterval, 1, 'omitnan');


%% Finally start the calculation

% Loop over averaging intervals
for i1 = 1:length(ind_all_intervals)  % == ntspec (?)


    %% Get the data, breaking it apart in chunks for computing fft
    
    %
    xdata_chunks = reshape(xdata_perinterval(:, i1), nfft, (npts_ininterval/nfft));

    %
    if length(indstart_firstchunk)~=1

        %
% %         pres_chunks = [pres_chunks, NaN(???)]

        %
        for i2 = 2:length(indstart_firstchunk)

            %
            inds_sub_aux = indstart_firstchunk(i2):indend_lastchunk(i2);

            %
            x_chunks_overlap_aux = reshape(xdata_perinterval(inds_sub_aux, i1), nfft, (npts_ininterval/nfft));

% %             %
% % % %             pres_chunks(:, indcols(i2)) = pres_chunks_overlap_aux;
        end

    end
    
    

    %% Preliminary calculations to compute fft

    % Compute mean pressure in each chunk
    mean_xdata_chunks = mean(xdata_chunks, 1);

    % Remove the mean
    xdata_anomaly_chunks = xdata_chunks - repmat(mean_xdata_chunks, nfft, 1);

    % Detrend each column with a 1st degree polynomial (i.e. a line)
    xdata_forfft_chunks = detrend(xdata_anomaly_chunks, 1);

    % Apply window
%     window_vec = 
%     window_matrix = repmat(window_vec, 1, )
    xdata_forfft_chunks = xdata_forfft_chunks .* window_matrix;

% %         % Compute variance of data and detrended data
% %         varPdetrend(i) = sum((pres_anomaly_chunks .* pres_anomaly_chunks), 1)./nfft;
% %         varPdetrend(i) = sum((pres_forfft_chunks .* pres_forfft_chunks), 1)./nfft;


    % Compute variance of the full interval (to compute the
    % normalization factor so that the average spectrum has
    % the same variance as the full interval)
    var_data_aux = sum((xdata_forfft_chunks(:).^2)) ./ numel(xdata_forfft_chunks);

    %
    dof(i1) = 2*length(~isnan(mean_xdata_chunks));


    %% Compute power spectra

    % Compute fft of each column
    fourier_coefs = fft(xdata_forfft_chunks, [], 1);
    
    % Select only coefficients associated with the mean
    % (0 frequency) and positive frequencies (the first
    % coefficient should be zero because the data was detrended)
    fourier_coefs = fourier_coefs(1:nnyq, :);

    % Normalize by frequency resolution
    x_all_Spec = (fourier_coefs .* conj(fourier_coefs)) ./ (nfft * df);

    % Factor of two to account for negative frequencies
    % don't include the last one because there is no
    % negative frequency counterpart for even nfft)
    x_all_Spec(2:end-1) = 2*x_all_Spec(2:end-1);

    % Average across all chunks in the interval
    x_avg_Spec = mean(x_all_Spec, 2);

    % Normalize spectrum to maintain variance of
    % the data before it was detrended (doesn't
    % include first coefficient, though it's zero
    % after detrending).
    var_Spec_tmp = sum(x_avg_Spec(2:end))*df;

    % Renormalization factor
    factor_renormalize = var_data_aux/var_Spec_tmp;

    % Renormalize average spectrum
    x_avg_Spec_renormalized = x_avg_Spec;
    x_avg_Spec_renormalized(2:end) = x_avg_Spec(2:end) * factor_renormalize;

% %         % Check Parseval's theorem
% %         [var_data_aux,   (sum(x_avg_Spec_renormalized(2:end)).*df)]

    %%

    % Assign to output
    Sxx(:, ind_all_intervals(i2)) = x_avg_Spec_renormalized;
    
%         FALK'S CODE DOES NOT USE THE NORMALIZED SPECTRUM!!! IT IS
%         COMPUTED, BUT IT IS NOT PASSED TO OUTPUT OR TO COMPUTATION
%         OF SEA SURFACE HEIGHT SPECTRUM



end



