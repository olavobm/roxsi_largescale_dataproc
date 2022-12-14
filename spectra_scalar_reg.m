function [Sxx, timespec, freqvec, dof, avgx] = spectra_scalar_reg(xtime, xdata, windowfft, windowavg, timespeclims, dt_step)
%% [Sxx, timespec, freqvec, dof, avgx] = SPECTRA_SCALAR_REG(xtime, xdata, windowfft, windowavg, timespeclims, dt_step)
%
%   inputs
%       - xtime: time grid (in date time format) of the data xdata.
%       - xdata: a vector with data to compute spectra.
%       - windowfft: time window (in seconds) to compute fft.
%       - windowavg: time window (in seconds) for averaging individual
%                    spectra (MUST BE MULTIPLE OF windowfft).
%       - timespeclims (optional): datetime limits.
%       - dt_step (optional):
%
%   outputs
%       - Sxx: spectra of xdata. A matrix where row-dimension if frequency
%              and comlumn-dimension is time.
%       - timespec: time (datetime) vector with same number
%                   of columns as Sxx.
%       - freqspec: frequency vector (in Hz) with same number
%                   of rows as Sxx.
%       - dof: vector with degrees of freedom for each spectrum in Sxx.
%       - avgx: average of xdata within each window windowavg.
%
%
% SPECTRA_SCALAR_REG.m computes spectra from timeseries xdata. One spectrum
% is computed for each averaging window windowavg and the resulting spectra
% are given by the columns of Sxx. The frequencies (row-dimension of Sxx)
% are given by output freqvec. The times associated with each column of
% column of Sxx are given by output timespec. The data for each spectrum
% in Sxx is centered at timespec (i.e. 50% of windowavg is taken to the
% left and to right of timespec).
%
% Within each averaging window, the timeseries is split in chunks with
% length windowfft. This function can handle either no overlap or 50%
% overlap between these chunks. The averaged spectrum from all of these chunks
% is computed by ignoring chunks with any gaps. The spectra computed from
% gappy data have smaller degrees of freedom (given in output dof).
%
% Spectra are computed between time limits given by timespeclims. If this
% is not given in the input, the function neglects a bit of the edges,
% but uses all the rest of the data. If timespeclims are beyond the times
% when there are data, then the corresponding edge columns of Sxx will be
% NaNs. The default time stepping is calculated from windowavg to use all
% of the data, but this can be given as input dt_step (though I need to
% test the function in this case!!!).
%
% The idea behind SPECTRA_SCALAR_REG.m is to compute spectra for a
% regular (_reg / normal / time-gridded) timeseries. Therefore, xtime
% must be a vector with constant grid step (dt). This idea is also 
% associated with having input timespeclims instead of a vector of
% times for when spectra should be computed. Therefore, this code
% is able to select data and compute fft's in an efficient way.
%
% Olavo Badaro Marques, Oct/2022.


%%
% ------------------------------------------------------------------
% ---------------------------- PREAMBLE ----------------------------
% ------------------------------------------------------------------

%% Check inputs

%
if ~isdatetime(xtime)
    error('Time vector of the data is NOT in datetime format.')
end

% #########################
% check data is time-gridded!!!
% #########################


%%

%
if ~exist('fracoverlap', 'var')
    fracoverlap = 0;
%     fracoverlap = 0.5;    % EITHER 0 OR 0.5. NO OTHER OPTION IMPLEMENTED
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

else

    %
    if isempty(timespeclims.TimeZone)
        %
        timespeclims.TimeZone = xtime.TimeZone;
    else
        %
        if ~isequal(timespeclims.TimeZone, xtime.TimeZone)
            error("Timezones don't match")
        end
    end
end


%% Windows for fft

% % window_handle = @rectwin;
window_handle = @hann;


%%
% ------------------------------------------------------------------
% ----------- DEFINE BASIC PARAMETERS OF THE CALCULATION -----------
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

%
window_weights = window(window_handle, nfft);    % this is already a column vector

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
% --------- DEFINE VARIABLES TO TAKE CARE OF SPLITTING THE ---------
% ------ DATA IN CHUNKS (INCLUDING THE OVERLAP) TO COMPUT FFT ------
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


%% If ind_data_wholeintervals is empty, there is no or little data,
% so function should stop here.
% 
% HOWEVER THERE ARE OTHER WAYS WHERE
% GAPS MAY NOT BE CAPTURED AND THE FUNCTION WILL BREAK LATER BELOW

%
if isempty(ind_data_wholeintervals)
    %
    Sxx = NaN;
    timespec = NaN;
    freqvec = NaN;
    dof = NaN;
    avgx = NaN;
    %
    return
end


%% Reshape data within whole intervals from a vector to a matrix

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
% (the coding here only allows for 0 or 50% overlap
% and the coding assumes 1/overlap is an interget)

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

% Compute upper bound of DOF (as if the chunks were fully independent)
% % 2*inds_cumsum(end)




%%
% ------------------------------------------------------------------
% ---------------------- NOW COMPUTE SPECTRA -----------------------
% ------------------------------------------------------------------

%%

%
freqvec = fm;

%
Sxx = NaN(length(fm), ntspec);
dof = NaN(1, ntspec);


%%

%
avgx = NaN(1, ntspec);
avgx(ind_all_intervals) = mean(xdata_perinterval, 1, 'omitnan');


%% Finally start the calculation

% Loop over averaging intervals
for i1 = 1:length(ind_all_intervals)

    %% Get the data, breaking it apart in chunks for computing fft
    % (note that overlapping segments are concatenated at the end
    % of matrix and not adjacent)
    
    %
    xdata_chunks = reshape(xdata_perinterval(:, i1), nfft, (npts_ininterval/nfft));

    % If there is overlap, then get the other segments
    if length(indstart_firstchunk)~=1

        %
        xdata_chunks = [xdata_chunks, NaN(nfft, (inds_cumsum(end) - inds_cumsum(1)))];

        %
        for i2 = 2:length(indstart_firstchunk)

            %
            inds_sub_aux = indstart_firstchunk(i2):indend_lastchunk(i2);

            %
            xdata_chunks(:, indcols(i2, 1):indcols(i2, 2)) = ...
                        reshape(xdata_perinterval(inds_sub_aux, i1), nfft, (length(inds_sub_aux)/nfft));
        end
    end
    

    %% Preliminary calculations to compute fft

    % Compute mean pressure in each chunk
    mean_xdata_chunks = mean(xdata_chunks, 1);

    % Remove the mean
    xdata_anomaly_chunks = xdata_chunks - repmat(mean_xdata_chunks, nfft, 1);

    % Detrend each column with a 1st degree polynomial (i.e. a line)
    xdata_forfft_chunks = detrend(xdata_anomaly_chunks, 1);

    % Compute variance of the full interval (to compute the  
    % normalization factor so that the average spectrum has
    % the same variance as the full interval) [not exactly right when
    % there are gaps, but shouldn't be a problem]
    var_data_aux = var(xdata_forfft_chunks(:), 'omitnan');

    % Apply window
    window_matrix = repmat(window_weights, 1, size(xdata_forfft_chunks, 2));
    xdata_forfft_chunks = xdata_forfft_chunks .* window_matrix;

    % Compute upper bound of degrees of freedom (i.e. in the sense
    % % that this upper bound would be the true dof if segments were
    % fully independent).
    % Gaps in chunks will decrease the degrees of freedom
    dof(ind_all_intervals(i1)) = 2*length(mean_xdata_chunks(~isnan(mean_xdata_chunks)));


    %% Compute power spectra

    % Compute fft of each column
    fourier_coefs = fft(xdata_forfft_chunks, [], 1);
    
    % Select only coefficients associated with the mean
    % (0 frequency) and positive frequencies (the first
    % coefficient should be zero because the data was detrended)
    fourier_coefs = fourier_coefs(1:nnyq, :);

    % Normalize by frequency resolution
    x_all_Spec = (fourier_coefs .* conj(fourier_coefs)) ./ ((nfft^2) * df);

    % Factor of two to account for negative frequencies
    % don't include the last one because there is no
    % negative frequency counterpart for even nfft) -- CHECK WHETHER IT'S THE OPPOSITE!!!!!
    x_all_Spec(2:end-1) = 2*x_all_Spec(2:end-1);

    % Average across all chunks in the interval and ignore chunks
    % with gaps (look at output dof to see which intervals are
    % complete or have gaps)
    x_avg_Spec = mean(x_all_Spec, 2, 'omitnan');

    % Normalize spectrum to maintain variance of
    % the data before it was detrended (doesn't
    % include first coefficient, though it's zero
    % after detrending).
    var_Spec_tmp = sum(x_avg_Spec(2:end))*df;

    % Renormalization factor due to windowing
    factor_renormalize = var_data_aux/var_Spec_tmp;

    % Renormalize average spectrum
    x_avg_Spec_renormalized = x_avg_Spec;
    x_avg_Spec_renormalized(2:end) = x_avg_Spec(2:end) * factor_renormalize;

% %         % Check Parseval's theorem
% %         [var_data_aux,   (sum(x_avg_Spec_renormalized(2:end)).*df)]

    
    %%

    % Assign to output
    Sxx(:, ind_all_intervals(i1)) = x_avg_Spec_renormalized;
    


end



