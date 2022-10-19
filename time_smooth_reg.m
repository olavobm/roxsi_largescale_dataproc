function [avgx, dtimeavg, nptsavg] = time_smooth_reg(xtime, xdata, windowavg, timelims, dt_step)
%% [avgx, dtimeavg, nptsavg] = TIME_SMOOTH_REG(xtime, xdata, windowavg, timelims, dt_step)
%
%   inputs
%       -
%       -
%       -
%       -
%       -
%
%   outputs
%       -
%       -
%       -
%
%
%
%
%
%
%
% See also: SPECTRA_SCALAR_REG.m


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
if ~exist('dt_step', 'var')
    %
    dt_step = windowavg/(24*3600);    % window in units of days/datenum
end


%
if ~exist('timelims', 'var')
    
    %
    tstart_dnum = xtime(1) + dt_step;
    tend_dnum = xtime(end) - dt_step;
    %
    timelims = [tstart_dnum, tend_dnum];

else

    %
    if isempty(timelims.TimeZone)
        %
        timelims.TimeZone = xtime.TimeZone;
    else
        %
        if ~isequal(timelims.TimeZone, xtime.TimeZone)
            error("Timezones don't match")
        end
    end
end


%% Window for averaging???

% % % % window_handle = @rectwin;
% % window_handle = @hann;


%%
% ------------------------------------------------------------------
% ----------- DEFINE BASIC PARAMETERS OF THE CALCULATION -----------
% ------------------------------------------------------------------


%%

%
dtimeavg = timelims(1): dt_step : timelims(2);
%
nptsgrid = length(dtimeavg);
half_t_interval = dt_step/2;

%
time_bounds = [(dtimeavg - half_t_interval); ...
               (dtimeavg + half_t_interval)];


%%

% Time difference, in seconds, of the data time grid
% (the code assumes it is a equally time grid)
dt = (xtime(2) - xtime(1));
dt = seconds(dt);

%
fs = 1/dt;     % sampling frequency in Hertz


%% The number of points per averaging window

%
npts_ininterval = fs*windowavg;


%%

% % % % % Number of points per chunk
% % % % nfft = windowfft*fs;
% % 
% % %
% % window_weights = window(window_handle, nfft);    % this is already a column vector
% % 
% % % Only allow for even nfft (you could want otherwise
% % % is sampling period is 1.5s, but that seems unlikely)
% % if mod(nfft, 2)~=0
% %     error('nfft is odd. It must be even.')
% % end
% % 
% % %
% % df = fs/(nfft-1);   % This is the frequency resolution (in Hz)
% % nnyq = nfft/2 + 1;
% % 
% % % Positive frequency vector (in Hz)
% % fm = (0:(nnyq-1)) * df;



%%
% ------------------------------------------------------------------
% ------- SPLIT DATA IN AVERAGING WINDOWS TO COMPUTE AVERAGES ------
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
    avgx = NaN;
    dtimeavg = NaN;
    nptsavg = NaN;
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

% THIS IMPORTANT IF THE TIMLE LIMITS IN THE OUTPUT
% GO BEYONG WHERE THERE IS ENOUGH DATA. BUT IN THIS CASE
% THE TIME AND X DATA IN OUTPUT WILL HAVE DIFFERENT NUMBER
% COLUMNS (AND I PROBABLY WANT THE SAME NUMBER OF COLUMNS)


%% Create indices to deal with overlap
% % % %
% % % % (the coding here only allows for 0 or 50% overlap
% % % % and the coding assumes 1/overlap is an interget)
% % % 
% % % %
% % % indstart_firstchunk = 1;
% % % indend_lastchunk = npts_ininterval;
% % % 
% % % %
% % % if (fracoverlap==0.5)
% % %     %
% % %     indstart_firstchunk = [indstart_firstchunk, (nfft/2)];
% % %     %
% % %     indend_lastchunk = [indend_lastchunk, (npts_ininterval - (1 + (nfft*fracoverlap)))];
% % % end
% % % 
% % % % Compute number of individual chunks for which fft will be computed
% % % indnumberchunks = (indend_lastchunk - indstart_firstchunk + 1)./nfft;
% % % %
% % % indcols = NaN(length(indnumberchunks), 2);
% % % indcols(1, 1) = 1;
% % % %
% % % inds_cumsum = cumsum(indnumberchunks);
% % % %
% % % indcols(2:end, 1) = inds_cumsum(1:end-1)+1;
% % % indcols(:, 2) = inds_cumsum;
% % % 
% % % % Compute upper bound of DOF (as if the chunks were fully independent)
% % % % % 2*inds_cumsum(end)

%%
% ------------------------------------------------------------------
% ---------------------- NOW COMPUTE AVERAGES ----------------------
% ------------------------------------------------------------------

%%

%
nptsavg = NaN(1, length(dtimeavg));
avgx = nptsavg;


%% Compute averages and count number of points per averaging window

%
avgx(ind_all_intervals) = mean(xdata_perinterval, 1, 'omitnan');

%
lgooddata = ~isnan(xdata_perinterval);
%
nptsavg(ind_all_intervals) = sum(lgooddata, 1);




