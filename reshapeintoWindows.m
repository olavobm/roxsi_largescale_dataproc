function [indsub, reshapeNdims] = reshapeintoWindows(dtime, dtimegrid, windowlen)
%% [indsub, reshapeNdims] = reshapeintoWindows(dtime, dtimegrid, windowlen)
%
%   inputs
%       -
%       -
%       -
%       -
%
%   outputs
%       - indsub:
%       - reshapeNdims:
%
%
%
%
%
%
%
%
% Olavo Badaro Marques, 21/02/2023.


%% Check dtime is gridded

%
dt = diff(dtime(1:2));

%
if any(diff(dtime) ~= dt)
    error('Data grid is not an equally spaced grid.')
end


%%

%
if isdatetime(dtime)
    %
    npts_window = (1/seconds(dt)) * seconds(windowlen);
%
else
    %
    npts_window = (1/dt) * (windowlen);
end


%%

%
time_bounds = [(dtimegrid(:) - (windowlen/2)).'; ...
               (dtimegrid(:) + (windowlen/2)).'];

%
ind_first_wholeinterval = find(time_bounds(1, :) >= dtime(1), 1, 'first');
ind_last_wholeinterval = find(time_bounds(2, :) < dtime(end), 1, 'last');


%%

%
ind_data_start_firstwindow = find(dtime == time_bounds(1, ind_first_wholeinterval));
%
ind_data_end_lastwindow = find(dtime == time_bounds(2, ind_last_wholeinterval));
ind_data_end_lastwindow = ind_data_end_lastwindow - 1;

%
indsub = ind_data_start_firstwindow : ind_data_end_lastwindow;


%%

%
Nindsub = length(indsub);

%
if mod(Nindsub, npts_window)~=0
    error('Integer not obtained. Reshape failed.')
end

%
npts_intervals = Nindsub/npts_window;

%
reshapeNdims = [npts_window, npts_intervals];


