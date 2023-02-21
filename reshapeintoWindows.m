function [indsub, reshapeNdims] = reshapeintoWindows(tdata, tgrid)
%% [indsub, reshapeNdims] = RESHAPEINTOWINDOWS(tdata, tgrid)
%
%   inputs
%       - tdata: independent variable of the gridded data.
%       - tgrid: equally spaced grid for windowed calculations (where the
%                window length is diff(tgrid)).
%
%   outputs
%       - indsub:
%       - reshapeNdims:
%
%
% RESHAPEINTOWINDOWS.m
%
%
%
%
%
% Olavo Badaro Marques, 21/02/2023.


%% Check dtime is gridded

%
dt = diff(tdata(1:2));

%
if any(diff(tdata) ~= dt)
    error('Data grid is not an equally spaced grid.')
end


%%

%
windowlen = diff(tgrid(1:2));

%
if isdatetime(tdata)
    %
    npts_window = (1/seconds(dt)) * seconds(windowlen);
%
else
    %
    npts_window = (1/dt) * (windowlen);
end


%%

%
time_bounds = [(tgrid(:) - (windowlen/2)).'; ...
               (tgrid(:) + (windowlen/2)).'];

%
ind_first_wholeinterval = find(time_bounds(1, :) >= tdata(1), 1, 'first');
ind_last_wholeinterval = find(time_bounds(2, :) < tdata(end), 1, 'last');


%%

%
ind_data_start_firstwindow = find(tdata == time_bounds(1, ind_first_wholeinterval));
%
ind_data_end_lastwindow = find(tdata == time_bounds(2, ind_last_wholeinterval));
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


