function time_true_shift = clockdrift_linear_correction(time_inst, clock_set_time, clock_end_time, clock_drift)
%% time_true_shift = CLOCKDRIFT_LINEAR_CORRECTION(time_inst, clock_set_time, clock_end_time, clock_drift)
%
%   inputs:
%       - time_inst: instrument time vector (in datenum format, 
%                    or where units are in days).
%       - clock_set_time: time when instrument was synced to reference.
%       - clock_end_time: time when clock drift was recorded.
%       - clock_drift: number of seconds sensor clock has drifted
%                      (negative means sensor is behind reference time).
%
%   outputs:
%       - time_true_shift: instrument time vector corrected by drift.
%
%
% CLOCKDRIFT_LINEAR_CORRECTION.m correct the time vector of an
% instrument assuming a linear clock drift relative to a
% reference time source.


%%

% instrument time since programmed (in sec)
time_sec = (time_inst - clock_set_time)*86400;

% Number of seconds measured by the instrument
% between set and end times
dt = (clock_end_time - clock_set_time)*86400;

% True number of seconds between set and end times
dt_true = dt - clock_drift;


%% Slope of linear time drift correction

%
time_factor2 = dt_true/dt;

%
fprintf('Time dilation factor is %.8f \n', time_factor2)

% True number of seconds since clock was set
time_true = time_sec * time_factor2;


%% Timestamps of the instrument corrected by clock drift

%
time_true_shift = time_true/86400 + clock_set_time;

