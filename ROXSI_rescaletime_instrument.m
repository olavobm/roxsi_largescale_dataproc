function time_vec_corrected = ROXSI_rescaletime_instrument(table_correction, SN, time_vec)
%% time_vec_corrected = ROXSI_RESCALETIME_INSTRUMENT(table_correction, SN, time_vec)
%
%   inputs
%       - table_correction: table defined by a separate script.
%       - SN: serial number of instrument of interest.
%       - time_vec: time vector (in datenum format) to be corrected.
%
%   outputs
%       - time_vec_corrected: time vector corrected with a
%                             linear drift assumption.
%
% ROXSI_RESCALETIME_INSTRUMENT.m is a higher level function that calls
% the function that corrects for clock drift.
%
%
%


%% Find row of the table that matches the SN input

%
ind_row_match = find(strcmp(table_correction.SN, SN));

%
if isempty(ind_row_match)
    error(['Instrument with SN = ' char(SN) ' not found in clock correction table.'])
end


%% Apply linear clock drift correction

% Apply clock drift correction if clock drift info is available
%
if isnan(table_correction.clockdrift(ind_row_match))
    % If clock drift is NaN, then there is no clock correction
    % to be applied. The output is 

    %
    warning(['No clock drift correction available for SN = ' char(SN)])

    %
    time_vec_corrected = time_vec;

else

    % Convert time from string to datenum
    str_format = 'yyyy/mm/dd HH:MM:SS';
    %
    time_set_datenum = datenum(table_correction.time_set_clockdrift(ind_row_match), str_format);
    time_end_datenum = datenum(table_correction.time_end_clockdrift(ind_row_match), str_format);
    
    % Apply linear clock drift correction.
    % NOTE THE MINUS SIGN (!). The function defines clock drift "properly",
    % that is, in the sense that a negative clock drift means that the
    % instrument time is behind the actual (reference) time. However,
    % Aquadopps and Signature1000s adopt the opposite definition.
    time_vec_corrected = clockdrift_linear_correction(time_vec, ...
                                                      time_set_datenum, time_end_datenum, ...
                                                      -table_correction.clockdrift(ind_row_match));

end