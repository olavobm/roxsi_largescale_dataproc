function s = Spotter_parsed_to_mat(dirparsed, timezone_UTCtoLocal)
%% s = SPOTTER_PARSED_TO_MAT(dirparsed, timezone_UTCtoLocal)
%
%   inputs
%       - dirparsed: directory with Spotter data (as parsed by
%                    Sofar's parsing Python script).
%       - timezone_GMTtoLocal (optional): an integer that converts from
%                                         UTC to local time.
%
%   outputs
%       - s: structure variable with Spotter buoy variables.
%
%
% SPOTTER_PARSED_TO_MAT.m
%
% Olavo Badaro Marques


%% Check if time zone correction 

%
if ~exist('timezone_UTCtoLocal', 'var')
    timezone_UTCtoLocal = 0;
    %
    warning(['Time zone correction not given in input. ' ...
             'No change in timezone has been applied.'])
end


%% Loop over *.csv files in put all of them in a structure variable

%
csv_fileinfo = dir(fullfile(dirparsed, '*.csv'));

% Loop over csv files
for i = 1:length(csv_fileinfo)

    % File name without csv extension
    filename_aux = csv_fileinfo(i).name(1:end-4);
    
    % Read all files, except for system.csv
    if filename_aux~="system"

        %
        file_fullpath_aux = fullfile(dirparsed, csv_fileinfo(i).name);
        % Print on the screen the file that is being read
        fprintf('Reading %s\n', file_fullpath_aux)

        %
        a = readtable(file_fullpath_aux, 'ReadVariableNames', true, ...
                                         'VariableNamingRule', 'preserve');
        a = renamevars(a, "# year", "year");

        % Rename vars in displacement so that x, y, and z are consistent
        if filename_aux=="displacement"
            % Check if the variable name appears as I expect (Sofar
            % may update their parsing script at some point)
            if any(strcmp(a.Properties.VariableNames, 'y(m)'))
                a = renamevars(a, "y(m)", "y (m)");
            end
            if any(strcmp(a.Properties.VariableNames, 'z(m)'))
                a = renamevars(a, "z(m)", "z (m)");
            end
        end

        % Compute datetime
        if filename_aux=="location" || filename_aux=="displacement" || filename_aux=="sst"
            var_milisecond = a.msec;
        else
            var_milisecond = a.milisec;
        end
        %
        time = datetime(a.year, a.month, a.day, a.hour + timezone_UTCtoLocal, a.min, a.sec, var_milisecond);
            
        % Add datetime before the variable/column year
        a = addvars(a, time, 'Before','year');
        s.(filename_aux) = a;
    end

end
