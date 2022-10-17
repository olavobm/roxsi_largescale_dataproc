%% Time grid pressure data from all smart moorings
% in ROXSI 2022 to a common time grid.

%
clear
close all


%% Data directory

%
% dir_data_L1 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/code_proc';
% dir_data_L1 = pwd;
%
dir_data_L1 = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Smart_Level1/not_gridded/';


%%

dir_output = fullfile(dir_data_L1, '..', 'gridded');


%% Get all L1 file names in folder

%
dir_all_L1 = dir(fullfile(dir_data_L1, '*_notgridded.mat'));
%
list_files = cell(1, length(dir_all_L1));

% Put file names in a cell array
for i = 1:length(list_files)
    %
    list_files{i} = dir_all_L1(i).name;
end

%
disp('------------ Will time grid pressure for the following smart moorings: ------------ ')
%
list_files


%% Loop over file names and load data

for i = 1:length(list_files)

    %
    data_loaded = load(fullfile(dir_data_L1, list_files{i}));

    %
    smartMoorlvl1.([list_files{i}(15:17) '_' list_files{i}(21:24)]) = data_loaded.spotsmartL1;

end

%
list_fields = fieldnames(smartMoorlvl1);


%% Get min/max date-times where observations are available

% --------------------------------------------------------
% Min/max time grabbing all data

% % %
% % min_max_times = [min(smartMoorlvl1.(list_fields{1}).dtime), ...
% %                  max(smartMoorlvl1.(list_fields{1}).dtime)];
% % 
% % % Get min/max time
% % for i = 2:length(list_fields)
% % 
% %     %
% %     min_max_times(1) = min([min_max_times(1), min(smartMoorlvl1.(list_fields{i}).dtime)]);
% %     min_max_times(2) = max([min_max_times(2), max(smartMoorlvl1.(list_fields{i}).dtime)]);
% % 
% % end
% % 
% % % Round time edges so that there is no fractional second
% % min_max_times(1) = datetime(min_max_times(1).Year, min_max_times(1).Month, min_max_times(1).Day, ...
% %                             min_max_times(1).Hour, min_max_times(1).Minute, ceil(min_max_times(1).Second), ...
% %                             'TimeZone', min_max_times.TimeZone);
% % min_max_times(2) = datetime(min_max_times(2).Year, min_max_times(2).Month, min_max_times(2).Day, ...
% %                             min_max_times(2).Hour, min_max_times(2).Minute, floor(min_max_times(2).Second), ...
% %                             'TimeZone', min_max_times.TimeZone);


% --------------------------------------------------------
% Min/max time chosen mannually (and removing some data at the edges)

% Starting on 06/17 17:45:00 should remove gaps with changing
% spotter mode for those deployed on the first day

%
min_max_times = [datetime(2022, 06, 17, 17, 45, 00), ...
                 datetime(2022, 07, 20, 06, 00, 00)];

min_max_times.TimeZone = 'America/Los_Angeles';


%% Create time grid

%
dt_grid = seconds(0.5);
%
dtime = min_max_times(1) : dt_grid : min_max_times(2);


%% Loop over data and grid pressure
%
% Takes less than 1 minute

%
gap_TH = 5;    % gap threshold in seconds


% Loop over Smart Moorings
tic
for i1 = 1:length(list_fields)

    %
    disp(' ')
    disp(' ')
    %
    disp(['--- Start to grid pressure from smart mooring ' list_fields{i1} ' ---'])

    % Initialize logical array to tell where we want
    % (and don't want to interpolate pressure)
    linterp_ok = true(1, length(dtime));

    % Find whether there are data points before or after the beginning
    % and end of the time grid
    lbefore_data = (dtime - smartMoorlvl1.(list_fields{i1}).dtime(1)) < 0;
    lafter_data = (dtime - smartMoorlvl1.(list_fields{i1}).dtime(end)) > 0;

    % Remove edges from interpolation if there are no data points
    linterp_ok(lbefore_data) = false;
    linterp_ok(lafter_data) = false;

    % Find where the are gaps longer than the threshold value
    lgaps = (diff(smartMoorlvl1.(list_fields{i1}).dtime) > seconds(gap_TH));

    %
    if any(lgaps)

        % Indices of gaps
        ind_gaps = find(lgaps);
        %
        Nind_gaps = length(ind_gaps);

        % Initialize an array for the time limits where there
        % are continuous data
        time_edges = NaT((Nind_gaps + 1), 2);
        time_edges.TimeZone = 'America/Los_Angeles';
        % Add the beginning and end
        time_edges(1, 1) = smartMoorlvl1.(list_fields{i1}).dtime(1);
        time_edges(end, 2) = smartMoorlvl1.(list_fields{i1}).dtime(end);

        % Add edges of time segments
        for i2 = 1:Nind_gaps

            %
            time_edges(i2, 2) = smartMoorlvl1.(list_fields{i1}).dtime(ind_gaps(i2));
            time_edges(i2+1, 1) = smartMoorlvl1.(list_fields{i1}).dtime(ind_gaps(i2)+1);

        end

        % Set to false the time grid points that don't have nearby
        % data points
        for i2 = 1:Nind_gaps
            
            %
            lingap_aux = (dtime > time_edges(i2, 2)) & ...
                         (dtime < time_edges((i2+1), 1));

            %
            linterp_ok(lingap_aux) = false;

        end
    end


    %% Time-grid pressure

    %
    smartMoorlvl1.(list_fields{i1}).pressureinterp = NaN(1, length(dtime));

    % Interpolate where distance from data points is sufficiently short
    smartMoorlvl1.(list_fields{i1}).pressureinterp(linterp_ok) = ...
                        interp1(smartMoorlvl1.(list_fields{i1}).dtime, ...
                                smartMoorlvl1.(list_fields{i1}).pressure, ...
                                dtime(linterp_ok));

    %
    disp(['--- Done gridding pressure from smart mooring ' list_fields{i1} ' ---'])


end

toc



%% Plot time-gridded data

%
figure
    hold on
    %
    for i = 1:length(list_fields)
        plot(dtime, ...
             smartMoorlvl1.(list_fields{i}).pressureinterp + 2.*(i-1), '.-')
    end

    %
    ylabel('Water pressure [dbar]', 'Interpreter', 'Latex', 'FontSize', 22)
    %
    grid on
    set(gca, 'FontSize', 16)
    %
    set(gcf, 'units', 'normalized')
    set(gcf, 'Position', [0.2, 0.2, 0.6, 0.4])


%% Now put the interpolated pressure data in individual structures
% and save each of them separately


%
for i = 1:length(list_fields)

    %
    spotsmartL1.mooringID = smartMoorlvl1.(list_fields{i}).mooringID;
    spotsmartL1.SN = smartMoorlvl1.(list_fields{i}).SN;
    %
    spotsmartL1.latitude = smartMoorlvl1.(list_fields{i}).latitude;
    spotsmartL1.longitude = smartMoorlvl1.(list_fields{i}).longitude;
    %
    spotsmartL1.zhab = smartMoorlvl1.(list_fields{i}).zhab;

    %
    spotsmartL1.dt = dt_grid;
    spotsmartL1.timelims = min_max_times;
    spotsmartL1.gapTH = seconds(gap_TH);
    %
    spotsmartL1.dtime = dtime;
    %
    spotsmartL1.pressure = smartMoorlvl1.(list_fields{i}).pressureinterp;
    
    % Turn all vectors into column vectors (so
    % that Matlab displays structure variable quickly)
    list_somefields = fieldnames(spotsmartL1);
    %
    for i2 = 1:length(list_somefields)
        %
        if isvector(spotsmartL1.(list_somefields{i2})) && ...
           ~isstruct(spotsmartL1.(list_somefields{i2})) && ...
           ~ischar(spotsmartL1.(list_somefields{i2}))
            %
            spotsmartL1.(list_somefields{i2}) = spotsmartL1.(list_somefields{i2})(:);
        end
    end

    %
    time_dataproc = datestr(datetime('now', 'TimeZone', 'America/Los_Angeles'), 'yyyy/mm/dd HH:MM:SS');
    %
    spotsmartL1.README = ['Level 1 (gridded) smart mooring data from ROXSI ' ...
                        '2022. Data gridded by script ' mfilename() '.m on ' ...
                        time_dataproc ' (PDT). Data gaps shorter or equal than ' ...
                        'the gap threshold (gapTH) are interpolated over. ' ...
                        'Pressure is the water pressure. ' ...
                        'Longitude and latitude (for the pressure measurement) ' ...
                        'were computed from the mean watch circle obtained ' ...
                        'from Spotter coordinates over the whole deployment. ' ...
                        'zhab is the the height in meters of the pressure ' ...
                        'sensor above the bottom.'];

    %
    disp(' ')
    disp(' ')
    %
    disp('---- Saving smart mooring level 1 data ---- ')
    % Save Level 1 gridded data
% %     save(fullfile(repo_dirpath(), ['smart_mooring_' spotsmartL1.mooringID '_' spotsmartL1.SN '_L1_gridded.mat']), 'spotsmartL1');
    save(fullfile(dir_output, ['roxsi_smartmooring_L1_' spotsmartL1.mooringID '_' spotsmartL1.SN '_gridded.mat']), 'spotsmartL1');
    %
    disp(['---- Done with level 1 time-gridding pressure from smart mooring ' spotsmartL1.mooringID ' - SN ' spotsmartL1.SN ' ----'])


    %
    clear spotsmartL1
end
