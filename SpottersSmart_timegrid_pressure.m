%% Time grid pressure data from all smart moorings
% in ROXSI 2022 to the same time grid

%
clear
close all


%% Data directory

%
dir_data_L1 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/code_proc';


%% Get all L1 file names in folder

%
dir_all_L1 = dir(fullfile(dir_data_L1, '*_L1.mat'));
%
list_files = cell(1, length(dir_all_L1));

% Put file names in a cell array
for i = 1:length(list_files)
    %
    list_files{i} = dir_all_L1(i).name;
end


%% Loop over file names and load data

for i = 1:length(list_files)

    %
    data_loaded = load(fullfile(dir_data_L1, list_files{i}));

    %
    smartMoorlvl1.([list_files{i}(15:17) '_' list_files{i}(21:24)]) = data_loaded.spotsmart;

end

%
list_fields = fieldnames(smartMoorlvl1);


%% Get min/max date-times where observations are available

%
min_max_times = [min(smartMoorlvl1.(list_fields{1}).dtime), ...
                 max(smartMoorlvl1.(list_fields{1}).dtime)];

% Get min/max time
for i = 2:length(list_fields)

    %
    min_max_times(1) = min([min_max_times(1), min(smartMoorlvl1.(list_fields{i}).dtime)]);
    min_max_times(2) = max([min_max_times(2), max(smartMoorlvl1.(list_fields{i}).dtime)]);

end

% Round time edges so that there is no fractional second


% Or instead set something mannually


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
