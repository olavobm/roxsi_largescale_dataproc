%% Script that loads through all Signature1000 files
% and gets the min/max datenum (i.e. from Data.Burst_Time)
% from each file. This will make it more manageable to
% process all of the data.
%
% The table with time limits for each Signature1000
% are ordered in chronological order of the data
% (and NOT necessarily in terms of alphabetical order
% of file names).


clear
close all

%% Set data directory

%
% % data_parentdir = fullfile(data_dirpath(), 'RAW', 'Signature1000');
data_parentdir = fullfile('/project/CSIDE/ROXSI/LargeScale_Data_2022/', 'RAW', 'Signature1000');


%
list_Sig1000 = {'A01_103043', 'B10_103045', ...
                'B13_103046', 'B15_103056', ...
                'B17_101923', 'C01_102128', ...
                'X05_100231', 'X11_101941'};


%% Start a diary

%
repo_dirfull = pwd;

%
log_file_name = ['log_Sig1000_timelimits_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(repo_dirfull, log_file_name))
%
totalRunTime = tic;

% Display on the screen:
disp(' '), disp(' ')
disp('------------------------------ Getting time limits for all Signature1000 files ------------------------------')
%
list_Sig1000


%% Get all time limits

% Loop over Signature1000
for i1 = 1:length(list_Sig1000)

    % Add SN and directory path to output structure
    sig1000timelims.(list_Sig1000{i1}).SN = convertCharsToStrings(list_Sig1000{i1}(5:end));
    sig1000timelims.(list_Sig1000{i1}).dirdata = fullfile(data_parentdir, list_Sig1000{i1}, 'converted');

    % Get all *.mat files
    dirfiles = dir(fullfile(sig1000timelims.(list_Sig1000{i1}).dirdata, '*.mat'));

    %
    Nfiles = length(dirfiles);
    

    % Pre-allocate
    filenames = strings(Nfiles, 1);
    time_begin = NaN(Nfiles, 1);
    time_end = NaN(Nfiles, 1);
    %
    time_begin_string = strings(Nfiles, 1);
    time_end_string = strings(Nfiles, 1);

    % Loop over files
    for i2 = 1:5%Nfiles
    
        %
        filenames(i2) = dirfiles(i2).name;

        %
        data_aux = load(fullfile(dirfiles(i2).folder, dirfiles(i2).name));
    
        %
        time_begin(i2) = min(data_aux.Data.Burst_Time);
        time_end(i2) = max(data_aux.Data.Burst_Time);

        %
        time_begin_string(i2) = convertCharsToStrings(datestr(time_begin(i2)));
        time_end_string(i2) = convertCharsToStrings(datestr(time_end(i2)));
    end

    %
    sig1000timelims.(list_Sig1000{i1}).Nfiles = Nfiles;

    %
    sig1000timelims.(list_Sig1000{i1}).timelims = ...
                                 table(filenames, time_begin, time_end, ...
                                       time_begin_string, time_end_string);
    %
    disp(['Done with Signature1000 ' list_Sig1000{i1} ' (' num2str(i1) ' out of ' num2str(length(list_Sig1000)) ')'])
    
end


%% Sort tables in terms of ascending time_begin

% Loop over Signature1000
for i1 = 1:length(list_Sig1000)

    %
    [~, ind_sort] = sort(sig1000timelims.(list_Sig1000{i1}).timelims.time_begin);

    %
    sig1000timelims.(list_Sig1000{i1}).timelims = sig1000timelims.(list_Sig1000{i1}).timelims(ind_sort, :);

end


%% Now save the structure in the current folder

%
% % save(fullfile(dir_output_data_L1, [str_filename '.mat']), 'sig1000timelims', '-v7.3')

%
save('Signature1000_time_limits.mat', 'sig1000timelims')


%% Close diary

%
disp('------------------------------ Done with getting time limits ------------------------------')

%
disp(' '), disp(' ')
disp('*** The total time to run this script was:')
%
toc(totalRunTime)

% Close the log file
diary('off');