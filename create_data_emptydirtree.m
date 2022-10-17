%% Create directory tree with no files for temporary data.
% It starts by removing an existing dir tree. !!!BE CAREFUL WITH THIS!!!

clear
close all


%%

%
dir_parent_location = '/home/omarques/Documents/obm_ROXSI/';
%
dir_parent_name = 'obm_DataLocal';

%
dir_data_parent = fullfile(dir_parent_location, dir_parent_name);


%% Remove (or rename) an existing dir tree that
% matches with what this script creates


if isfolder(dir_data_parent)

    % % %
    % % command_remove = ['rm -rf ' dir_data_parent];
    % % system(command_remove)

    % Or change name
    %
    time_dataproc = datetime('now', 'TimeZone', 'Local');
    time_dataproc_char = datestr(time_dataproc, 'yyyy_mm_dd_HH_MM');
    %
    command_rename_old = ['mv ' dir_data_parent ' ' dir_data_parent '_' time_dataproc_char];
    system(command_rename_old)

end


%% Create the data directory tree

%
mkdir(dir_data_parent)

%
mkdir(fullfile(dir_data_parent, 'Level1_Data'))
mkdir(fullfile(dir_data_parent, 'Level2_Data'))


%
mkdir(fullfile(dir_data_parent, 'Level1_Data', 'Spotter_Level1'))
mkdir(fullfile(dir_data_parent, 'Level1_Data', 'Spotter_Smart_Level1'))

% --------------------
%
mkdir(fullfile(dir_data_parent, 'Level1_Data', 'Spotter_Smart_Level1', 'not_gridded'))
mkdir(fullfile(dir_data_parent, 'Level1_Data', 'Spotter_Smart_Level1', 'gridded'))

%
mkdir(fullfile(dir_data_parent, 'Level1_Data', 'Spotter_Smart_Level1', 'not_gridded', 'figs_QC'))
mkdir(fullfile(dir_data_parent, 'Level1_Data', 'Spotter_Smart_Level1', 'gridded', 'figs_QC'))


% --------------------
mkdir(fullfile(dir_data_parent, 'Level1_Data', 'Spotter_Level1', 'figs_QC'))



% --------------------------------------------------------

%
mkdir(fullfile(dir_data_parent, 'Level2_Data', 'Spotter_Level2'))
mkdir(fullfile(dir_data_parent, 'Level2_Data', 'Spotter_Smart_Level2'))


%
mkdir(fullfile(dir_data_parent, 'Level2_Data', 'Spotter_Level2', 'figs_QC'))
mkdir(fullfile(dir_data_parent, 'Level2_Data', 'Spotter_Smart_Level2', 'figs_QC'))
