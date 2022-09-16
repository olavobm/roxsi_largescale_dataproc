%% Loop over Spotters of interest and convert
% data from csv files to *.mat files.

clear
close all

%%

%
dir_parentdata = fullfile(data_dirpath(), 'RAW');

% Output directories
% Same directory for Spotters and SmartSpotters
dir_output = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';


%%

%
dir_Spotters = fullfile(dir_parentdata, 'Spotters');
%
list_Spotters = {'B01_spot1150', 'B01_spot1158', ...
                 'B03_spot1152', 'B05_spot1153', ...
                 'X01_spot1151', 'X03_spot1157', 'X04_spot1155'};

%
dir_SmartSpotters = fullfile(dir_parentdata, 'Spotters_Smart');
%
list_SmartSpotters = {'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
                      'E07_spot1855', 'E07_spot1857', ...
                      'E08_spot1852', ...
                      'E09_spot1850', 'E09_spot1856', ...
                      'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};


%%

%
timezone_UTCtoLocal = -7;


%% Loop over Spotters and convert data from csv to *.mat files

% Standard Spotters
for i = 1:length(list_Spotters)
    %
    disp(' '), disp(' ')
    %
    s = Spotter_parsed_to_mat(fullfile(dir_Spotters, 'SDcards', list_Spotters{i}, 'parsed'), ...
                              timezone_UTCtoLocal);
    
    % Save output data
    save('-v7.3', fullfile(dir_output, [list_Spotters{i} '.mat']), 's')
    %
    clear s
end

% Smart Spotters
for i = 1:length(list_SmartSpotters)
    %
    disp(' '), disp(' ')
    %
    s = Spotter_parsed_to_mat(fullfile(dir_SmartSpotters, 'SDcards', list_SmartSpotters{i}, 'parsed'), ...
                              timezone_UTCtoLocal);
    
    % Save output data
    save('-v7.3', fullfile(dir_output, [list_SmartSpotters{i} '.mat']), 's')
    %
    clear s
end


