%% Convert a *.csv file with mooring locations into a matlab table.

clear
close all

%%

%
full_filename = mfilename('fullpath');

%
indslash = strfind(full_filename, filesep);

%
dir_table = full_filename(1:indslash(end));

%
table_filename = 'ROXSI2022_mooring_locations.csv';


%% Parse columns into variables

%
fileID_aux = fopen(fullfile(dir_table, table_filename));

%
data_aux = textscan(fileID_aux, '%s %s %f%f %f%f %s %s', 'Headerlines', 1, 'Delimiter', ',');
        
%
mooringID = convertCharsToStrings(data_aux{1});

%
roxsiarray = convertCharsToStrings(data_aux{2});

%
planned_latitude = data_aux{3};
planned_longitude = data_aux{4};

%
latitude = data_aux{5};
longitude = data_aux{6};

%
instrument = convertCharsToStrings(data_aux{7});

%
clear data_aux


%% It turns out that the "actual"/GPS location of mooring
% C03ap is awardly close to C02p. Based Olavo swimming next
% to the floats with the DiveJet, and the procedute of
% deploying an ADCP from the boat, Olavo thinks the GPS
% location is wrong. Thus, I will replace the "actual" coordinate
% with the planned value.

% % %
% % indmatchC03 = find(strcmp(mooringID, "C03ap"));
% % %
% % if length(indmatchC03)~=1;    error('unexpected behavior');    end
% % 
% % %
% % latitude(indmatchC03) = planned_latitude(indmatchC03);
% % longitude(indmatchC03) = planned_longitude(indmatchC03);


%%

% Replace flag (999999) with NaN
%
flag_nopoint = 9999999999;

%
lmatchflag = latitude == flag_nopoint;

%
latitude(lmatchflag) = NaN;
longitude(lmatchflag) = NaN;
%
planned_latitude(lmatchflag) = NaN;
planned_longitude(lmatchflag) = NaN;


%% Create table and save in the directory dir_table

%
mooringtable = table(mooringID, roxsiarray, instrument, ...
                     latitude, longitude, ...
                     planned_latitude, planned_longitude);


%
save(fullfile(dir_table, 'ROXSI2022_mooringtable.mat'), 'mooringtable');