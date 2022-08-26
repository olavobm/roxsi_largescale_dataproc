function dataAll = Spotter_readmulti_SMD(dir_SDcard, list_files)
%% dataAll = SPOTTER_READMULTI_SMD(dir_SDcard, list_files)
%
%   inputs
%       - dir_SDcard: directory where the content of the Spotter SD card
%                     has been saved.
%       - list_files (optional): a cell or string array with the filenames
%                                in the order that they will be read. If
%                                this input is not given, then read all
%                                SMD files in the folder.
%
%   outputs
%       - dataAll: data structure with the data in the SMD files.
%
%
% SPOTTER_READMULTI_SMD.m is a higher level function reading
% pressure data from *_SMD.CSV files recorded in the Spotter
% (Smart Mooring) SD card. The data are concatenated together
% in the order of the output of dir.m (alphabetical/ascending
% order).
%
% NOTE(!!!) the time output is in UTC time (at least
% when I wrote this function).
%
% Olavo Badaro Marques


%% Get the list of all *_SMD.CSV files in the target directory

%
list_files = dir(fullfile(dir_SDcard, "*_SMD.CSV"));
%
Nfiles = length(list_files);


%% Loop over files and load all the data

%
for i = 1:Nfiles

    %
    dataAll.eachfile(i) = Spotter_read_SMD(fullfile(list_files(i).folder, list_files(i).name));

end


%% Concatenate all the data together

%
dataAll.allfiles.Nfiles = Nfiles;

%
list_fields = fieldnames(dataAll.eachfile(1));

% Remove nobservations, because this will
% be a vector with different dimension
% % list_fields = list_fields(~strcmp(list_fields, "nobservations"));
list_fields = setdiff(list_fields, {'filename', 'nobservations'}, 'stable');

% Preallocate variables
dataAll.allfiles.nobservations = NaN(1, Nfiles);
dataAll.allfiles.filename = strings(1, Nfiles);
% Create variables as empty arrays
for i = 1:length(list_fields)
    %
    dataAll.allfiles.(list_fields{i}) = [];
end

%
for i1 = 1:Nfiles

    %
    ndata_points = dataAll.eachfile(i1).nobservations;

    %
    dataAll.allfiles.nobservations(i1) = ndata_points;
    dataAll.allfiles.filename = dataAll.eachfile(i1).filename;

    %
    for i2 = 1:length(list_fields)
        %
        dataAll.allfiles.(list_fields{i2}) = [dataAll.allfiles.(list_fields{i2}); ...
                                              dataAll.eachfile(i1).(list_fields{i2})];
    end


end





