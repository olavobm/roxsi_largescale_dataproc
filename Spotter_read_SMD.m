function smd_data = Spotter_read_SMD(smdfile)
%% smd_data = SPOTTER_READ_SMD(smdfile)
%
%   inputs
%       - smdfile: file name (with extension) of a *_SMD file.
%
%   outputs
%       - smd_data: structure array with at least the time and
%                   pressure in the file smdfile.
%
%
% SPOTTER_READ_SMD.m reads data from only one *_SMD.CSV file
% that is saved in the SD card by the Spotter/Smart Mooring.
%
% The time variable is read as it is, which (as I write this)
% is in unixEpoch units, referenced to UTC time.
%
% Olavo Badaro Marques

%
[dir_file, file_name_noextension] = fileparts(smdfile);

%
loaded_SmartMooring = readmatrix(smdfile, 'NumHeaderLines', 1);

%
smd_data.filename = convertCharsToStrings([file_name_noextension '.CSV']);
%
smd_data.nobservations = size(loaded_SmartMooring, 1);
%
smd_data.dataflag = loaded_SmartMooring(:, 2);
smd_data.unixEpoch = loaded_SmartMooring(:, 1);
smd_data.pressure = loaded_SmartMooring(:, 6);

%
smd_data.dtime = 719529 + (smd_data.unixEpoch./86400);