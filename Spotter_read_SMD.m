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
smd_data.filename = convertCharsToStrings([file_name_noextension '.CSV']);


%% (Trying) to read with a higher level function

%
loaded_SmartMooring = readmatrix(smdfile, 'NumHeaderLines', 1);

%
smd_data.filename = convertCharsToStrings([file_name_noextension '.CSV']);
%
smd_data.nobservations = size(loaded_SmartMooring, 1);
%
smd_data.unixEpoch = loaded_SmartMooring(:, 1);
smd_data.link = loaded_SmartMooring(:, 2);
smd_data.pressure = loaded_SmartMooring(:, 6);

%
smd_data.dtime = 719529 + (smd_data.unixEpoch./86400);


%% Due to a few instancies that lead to incosistent formatting
% of the SMD file, it's more rigorous to use fgetl
%
% -- but then it takes A LOT LONGER!!!
% about ~6 min per Smart Mooring with full data (> 10 times
% than using a high level function).

% % % -------------------------------------------------
% % % Get total number of lines
% % [~, wc_output] = system(['wc -l ' smdfile]);
% % %
% % ind_spaces = strfind(wc_output, ' ');
% % %
% % nlines = str2double(wc_output(3:(ind_spaces(end)-1)));
% % %
% % nlines_observations = nlines - 1;
% % 
% % % -------------------------------------------------
% % % Preallocate
% % %
% % smd_data.nobservations = nlines_observations;
% % %
% % smd_data.unixEpoch = NaN(nlines_observations, 1);
% % smd_data.link = NaN(nlines_observations, 1);
% % smd_data.logtype = strings(nlines_observations, 1);
% % %
% % smd_data.numbertimeweird = NaN(nlines_observations, 1);
% % smd_data.pressure = NaN(nlines_observations, 1);
% % 
% % % -------------------------------------------------
% % % Read the data
% % %
% % fid = fopen(smdfile);
% % 
% % % Read this to skip the 1-line header
% % tline = fgetl(fid);
% % 
% % % Counter to assign to preallocated output arrays
% % ind_count = 0;
% % 
% % %
% % while 1
% % 
% %     % Read new line and terminate loop at the end of file
% %     tline = fgetl(fid);
% %     if ~ischar(tline), break, end
% % 
% %     % Good to check after the loop that ind_count matches nlines_observations
% %     ind_count = ind_count + 1;
% % 
% %     % Find the commas (5 is the correct one, apart from the bad line)
% %     ind_commas = strfind(tline, ',');
% % 
% %     % Get unix epoch, link, and logtype variables
% %     smd_data.unixEpoch(ind_count) = str2double( tline(1:(ind_commas(1)-1)) );
% %     smd_data.link(ind_count) = str2double( tline((ind_commas(1)+1)) );
% %     %
% %     smd_data.logtype(ind_count) = tline((ind_commas(2)+1):(ind_commas(3)-1));
% % 
% %     % Depending on logtype, get other variables in different ways
% %     %
% %     if strcmp( smd_data.logtype(ind_count), 'DATA')
% %         % If DATA then read, normally
% %          smd_data.numbertimeweird(ind_count) = str2double( tline((ind_commas(4)+1):(ind_commas(5)-1)) );
% %          smd_data.pressure(ind_count) = str2double( tline((ind_commas(5)+1):end) );
% %     %
% %     elseif strcmp( smd_data.logtype(ind_count), 'BSYS')
% %         % Otherwise only read other variable. There (at least)
% %         % 2 instances of BSYS. In one, unixEpoch is zero and
% %         % numbertimeweird is also bogus. In another, they
% %         % are both fine, but it's the pressure that is not
% %         % available
% %          smd_data.numbertimeweird(ind_count) = str2double( tline((ind_commas(3)+1):(ind_commas(4)-1)) );
% %     %
% %     else
% %         error(['Unexpected string in ' smdfile])
% %     end
% % 
% % end
% % %
% % fclose(fid);
% % 
% % % Convert from Unix epoch to datenum
% % smd_data.dtime = 719529 + (smd_data.unixEpoch./86400);