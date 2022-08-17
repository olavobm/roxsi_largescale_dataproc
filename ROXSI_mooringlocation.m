function mooringtableout = ROXSI_mooringlocation(mooringID, strinstrument, strmooringline)
%% mooringtableout = ROXSI_MOORINGLOCATION(mooringID, strinstrument, strmooringline)
%
%   inputs
%       - mooringID (optional): ID (or IDs) of the desired mooring
%                               locations. This can be both the ID
%                               with or without the lower case letters
%                               specifying the instruments (but don't
%                               mix and match formats).
%       - strinstrument (optional): string to restrict which mooring
%                                   instruments to load (note that ADCP
%                                   moorings with T-chains have two
%                                   locations for the same mooring ID).
%       - strmooringline (optional):
%
%
%   outputs
%       - mooringtableout: a Matlab table with mooring locations.
%
%
% If no inputs are given, then the output includes all mooring data.
%
%
% Still work in progress.


%% Load mooring table that must be saved in the
% same directory as this function

%
full_filename = mfilename('fullpath');

%
indslash = strfind(full_filename, filesep);

%
dir_table = full_filename(1:indslash(end));

%
load(fullfile(dir_table, 'ROXSI2022_mooringtable.mat'), 'mooringtable')


%%

%
if nargin==0

    %
    mooringtableout = mooringtable;

% ------------------------
% Could include an ifelse option without
% mooringID, but with the other inputs
% ------------------------

%
else

    % ------------------------------------------------
    % Convert mooringID to cell array so that it can
    % be compared with the mooringID in the table
    % without looking at lower case letters at the end

    % If it's a character array, convert to cell array
    if ischar(mooringID)
        mooringID = {mooringID};
    end

    % If mooring ID is a string array, convert to cell array
    if isstring(mooringID)
        mooringID = convertStringsToChars(mooringID);
        mooringID = {mooringID};
    end

    % ------------------------------------------------
    %
    if exist('strinstrument', 'var') && ~isempty(strinstrument)
        %
        lsubinstrument = true;

        % If it's a character array, convert to cell array
        if ischar(strinstrument)
            strinstrument = {strinstrument};
        end
        % If mooring ID is a string array, convert to cell array
        if isstring(strinstrument)
            strinstrument = convertStringsToChars(strinstrument);
            strinstrument = {strinstrument};
        end
        
    %
    else
        lsubinstrument = false;
        
    end

% %     %
% %     if exist('strinstrument', 'var') && isempty(strinstrument)
% %         %
% %         lsubline = true;
% %     %
% %     else
% %         lsubline = false;
% %         
% %     end


    % ------------------------------------------------
    % Pre-allocate logical variable to subset the mooring table
    loutput = false(size(mooringtable, 1), 1);

    % Loop over moorings
    for i = 1:length(mooringID)
        %
        lmatch_aux = strncmp(mooringtable.mooringID, mooringID{i}, length(mooringID{i}));

        %
        if lsubinstrument
            %
            lmatch_instrument_aux = strncmp(mooringtable.instrument, strinstrument{i}, length(strinstrument{i}));

            %
            lmatch_aux = lmatch_aux & lmatch_instrument_aux;
        end

        %
        if any(lmatch_aux)
            %
            loutput(lmatch_aux) = true;
        end


% %         %
% %         if exist('strmooringline', 'var')
% % 
% %         end

    end

    % ------------------------------------------------
    mooringtableout = mooringtable(loutput, :);


end

