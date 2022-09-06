function loginfo = RBR_extract_loginfo(log_filename, day_info, SN)
%% loginfo = RBR_EXTRACT_LOGINFO(log_filename, day_info, SN)
%
%   inputs
%       - log_filename: file name (with extension).
%       - day_info: in the format yyyy-mm-dd.
%       - SN: serial number of the sensor of interest
%
%   outputs
%       - loginfo:
%
%
%
%
%
%
%

%% 

%
fid = fopen(log_filename);


%
lfound_SN = false;

%
loginfo = [];

%
% while 1 && (~lfound_SN)    % this assumes that after finding one instance
                             % of instrument SN and reaching the following instrument with a different
                             % serial number, there is no more informatio about the first instrument

% Skip a few lines
disp(' ')

% Read the full file
while 1
    
    % Read new line and terminate loop at the end of file
    tline = fgetl(fid);
    if ~ischar(tline), break, end

    
    % Only look for information if it's on the day
    % of the programming
    if strncmp(tline, day_info, length(day_info))

        %
        if contains(tline, 'RESPONSE: id model = RBRsolo,') && contains(tline, SN)

            % Print message to the screen that log information
            % has been found for the specific RBR instrument
            disp(['--- log info from RBR SN ' SN ' has been found ---'])

            % If the instrument SN has NOT been found yet
            if ~lfound_SN
                %
                lfound_SN = true;

            % If the instrument SN has already been found, then there 
            % is a second instance with information on instrument SN.
            else
                % Throw a warning and break the loop/end the function
                warning(['A second instance of RBR SN ' SN ' was found!!!'])

                %
                break
            end

% %             lfound_SN = true;
% %             keyboard

            %
            while lfound_SN


                %
                if contains(tline, 'RESPONSE: now =') || ...
                   contains(tline, 'Apparent time difference between logger and host:')

% % %                         % Output line information
% % %                         tline
% %                         keyboard

                    % First iteration/line
                    if isempty(loginfo)
                        %
                        loginfo = string(tline);
                    % Append other lines in a string array
                    else
                        loginfo = [loginfo; string(tline)];
                    end
                    
                end
                
                %
                if contains(tline, 'RESPONSE: id model = RBRsolo,')
                    if ~contains(tline, SN)
%                         %
%                         lfound_SN = false;
                        %
                        disp(['--- Done with getting log info from RBR SN ' SN ' ---'])

                        %
                        break
                    end
                end
                
                % Read the next line for the next iteration
                % of the while loop
                tline = fgetl(fid);
                %
                % (note there is a potential, but unlikely problem
                % here). When . 

                % If the instrument SN is the last one in the
                % file, then the inner loop needs to be told
                % that then end of the file has been reached
                if ~ischar(tline), break, end
            end

        end
    end


end

%
fclose(fid);

%
if isempty(loginfo)
    %
    warning(['--- NO log on RBR SN ' SN ' was found!!! ---'])
end


% Skip a few lines
disp(' ')



