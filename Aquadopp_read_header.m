function headerAquadopp = Aquadopp_read_header(file_name)
%% headerAquadopp = AQUADOPP_READ_HEADER(file_name)
%
%   inputs
%       - file_name: file name of the Aquadopp *.hdr file (with extension).
%
%   outputs
%       - headerAquadopp: structure variable.
%
%
% AQUADOPP_READ_HEADER.m reads an Aquadopp *.hdr file and extract
% some of its information.
%
% Olavo Badaro Marques.



%% The sections of the header that are extracted separately

%
line_overview_datapoints = 'Number of measurements';

%
line_UserSetup = 'User setup';

%
line_hardwareConfig = 'Hardware configuration';
%
line_headConfig = 'Head configuration';


%
line_cellcenters = 'Current profile cell center distance from head (m)';


%% Open the header file and extract the information

%
fidHeaderFile = fopen(file_name);

%
while 1
    
    %
    line_aux = fgetl(fidHeaderFile);

    % Breaks while loop when reaches the end of file (where
    % fgetl outputs the number -1)
    if ~ischar(line_aux)
        break
    end

    %% ------------------------------------------------------------
    % Overview *.prf file
    %
    if strncmp(line_aux, line_overview_datapoints, length(line_overview_datapoints))

        %
        while length(line_aux)>1

            % ---------------------------------
            % Save new line in output
            if ~exist('data_Overview', 'var')
                data_Overview = sprintf([line_aux '\n']);
            else
                %
                line_append_aux = regexprep(line_aux, '%', '%%');
                %
                data_Overview = [data_Overview, sprintf([line_append_aux '\n'])];
            end

            % ---------------------------------
            % Next line
            line_aux = fgetl(fidHeaderFile);
            
        end

        %
        headerAquadopp.header.Overview = data_Overview;

    end


    %% ------------------------------------------------------------
    % User setup section
    %
    if strncmp(line_aux, line_UserSetup, length(line_UserSetup))

        %
        while length(line_aux)>1

            % ---------------------------------
            % Save new line in output
            if ~exist('data_UserSetup', 'var')
                data_UserSetup = sprintf([line_aux '\n']);
            else
                %
                line_append_aux = regexprep(line_aux, '%', '%%');
                %
                data_UserSetup = [data_UserSetup, sprintf([line_append_aux '\n'])];
            end
            
            % ---------------------------------
            % When reaching Profile interval, extract to a number (here we assume this is in seconds)
            if strncmp(line_aux, 'Profile interval', length('Profile interval'))
                %
                ind_spaces = strfind(line_aux, ' ');

                %
                info_aux = line_aux((ind_spaces(end-1)+1):(ind_spaces(end)-1));
                %
                headerAquadopp.profileInterval = str2double(info_aux);

            end

            % ---------------------------------
            % When reaching Cell size, extract to a number (convert from cm to m)
            if strncmp(line_aux, 'Cell size', length('Cell size'))
                %
                ind_spaces = strfind(line_aux, ' ');

                %
                info_aux = line_aux((ind_spaces(end-1)+1):(ind_spaces(end)-1));
                %
                headerAquadopp.binsize = str2double(info_aux)./100;

            end

            % ---------------------------------
            % When reaching Average interval, extract to a number (here we assume this is in seconds)
            if strncmp(line_aux, 'Average interval', length('Average interval'))
                %
                ind_spaces = strfind(line_aux, ' ');

                %
                info_aux = line_aux((ind_spaces(end-1)+1):(ind_spaces(end)-1));
                %
                headerAquadopp.averageInterval = str2double(info_aux);

            end

            % ---------------------------------
            % When reaching number of cells, extract to a number
            if strncmp(line_aux, 'Number of cells', length('Number of cells'))

                %
                info_aux = line_aux((length('Number of cells') + 1):length(line_aux));
                %
                headerAquadopp.numberofcells = str2double(strip(info_aux));

            end

            % ---------------------------------
            % Next line
            line_aux = fgetl(fidHeaderFile);
            
        end

        %
        headerAquadopp.header.UserSetup = data_UserSetup;
    end


    %% ------------------------------------------------------------
    % Hardware configuration
    if strncmp(line_aux, line_hardwareConfig, length(line_hardwareConfig))
        %
        while length(line_aux)>1

            % ---------------------------------
            % Save new line in output
            if ~exist('hardware_Config', 'var')
                hardware_Config = sprintf([line_aux '\n']);
            else
                %
                line_append_aux = regexprep(line_aux, '%', '%%');
                %
                hardware_Config = [hardware_Config, sprintf([line_append_aux '\n'])];
            end

            % ---------------------------------
            % Next line
            line_aux = fgetl(fidHeaderFile);
            
        end

        %
        headerAquadopp.header.hardwareConfiguration = hardware_Config;
    end


    %% ------------------------------------------------------------
    % Head configuration
    if strncmp(line_aux, line_headConfig, length(line_headConfig))

        %
        while length(line_aux)>1

            % ---------------------------------
            % Save new line in output
            if ~exist('head_Config', 'var')
                head_Config = sprintf([line_aux '\n']);
            else
                %
                line_append_aux = regexprep(line_aux, '%', '%%');
                %
                head_Config = [head_Config, sprintf([line_append_aux '\n'])];
            end

            % ---------------------------------
            % Next line
            line_aux = fgetl(fidHeaderFile);
            
        end

        %
        headerAquadopp.header.headConfiguration = head_Config;
    end


    %% ------------------------------------------------------------
    % When reaching center of cells section, extract the cell centers
    if strncmp(line_aux, line_cellcenters, length(line_cellcenters))

        % Skip the line with a bunch of dashes
        fgetl(fidHeaderFile);

        % Pre-allocate space in the output variable
        headerAquadopp.cellcenter = NaN(headerAquadopp.numberofcells, 1);

        %
        for i = 1:headerAquadopp.numberofcells

            %
            line_aux = fgetl(fidHeaderFile);

            %
            ind_spaces = strfind(line_aux, ' ');

            %
            headerAquadopp.cellcenter(i) = str2double(line_aux((ind_spaces(end)+1):end));

        end
    end


end



