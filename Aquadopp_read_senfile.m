function senAQDP = Aquadopp_read_senfile(file_name, list_vars)
%% senAQDP = AQUADOPP_READ_SENFILE(file_name, list_vars)
%
%   inputs
%       - file_name: file name of the Aquadopp *.sen file.
%       - list_vars (optional): 
%
%   outputs
%       - senAQDP: structure variable with fields containing the
%                  variables read from the *.sen file.
%
%
%
%
%
% Olavo Badaro Marques, 04/Aug/2022.


%%

%
if ~exist('list_vars', 'var')
    %
    list_vars = ["time", "pressure"];
end

%
Nvariables = length(list_vars);


%% 



% From *.hdr (header file):
%
% month day year hour min sec / error code and status code (unread) /
% voltage / soundspeed / heading, pitch, roll /
% pressure / temperature / analog input 1 and 2 (unread )

%
read_format_cell = {'%*s%*s%*s%*s%*s%*s', '%*f%*f', '%*f', '%*f', '%*f', '%*f', '%*f', '%*f', '%*f', '%*f%*f'};
% % read_format_cell = {'%*d%*d%*d%*d%*d%*d', '%*f%*f', '%*f', '%*f', '%*f', '%*f', '%*f', '%*f', '%*f', '%*f%*f'};
vars_cell = {'time', 'codes', 'voltage', 'soundspeed', 'heading', 'pitch', 'roll', 'pressure', 'temperature', 'analoginputs'};
%
varsncols = [6, 2, 1, 1, 1, 1, 1, 1, 1, 2];


%% I NEED TO SORT list_vars S0 THAT IT APPEARS IN THE SAME ORDER AS vars_cel!!!


%% For the variables you want to read, edit the string
% format such that it removes the string * of the variables
% that will be read

%
for i = 1:Nvariables
    %
    ind_match_var = find(strcmp(vars_cell, list_vars(i)));

    %
    read_format_cell(ind_match_var) = erase(read_format_cell(ind_match_var), '*');

end

%
str_read_data = [read_format_cell{:}];


%% Read just the first line to check whether the number
% of variables matches what this function expects

%
fid_aux = fopen(file_name);

%
Ntotalcols = 17;

%
first_line = fgetl(fid_aux);

%
inds_spaces = strfind(first_line, ' ');

% Remove spaces that appear continuously
diff_inds_spaces = diff(inds_spaces);

%
keep_space = false(1, length(inds_spaces));
%
keep_space([true, diff_inds_spaces~=1]) = true;

%
if (length(find(keep_space))+1) ~= Ntotalcols
    %
    error(['The number of variables (columns) in ' file_name ' is different than what the data reading function expects.'])
end


%%

% Close the file and open it again, because otherwise
% it skips the first line that I read above with fgetl
%
fclose(fid_aux);
%
fid_aux = fopen(file_name);

%
data_Aquadopp = textscan(fid_aux, str_read_data);


%% Match the variables that have been read with the columns where they are located

% Cell array because I'm not treating all columns as independent
indcolsread = NaN(1, Nvariables);

%
for i = 1:Nvariables
    %
    indcolsread(i) = varsncols(strcmp(vars_cell, list_vars(i)));    
end

%
indcolsread = cumsum(indcolsread);


%% Organize output (convert date to datetime)

%
for i = 1:Nvariables
    %
    if strcmp(list_vars(i), "time")

        %
        Nobs = length(data_Aquadopp{1});

% %         %
% %         senAQDP.time = repmat(" ", Nobs, 1);
% %         keyboard
% %         %
% %         for i2 = 1:Nobs
% % 
% %             %
% %             senAQDP.time(i2) = [data_Aquadopp{3}{i2} '/' data_Aquadopp{1}{i2} '/' data_Aquadopp{2}{i2} ' ' ...
% %                                 data_Aquadopp{4}{i2} ':' data_Aquadopp{5}{i2} ':' data_Aquadopp{6}{i2}];
% % 
% %         end
        %
        cell_dates_cat = strcat(data_Aquadopp{3}, repmat({'/'}, Nobs, 1), data_Aquadopp{1}, repmat({'/'}, Nobs, 1), data_Aquadopp{2}, ...
                                repmat({' '}, Nobs, 1), ...
                                data_Aquadopp{4}, repmat({':'}, Nobs, 1), data_Aquadopp{5}, repmat({':'}, Nobs, 1), data_Aquadopp{6});
        senAQDP.time = string(cell_dates_cat);
        
    %
    else

        %
        senAQDP.(list_vars(i)) = data_Aquadopp{indcolsread(i)};

    end
end






