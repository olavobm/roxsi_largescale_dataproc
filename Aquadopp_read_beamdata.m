function beamAquadopp = Aquadopp_read_beamdata(file_name, list_vars)
%% beamAquadopp = AQUADOPP_READ_BEAMDATA(file_name, list_vars)
%
%   inputs
%       - file_name: file name of the Aquadopp *.v* and *.c* files.
%       - list_vars (optional): 
%
%   outputs
%       - beamAquadopp: structure variable with fields containing the
%                       read variables.
%
%
%
%
%
% Olavo Badaro Marques, 12/Aug/2022.


%%

%
if ~exist('list_vars', 'var')
    %
    list_vars = ["v1", "v2", "v3", "a1", "a2", "a3"];
end

%
Nvariables = length(list_vars);


%% Check first line of one file to see how
% many columns (= number of bins) the file has

%
file_read_aux = [file_name '.' char(list_vars(1))];
%
fid_aux = fopen(file_read_aux);

%
first_line = fgetl(fid_aux);

% Remove potential spaces before and after the
% first and last numbers (respectively)
first_line = strtrim(first_line);

%
inds_spaces = strfind(first_line, ' ');

% Remove spaces that may appear continuously
diff_inds_spaces = diff(inds_spaces);
%
keep_space = false(1, length(inds_spaces));
%
keep_space([true, diff_inds_spaces~=1]) = true;

%
Nbins = length(find(keep_space)) + 1;

% Close the file (otherwise textscan below skips
% the line read with fgetl)
fclose(fid_aux);


%% Create string that specifies the format to read the data

%
str_read_data = repmat('%f', 1, Nbins);


%% Reads the data


% Loop over files/variables that will be read
for i1 = 1:Nvariables

    %
    file_read_aux = [file_name '.' char(list_vars(i1))];

    %
    fid_aux = fopen(file_read_aux);
    
    %
    data_Aquadopp = textscan(fid_aux, str_read_data);

    % Pre-allocate space
    beamAquadopp.(list_vars(i1)) = NaN(Nbins, length(data_Aquadopp{1}));

    %
    for i2 = 1:Nbins
        %
        beamAquadopp.(list_vars(i1))(i2, :) = data_Aquadopp{i2};
    end


end

