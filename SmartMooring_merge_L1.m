function mergedL1 = SmartMooring_merge_L1(dirbuoy, dirpressure, SN)
%% mergedL1 = SmartMooring_merge_L1(dirbuoy, dirpressure, SN)
%
%
%
%
%
%
%
%
%


%%

%
list_buoydata_files = dir(fullfile(dirbuoy, '*.mat'));
list_pressure_files = dir(fullfile(dirpressure, '*.mat'));


%% Find the files that match with the SN in the input

%
ind_match_buoy = 0;
ind_match_pressure = 0;

%
lloop_aux = true;
%
while lloop_aux && (ind_match_buoy < length(list_buoydata_files))
    %
    ind_match_buoy = ind_match_buoy + 1;
    %
    indfind_aux = strfind(list_buoydata_files(ind_match_buoy).filename, SN);
    
    %
    if ~isempty(indfind_aux)
        lloop_aux = false;
    end
end

%
lloop_aux = true;
%
while lloop_aux && (ind_match_pressure < length(list_pressure_files))
    %
    ind_match_pressure = ind_match_pressure + 1;
    %
    indfind_aux = strfind(list_pressure_files(ind_match_pressure).filename, SN);
    
    %
    if ~isempty(indfind_aux)
        lloop_aux = false;
    end
end



%%

%
filefullname_buoy = fullfile(list_buoydata_files(ind_match_buoy).folder, list_buoydata_files(ind_match_buoy).filename);
filefullname_pres = fullfile(list_pressure_files(ind_match_pressure).folder, list_pressure_files(ind_match_pressure).filename);

%
varBuoy = who('-file', filefullname_buoy);
varPres = who('-file', filefullname_pres);

%
data_buoy_aux = load(filefullname_buoy);
data_pres_aux = load(filefullname_pres);

%
data_buoy_aux = data_buoy_aux.(varBuoy{1});
data_pres_aux = data_pres_aux.(varPres{1});


%% Copy buoy data

%
mergedL1 = data_buoy_aux;

% Remove stuff???


%% Copy pressure data

%
list_fields = {'latitude', 'longitude', 'zhab', 'dt', ...
               'dtime', 'pressure'};

%
for i = 1:length(list_fields)
    mergedL1.pressuredata.(list_fields{i}) = data_pres_aux.(list_fields{i});
end




%% Final stuff

