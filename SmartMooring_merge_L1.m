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
ind_match_buoy = 1;
ind_match_pressure = 1;

%
for i = 1:length(list_buoydata_files)
end

%
for i = 1:length(list_pressure_files)
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



% % %
% % data_buoy_aux
% % data_pres_aux

%% Copy buoy data

%
mergedL1 = filefullname_pres;

% Remove stuff???


%% Copy pressure data

%
mergedL1.pressuredata.latitude = 
mergedL1.pressuredata.longitude = 

%
mergedL1.pressuredata.zhab = 
mergedL1.pressuredata.dt = 

%
mergedL1.pressuredata.dtime = 
mergedL1.pressuredata.pressure = 


%% Final stuff

