function mooringwavestats = consistent_stats_roxsi2022(mooringID, freq_lims)
%%
%
%
%
% based on elevation or pressure-based elevation
%
%
%
%
%


%% Load mooring table

%
mooringtable = load(['/home/omarques/Documents/MATLAB' ...
                     '/roxsi_largescale_dataproc/ROXSI2022_mooringtable.mat']);
mooringtable = mooringtable.mooringtable;

%
% % mooringtable


%% Define directories

%
dirparent_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/';
% % Level1_Data

%
datapointers.instrumentdir.spotter = fullfile(dirparent_data, 'Level1_Data', 'Spotter_Level1');
datapointers.instrumentdir.smartspotter
%
datapointers.instrumentdir.SoloD
datapointers.instrumentdir.aquadopp
datapointers.instrumentdir.signature

% % %
% % dirbuoy = fullfile(dirparent_data, 'Spotter_Level1');
% % dirpressure = fullfile(dirparent_data, 'Spotter_Smart_Level1', 'gridded');


%%

%
list_SoloD = {''};
list_aquadopp = {''};
list_Signature = {''};
%
list_Spotter = {'B01', 'B03', 'B05', 'X01', 'X03', 'X04'};
list_SmartMooring = {'E01', 'E02', 'E05', 'E07', 'E08', 'E09', 'E10', 'E11', 'E13'};

%
datapointers.permooring.SoloD.mooringID = list_SoloD;
datapointers.permooring.SoloD.datadir = datapointers.instrumentdir.SoloD;
%
datapointers.permooring.aquadopp.mooringID = list_aquadopp;
datapointers.permooring.aquadopp.datadir = datapointers.instrumentdir.aquadopp;


%%

%
datapointers.permooring.Spotter.B01.filename = 'roxsi_spotter_L1_B01_1158.mat';
datapointers.permooring.Spotter.B03.filename = 'roxsi_spotter_L1_B03_.mat';
datapointers.permooring.Spotter.B05.filename = 'roxsi_spotter_L1_B05_.mat';
datapointers.permooring.Spotter.X01.filename = 'roxsi_spotter_L1_X01_.mat';
datapointers.permooring.Spotter.X03.filename = 'roxsi_spotter_L1_X03_.mat';
datapointers.permooring.Spotter.X04.filename = 'roxsi_spotter_L1_X04_.mat';

%
datapointers.permooring.smartspotter.E08.filename1 = fullfile(datapointers.instrumentdir.smartspotter, 'roxsi_spotter_L1_X04_.mat');
datapointers.permooring.smartspotter.E08.filename2 = fullfile(datapointers.instrumentdir.spotter, '.mat');


%%

% % datapointers.permooring.SoloD

% Bottom depth
% varfieldnames.soloD.bottomdepth = 'mean_depth';
% varfieldnames.soloD.timespec = 'dtime';
% varfieldnames.soloD.frequency = 'freq';
% spectra
%
varfieldnames.smartspotter.bottomdepth = {'pressuredata', 'bottomdepth'};

% Time vector

% Frequency vector

% Spectra
varfieldnames.smartspotter.See = {'spectra', 'See'};


%%


% % % If no recomputation is necessary
% % varfieldnames.soloD.bottomdepth = 'mean_depth';
% % varfieldnames.soloD.timespec = 'dtime';
% % varfieldnames.soloD.frequency = 'freq';
% % varfieldnames.soloD.HsigSS = 'HsigSS';
% % varfieldnames.soloD.meanfreqSS = 'meanfreqSS';
% % 
% % 
% % % % 
% % % % v3_dynfld = {'averaged', char(uvw_names(indrow_aux, 3))};
% % % % v3plt_var = getfield(aquadoppL1, v3_dynfld{:});

%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------

%%

%
list_typesinstr = fieldnames(datapointers.permooring);

%
ind_matchmooring = 0;

%
lloop_aux = true;
%
while lloop_aux && (ind_matchmooring < length(list_typesinstr))
    %
    ind_matchmooring = ind_matchmooring + 1;
    %
    if any(strcmp(datapointers.permooring.(list_typesinstr{ind_matchmooring}).mooringID, mooringID))
        lloop_aux = false;
    end
end

%
if lloop_aux
    error(['Mooring ' mooringID ' not found!'])
end


%%



%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------

%%

% % %
% % dirdata_load = datapointers.instrumentdir.(list_typesinstr{ind_matchmooring});


%%

% ------------------------
%
if ~strcmp(list_typesinstr{ind_matchmooring}, 'smartspotter')

    %
    varMooring = who('-file', datapointers.permooring.(list_typesinstr{ind_matchmooring}).filename);

    %
    data_aux = load(datapointers.permooring.(list_typesinstr{ind_matchmooring}).filename);

    %
    datamooring = data_aux.(varMooring{1});

    %
    clear data_aux


% ------------------------
else

% %     %
% %     data_1 = load(datapointers.permooring.smartspotter.filename1);
% %     data_2 = load(datapointers.permooring.smartspotter.filename2);

    %
    varBuoy = who('-file', datapointers.permooring.smartspotter.filename2);
    varPres = who('-file', datapointers.permooring.smartspotter.filename1);
    
    %
    data_buoy_aux = load(datapointers.permooring.smartspotter.filename2);
    data_pres_aux = load(datapointers.permooring.smartspotter.filename1);
    
    %
    data_buoy_aux = data_buoy_aux.(varBuoy{1});
    data_pres_aux = data_pres_aux.(varPres{1});

    %
    datamooring = data_buoy_aux;    
    % Remove stuff???
    
    % Copy pressure data

    %
    list_fields = {'latitude', 'longitude', 'zhab', 'dt', ...
                   'dtime', 'pressure'};
    
    %
    for i = 1:length(list_fields)
        datamooring.pressuredata.(list_fields{i}) = data_pres_aux.(list_fields{i});
    end

    %
    clear data_buoy_aux data_pres_aux
end


%%
% ------------------------------------------------
% -------- MAY COMPUTE CERTAIN QUANTITIES --------
% ---------- DEPENDING ON THE INSTRUMENT ---------
% ------------------------------------------------

%%

g = 9.8;
rho0 = 1030;

%% For SoloDs

%% For ???

%% For smart moorings


%
if strcmp(list_typesinstr{ind_matchmooring}, 'smartspotter')

    % ------------------------------------------------

    %
    datamooring.pressuredata.bottomdepth = 1e4 * datamooring.pressuredata.pressure / (g*rho0);


    % ------------------------------------------------
    % 

    % Get appropriate parameters for averaging
% %     (windowavg, timelims, dt_step)
    windowavg = 3600;
    time_lims = datamooring.spectra.dtime([1, end]);
    

    %
    [datamooring.spectra.bottomdepth, ~, ~] = ...
                    time_smooth_reg(datamooring.pressuredata.dtime, ...
                                    datamooring.pressuredata.bottomdepth, ...
                                    windowavg, time_lims);

    %

end



%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------

%%

%
frequency = getfield(datamooring, varfieldnames.(list_typesinstr{ind_matchmooring}).frequency{:});
elevationSpectra = getfield(datamooring, varfieldnames.(list_typesinstr{ind_matchmooring}).See{:});


%%

%
linfreqlims = (frequency >= freq_lims(1)) & (frequency <= freq_lims(2));


%% Make dimensions consistent

%% Recompute bulk statistics (significant wave height and mean period)

%
frequency_matrix_aux = repmat(frequency(linfreqlims), 1, size(elevationSpectra, 2));

%
m0 = trapz(frequency(linfreqlims), elevationSpectra(linfreqlims, :), 1);
m1 = trapz(frequency(linfreqlims), elevationSpectra(linfreqlims, :) .* frequency_matrix_aux, 1);


%
Hsig = 4.*sqrt(m0);
T_mean = m0./m1;


%%
% ------------------------------------------------
% ------------ NOW COMPUTE GROUP VELOCITY ------------------------------------
% ------------------------------------------------
