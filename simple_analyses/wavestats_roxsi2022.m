function mooringwavestats = wavestats_roxsi2022(mooringID, freq_lims, lsmartmooringusepres)
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

%%

addpath(genpath('/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc/'))


%%

if ~exist('lsmartmooringusepres', 'var')
    lsmartmooringusepres = false;
end


%% Load mooring table

% % %
% % mooringtable = load(['/home/omarques/Documents/MATLAB' ...
% %                      '/roxsi_largescale_dataproc/ROXSI2022_mooringtable.mat']);
% % mooringtable = mooringtable.mooringtable;


%% Define directories

%
% % dirparent_data_bw = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
dirparent_data_local = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/';

%
datapointers.instrumentdir.SoloD = fullfile(dirparent_data_local, 'Level2_Data', 'SoloD_Level2');

%
datapointers.instrumentdir.spotter = fullfile(dirparent_data_local, 'Level1_Data', 'Spotter_Level1');
datapointers.instrumentdir.smartspotter = fullfile(dirparent_data_local, 'Level1_Data', 'Spotter_Smart_Level1', 'gridded/');
%
datapointers.instrumentdir.aquadopp = fullfile(dirparent_data_local, 'Level2_Data', 'Aquadopp_Level2');
datapointers.instrumentdir.signature = fullfile(dirparent_data_local, 'Level2_Data', 'Signature_Level2');

% % %
% % dirbuoy = fullfile(dirparent_data, 'Spotter_Level1');
% % dirpressure = fullfile(dirparent_data, 'Spotter_Smart_Level1', 'gridded');


%%

%
list_SoloD = {'X07', 'X08', 'X09', 'X10', 'X12', 'X14', ...
              'A02', 'A04', 'A05', 'A06', 'A07', 'A09', ...   % no data for 'A08'
              'B04', 'B06', 'B07', 'B09', 'B12', 'B14', 'B16', 'B18', ...
              'C02', 'C03', 'C05', 'C06', 'C07', 'C08', 'C09', ...   % no data for C04
              'D01', 'D02'};

%
list_Spotter = {'B01', 'B03', 'B05', 'X01', 'X03', 'X04'};
%
list_SmartMooring = {'E01', 'E02', 'E05', 'E07', 'E08', 'E09', 'E10', 'E11', 'E13'};

%
list_Aquadopp = {'B08', 'B11', ...    % A03 should have had a SoloD, but I missed it
                 'E03', 'E04', 'E06', 'E12', ...
                 'X06', 'X13', ...
                 'F01', 'F02', 'F03', 'F04', 'F05'};

%
list_Signature = {'A01', ...
                  'B10', 'B13', 'B15', 'B17', ...
                  'C01', ...
                  'X05', 'X11'};

%
datapointers.permooring.SoloD.mooringID = list_SoloD;
datapointers.permooring.SoloD.datadir = datapointers.instrumentdir.SoloD;

%
datapointers.permooring.Spotter.mooringID = list_Spotter;
datapointers.permooring.Spotter.datadir = datapointers.instrumentdir.spotter;
datapointers.permooring.SmartMooring.mooringID = list_SmartMooring;
datapointers.permooring.SmartMooring.datadir = datapointers.instrumentdir.smartspotter;

%
datapointers.permooring.Aquadopp.mooringID = list_Aquadopp;
datapointers.permooring.Aquadopp.datadir = datapointers.instrumentdir.aquadopp;

%
datapointers.permooring.Signature.mooringID = list_Signature;
datapointers.permooring.Signature.datadir = datapointers.instrumentdir.signature;


%%

%
for i = 1:length(list_SoloD)
    %
    datapointers.permooring.SoloD.(list_SoloD{i}).filename = fullfile(datapointers.instrumentdir.SoloD, ['roxsi_soloD_L2_' list_SoloD{i} '.mat']);
end

%
datapointers.permooring.Spotter.B01.filename = 'roxsi_spotter_L1_B01_1158.mat';
datapointers.permooring.Spotter.B03.filename = 'roxsi_spotter_L1_B03_1152.mat';
datapointers.permooring.Spotter.B05.filename = 'roxsi_spotter_L1_B05_1153.mat';
datapointers.permooring.Spotter.X01.filename = 'roxsi_spotter_L1_X01_1151.mat';
datapointers.permooring.Spotter.X03.filename = 'roxsi_spotter_L1_X03_1157.mat';
datapointers.permooring.Spotter.X04.filename = 'roxsi_spotter_L1_X04_1155.mat';

% Need to add other Smart Moorings later
datapointers.permooring.SmartMooring.E02.filename1 = fullfile(datapointers.instrumentdir.smartspotter, 'roxsi_smartmooring_L1_E02sp_1859_gridded.mat');
datapointers.permooring.SmartMooring.E02.filename2 = fullfile(datapointers.instrumentdir.spotter, 'roxsi_spotter_L1_E02_1859.mat');
datapointers.permooring.SmartMooring.E02.filename3 = fullfile(dirparent_data_local, 'Level2_Data', 'Spotter_Smart_Level2', 'roxsi_smartmooring_L2_E02sp_1859.mat');
%
datapointers.permooring.SmartMooring.E08.filename1 = fullfile(datapointers.instrumentdir.smartspotter, 'roxsi_smartmooring_L1_E08sp_1852_gridded.mat');
datapointers.permooring.SmartMooring.E08.filename2 = fullfile(datapointers.instrumentdir.spotter, 'roxsi_spotter_L1_E08_1852.mat');
datapointers.permooring.SmartMooring.E08.filename3 = fullfile(dirparent_data_local, 'Level2_Data', 'Spotter_Smart_Level2', 'roxsi_smartmooring_L2_E08sp_1852.mat');


% % roxsi_spotter_L1_B03_.mat  roxsi_spotter_L1_E07_1857.mat  roxsi_spotter_L1_E13_1849.mat
% %  roxsi_spotter_L1_B05_.mat  roxsi_spotter_L1_E08_.mat  roxsi_spotter_L1_X01_.mat
% %    roxsi_spotter_L1_E01_1851.mat  roxsi_spotter_L1_E09_1850.mat  roxsi_spotter_L1_X03_.mat
% %    roxsi_spotter_L1_E02_1859.mat  roxsi_spotter_L1_E09_1856.mat  roxsi_spotter_L1_X04_.mat
% % roxsi_spotter_L1_B01_1150.mat			    roxsi_spotter_L1_E05_1853.mat  roxsi_spotter_L1_E10_1848.mat
% % roxsi_spotter_L1_B01_1158.mat			    roxsi_spotter_L1_E07_1855.mat  roxsi_spotter_L1_E11_1860.mat

% Aquadopp
datapointers.permooring.Aquadopp.B08.filename = 'roxsi_aquadopp_L2_B08at_13288.mat';
datapointers.permooring.Aquadopp.B11.filename = 'roxsi_aquadopp_L2_B11a_12280.mat';
datapointers.permooring.Aquadopp.E03.filename = 'roxsi_aquadopp_L2_E03a_13300.mat';
datapointers.permooring.Aquadopp.E04.filename = 'roxsi_aquadopp_L2_E04a_13172.mat';
datapointers.permooring.Aquadopp.E06.filename = 'roxsi_aquadopp_L2_E06a_9736.mat';
datapointers.permooring.Aquadopp.E12.filename = 'roxsi_aquadopp_L2_E12a_11150.mat';
%
datapointers.permooring.Aquadopp.X06.filename = 'roxsi_aquadopp_L2_X06a_13290.mat';
datapointers.permooring.Aquadopp.X13.filename = 'roxsi_aquadopp_L2_X13a_9945.mat';
%
datapointers.permooring.Aquadopp.F01.filename = 'roxsi_aquadopp_L2_F01a_9995.mat';
datapointers.permooring.Aquadopp.F02.filename = 'roxsi_aquadopp_L2_F02a_5838.mat';
datapointers.permooring.Aquadopp.F03.filename = 'roxsi_aquadopp_L2_F03a_5384.mat';
datapointers.permooring.Aquadopp.F04.filename = 'roxsi_aquadopp_L2_F04a_5401.mat';
datapointers.permooring.Aquadopp.F05.filename = 'roxsi_aquadopp_L2_F05a_14032.mat';

% Signature
datapointers.permooring.Signature.A01.filename = 'roxsi_signature_L2_A01_103043.mat';
datapointers.permooring.Signature.B10.filename = 'roxsi_signature_L2_B10_103045.mat';
datapointers.permooring.Signature.B13.filename = 'roxsi_signature_L2_B13_103046.mat';
datapointers.permooring.Signature.B15.filename = 'roxsi_signature_L2_B15_103056.mat';
datapointers.permooring.Signature.B17.filename = 'roxsi_signature_L2_B17_101923.mat';
datapointers.permooring.Signature.C01.filename = 'roxsi_signature_L2_C01_102128.mat';
datapointers.permooring.Signature.X05.filename = 'roxsi_signature_L2_X05_100231.mat';
datapointers.permooring.Signature.X11.filename = 'roxsi_signature_L2_X11_101941.mat';


% --------------------------------------------
% Add directories to the file name
%
list_loop = {'Spotter', 'Aquadopp', 'Signature'};
%
for i1 = 1:length(list_loop)
    %
    list_mooringID = datapointers.permooring.(list_loop{i1}).mooringID;
    %
    for i2 = 1:length(list_mooringID)
        %
        datapointers.permooring.(list_loop{i1}).(list_mooringID{i2}).filename = ...
            fullfile(datapointers.permooring.(list_loop{i1}).datadir, ...
                     datapointers.permooring.(list_loop{i1}).(list_mooringID{i2}).filename);
    end
end



%%

% Time vector
varfieldnames.SoloD.dtime = {'dtime'};
varfieldnames.Spotter.dtime = {'spectra', 'dtime'};
varfieldnames.SmartMooring.dtime = {'spectra', 'dtime'};
varfieldnames.Aquadopp.dtime = {'dtime'};
varfieldnames.Signature.dtime = {'dtime'};

% Bottom depth
varfieldnames.SoloD.bottomdepth = {'mean_depth'};
varfieldnames.Spotter.bottomdepth = {'spectra', 'bottomdepth'};
varfieldnames.SmartMooring.bottomdepth = {'spectra', 'bottomdepth'};
varfieldnames.Aquadopp.bottomdepth = {'bottomdepthmean'};
varfieldnames.Signature.bottomdepth = {'bottomdepthmean'};

% Frequency vector
varfieldnames.SoloD.frequency = {'freq'};
varfieldnames.Spotter.frequency = {'spectra', 'frequency'};
varfieldnames.SmartMooring.frequency = {'spectra', 'frequency'};
varfieldnames.Aquadopp.frequency = {'frequency'};
varfieldnames.Signature.frequency = {'frequency'};

% Spectra
varfieldnames.SoloD.See = {'See'};
varfieldnames.Spotter.See = {'spectra', 'See'};
varfieldnames.SmartMooring.See = {'spectra', 'See'};
varfieldnames.Aquadopp.See = {'See'};
varfieldnames.Signature.See = {'See'};


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
% --------- FIND WHICH FILE MATCHES WITH ---------
% ------ THE MOORING ID NUMBER IN THE INPUT ------
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
% ------------------------------------------------
% ---------------- LOAD THE DATA -----------------
% ------------------------------------------------

%%

% % %
% % dirdata_load = datapointers.instrumentdir.(list_typesinstr{ind_matchmooring});


%% Load the data


% ------------------------
% If it's not a Smart Mooring
if ~strcmp(list_typesinstr{ind_matchmooring}, 'SmartMooring')

    %
    varMooring = who('-file', datapointers.permooring.(list_typesinstr{ind_matchmooring}).(mooringID).filename);

    %
    data_aux = load(datapointers.permooring.(list_typesinstr{ind_matchmooring}).(mooringID).filename);

    %
    datamooring = data_aux.(varMooring{1});

    %
    clear data_aux


% ------------------------
% If it is a Smart Mooring
else

    % If the elevation will be taken from pressure data
    if lsmartmooringusepres
        
        %
        varPres = who('-file', datapointers.permooring.SmartMooring.(mooringID).filename3);
        data_pres_aux = load(datapointers.permooring.SmartMooring.(mooringID).filename3);
        data_pres_aux = data_pres_aux.(varPres{1});
        
        %
        datamooring = data_pres_aux;
        
        % Put it in the same format as the other option
        datamooring.spectra.dtime = datamooring.dtime;
        datamooring.spectra.bottomdepth = datamooring.mean_depth;
        datamooring.spectra.frequency = datamooring.frequency;
        datamooring.spectra.See = datamooring.See;
        
        % Remove repetition
        datamooring = rmfield(datamooring, {'dtime', 'mean_depth', 'frequency', 'See'});
 
    % Or if elevation will be taken from the Spotter buoy
    else
        %
        varBuoy = who('-file', datapointers.permooring.SmartMooring.(mooringID).filename2);
        varPres = who('-file', datapointers.permooring.SmartMooring.(mooringID).filename1);

        %
        data_buoy_aux = load(datapointers.permooring.SmartMooring.(mooringID).filename2);
        data_pres_aux = load(datapointers.permooring.SmartMooring.(mooringID).filename1);

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
end


%%
% ------------------------------------------------
% -------- MAY COMPUTE CERTAIN QUANTITIES --------
% ---------- DEPENDING ON THE INSTRUMENT ---------
% ------------------------------------------------

%%

g = 9.8;
rho0 = 1030;

%% For SoloDs -- datetime from datenum

%
if strcmp(list_typesinstr{ind_matchmooring}, 'SoloD')
    datamooring.dtime = datetime(datamooring.time_dnum, 'ConvertFrom', 'datenum');
    datamooring.dtime.TimeZone = 'America/Los_Angeles';
end


%% For Spotters -- bottom depth from z_MSL

%
if strcmp(list_typesinstr{ind_matchmooring}, 'Spotter')

    %
    datamooring.location.bottomdepth = -datamooring.location.z_msl;

    % Get appropriate parameters for averaging
    windowavg = 3600;
    time_lims = datamooring.spectra.dtime([1, end]);
    

    %
    [datamooring.spectra.bottomdepth, ~, ~] = ...
                    time_smooth_reg(datamooring.location.dtime, ...
                                    datamooring.location.bottomdepth, ...
                                    windowavg, time_lims);
end


%% For smart moorings

%
if strcmp(list_typesinstr{ind_matchmooring}, 'SmartMooring')

    %
    if ~lsmartmooringusepres
        
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


    end
end


%%

% ##################################
% REMOVE IN THE FUTURE!!!!
if strcmp(list_typesinstr{ind_matchmooring}, 'Spotter') || ...
   strcmp(list_typesinstr{ind_matchmooring}, 'SmartMooring')
    
    %
    if ~lsmartmooringusepres
        datamooring.spectra.frequency = datamooring.spectra.frequency(1:end-1);
    end
end
% ##################################


%% Compute bottom depth for Aquadopps (remove this in the future)

%
if strcmp(list_typesinstr{ind_matchmooring}, 'Aquadopp')
    datamooring.bottomdepthmean = 1e4 * datamooring.pressuremean / (g*rho0);
end


%% May want to recompute elevation spectra for
% pressure sensors using a different depth


%%

%
disp(['---- Done with getting data from mooring ' mooringID ' ----'])


%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------

%%

%
mooringwavestats.mooringID = mooringID;

%
mooringwavestats.instrument = list_typesinstr{ind_matchmooring};


%%

%
mooringwavestats.dtime = getfield(datamooring, varfieldnames.(list_typesinstr{ind_matchmooring}).dtime{:});

%
mooringwavestats.bottomdepth = getfield(datamooring, varfieldnames.(list_typesinstr{ind_matchmooring}).bottomdepth{:});

%
mooringwavestats.frequency = getfield(datamooring, varfieldnames.(list_typesinstr{ind_matchmooring}).frequency{:});
mooringwavestats.df = mooringwavestats.frequency(2) - mooringwavestats.frequency(1);

%
mooringwavestats.See = getfield(datamooring, varfieldnames.(list_typesinstr{ind_matchmooring}).See{:});


%% Make sure vectors are column vectors

%
list_fields_aux = fieldnames(mooringwavestats);

%
for i = 1:length(list_fields_aux)
    if ~ischar(mooringwavestats.(list_fields_aux{i}))
        if isvector(mooringwavestats.(list_fields_aux{i}))
            mooringwavestats.(list_fields_aux{i}) = mooringwavestats.(list_fields_aux{i})(:);
        end
    end
end


%% Check dimensions!!!

if size(mooringwavestats.See, 1)~=length(mooringwavestats.frequency)
    mooringwavestats.See = mooringwavestats.See.';
end


%%

%
mooringwavestats.freq_lims = freq_lims;
mooringwavestats.linlims = (mooringwavestats.frequency >= freq_lims(1)) & ...
                           (mooringwavestats.frequency <= freq_lims(2));


%% Recompute bulk statistics (significant wave height and mean period)

%
frequency_inlims_aux = mooringwavestats.frequency(mooringwavestats.linlims);

%
frequency_matrix_aux = repmat(frequency_inlims_aux, 1, length(mooringwavestats.dtime));

%
m0 = trapz(frequency_inlims_aux, mooringwavestats.See(mooringwavestats.linlims, :), 1);
m1 = trapz(frequency_inlims_aux, mooringwavestats.See(mooringwavestats.linlims, :) .* frequency_matrix_aux, 1);

%
m0 = m0(:);
m1 = m1(:);

%
mooringwavestats.Hsig = 4.*sqrt(m0);

%
mooringwavestats.Tmean = m0./m1;

%
mooringwavestats.Hsig = mooringwavestats.Hsig(:);
mooringwavestats.Tmean = mooringwavestats.Tmean(:);

