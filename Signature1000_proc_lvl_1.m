%% Script that processes (from RAW to level 1) the Signature1000
% data from the large-scale array in the ROXSI 2022 experiment.
% 
% This script does all the processing, while
% Signature1000_proc_scalars_lvl_1.m and Signature1000_proc_velocity_lvl_1.m
% are adapted parts of the current script. This script creates
% different files. 
%
% All Signatures (ADCPs) were programmed to start
% sampling at 2022/06/21 18:00:00 (local). Signature
% SN 103045, at B10, was only deployed later due to
% weather conditions.
%
% This script DOES NOT process high-resolution or echosounder
% data from the Signatures.

clear
close all


%%
% --------------------------------------
% -------- SET DIRECTORY PATHS ---------
% --------------------------------------


%% Raw data directory

dirparent_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
dir_data_raw = fullfile(dirparent_data, 'RAW', 'Signature1000');


%% Output directory

dir_output_L1 = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Signature_Level1/';


%%
% -------------------------------------------
% --- DEFINE VARIABLES FOR DATA PROCESSING --
% -------------------------------------------

%% Load ADCP deployment information

% dir_coderepo = repo_dirpath();
dir_coderepo = '/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc';
load(fullfile(dir_coderepo, 'deploymentInfo_ROXSI2022.mat'), 'deploymentInfo_ROXSI2022')


%% List of Signatures that will be processed

% % % All Signatures
% % list_Signature = {'A01_103043', ...
% %                   'B10_103045', ...
% %                   'B13_103046', ...
% %                   'B15_103056', ...
% %                   'B17_101923', ...
% %                   'C01_102128', ...
% %                   'X11_101941', ...
% %                   'X05_100231'};
% % 
% % % Just a test
% % list_Signature = {'A01_103043'};
list_Signature = {'B10_103045'};

%
Nsignatures = length(list_Signature);


%% If you want to test the code/process a small
% amount of data, set this time limits accordingly.
% Otherwise, comment this out and all the data
% will be processed between deployment and 
% recovery

%
time_lims_proc = [datetime(2022, 07, 01, 00, 00, 00), ...
                  datetime(2022, 07, 03, 00, 00, 00)];
time_lims_proc.TimeZone = 'America/Los_Angeles';

%
Ndatasegments = size(time_lims_proc, 1);

%
if exist('time_lims_proc', 'var')
    lpredeflimits = true;
else
    lpredeflimits = false;
end
    

%% The list of ADCPs processed with either 4 or 5 beams

%
list_4beams = {'B10'};    % B10 has 5th beam, but on HR mode (up to 4 meters above the transducer). Others did the echosounder.
list_5beams = {'A01', 'B13', 'B15', 'B17', 'C01', 'X05', 'X11'};

% Check all of the Signatures in the list are in either
% list for Signatures with 4 or 5 beams
for i = 1:length(list_Signature)

    %
    linanylist = any(contains(list_4beams, list_Signature{i}(1:3))) || ...
                 any(contains(list_5beams, list_Signature{i}(1:3)));

    %
    if ~linanylist

        error(['Signature ' list_Signature{i} ' is not present in ' ...
               'any of the lists specifying velocity processing.'])

    end

end


%% List of all variables/field names in the data structure
% as data is first stored in data structure

%
list_rawdata = {'timedatenum', 'pressure', 'temperature', ...
                'heading', 'pitch', 'roll', ...
                'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
                'amp1', 'amp2', 'amp3', 'amp4', 'amp5', ...
                'cor1', 'cor2', 'cor3', 'cor4', 'cor5', ...
                'timedatenum5'};

%%
% -------------------------------------
% --------- LOAD ATM PRESSURE ---------
% -------------------------------------

%% Load atmospheric pressure

%
atm_pressure = load(fullfile(dirparent_data, 'RAW', ...
                      'noaa_mry_barometric_pressure', 'atm_pressure.mat'));
% Make correction for the spatial variability
atm_pressure.atm_pres = atm_pressure.atm_pres - (100*0.032);    % the correction is 0.032 dbar, and atmospheric pressure is in milibar

% Time limits when pressure offset was set for Aquadopps
time_lims = [datetime(2022, 06, 13, 08, 0, 0), ...
             datetime(2022, 06, 14, 15, 0, 0)];
time_lims.TimeZone = atm_pressure.time_vec.TimeZone;
%
lintimelims = (atm_pressure.time_vec >= time_lims(1)) & ...
              (atm_pressure.time_vec <= time_lims(2));

%
atmpresanomaly.dtime = atm_pressure.time_vec;
atmpresanomaly.atm_anomaly = (atm_pressure.atm_pres - mean(atm_pressure.atm_pres(lintimelims)))./100;
atmpresanomaly.units = 'dbar';

% Remove NaNs in atmospheric pressure (just a few, very sparse,
% such that interpolation over them is fine)
lgood = ~isnan(atmpresanomaly.atm_anomaly);
%
atmpresanomaly.dtime = atmpresanomaly.dtime(lgood);
atmpresanomaly.atm_anomaly = atmpresanomaly.atm_anomaly(lgood);


% % % % Check the atmospheric pressure correction
% % %
% % figure
% %     plot(atm_pressure.time_vec, atm_pressure.atm_pres)
% %     hold on
% %     overlayline('v', time_lims, '--k')
% % 
% % %
% % figure
% %     plot(atmpresanomaly.dtime, atmpresanomaly.atm_anomaly)
% %     hold on
% %     overlayline('v', time_lims, '--k')


%%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% -------------------------- DO DATA PROCESSING --------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_Signature_procL1_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_output_L1, log_file_name))

%
totalRunTime = tic;


%% Display on the screen:

%
disp(' '), disp(' ')
disp('------------------------------ Processing scalar variables from Signature1000s ------------------------------')
disp('List of Signature1000s being processed:')
%
for i = 1:Nsignatures
    disp([num2str(i) ' - ' list_Signature{i}])
end


%% Loop over Signatures in the list and process dat

% Loop over Signature1000's in the list
for i1 = 1:Nsignatures


    %%
    % ------------------------------------------
    % ----- PREAMBLE OF THE DATA PROCESSING ----
    % ------------------------------------------

    %% Display progress messahe on the screen
    %
    disp(' '), disp(' ')
    disp(['--- Start processing Signature1000 data from: ' list_Signature{i1} ' ---'])


    %% Get the Signature1000 files in the correct order

    %
    dir_data_aux = fullfile(dir_data_raw, list_Signature{i1}, 'converted');
    %
    list_dir_aux = Signature_orderfiles(dir_data_aux, list_Signature{i1}(1:3));

    %
    Nfiles = length(list_dir_aux);


    %% First add basic metadata to structure variable

    %
    sigL1.SN = convertCharsToStrings(list_Signature{i1}(5:end));
    sigL1.mooringID = convertCharsToStrings(list_Signature{i1}(1:3));

    % Latitude/longitude
    info_mooringtable = ROXSI_mooringlocation(sigL1.mooringID, "ADCP");
    %
    sigL1.latitude = info_mooringtable.latitude;
    sigL1.longitude = info_mooringtable.longitude;

    %
    sigL1.site = info_mooringtable.roxsiarray;

    %
    [sigL1.X, sigL1.Y] = ROXSI_lltoxy(sigL1.latitude, sigL1.longitude, sigL1.site);

    %
    sigL1.coordsystem = 'localXY';

    % In clockwise degrees from the true north
    sigL1.magdec = 12.86;


    %%

    %
    if any(contains(list_5beams, list_Signature{i1}(1:3)))
        sigL1.l5beams = true;
    else
        sigL1.l5beams = false;
    end

    
    %% Get trimming times for the deployment period of this Signature

    %
    ind_row_match = find(strcmp(deploymentInfo_ROXSI2022.SN, list_Signature{i1}(5:end)));
    
    %
    time_1 = deploymentInfo_ROXSI2022.time_begin_trim(ind_row_match);
    time_2 = deploymentInfo_ROXSI2022.time_end_trim(ind_row_match);

    %
    time_1 = datetime(datenum(time_1, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    time_2 = datetime(datenum(time_2, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');


    %% Define trimming edges as time processing edges

    %
    if ~lpredeflimits
        %
        time_lims_proc = [time_1, time_2];
        Ndatasegments = size(time_lims_proc, 1);
    end


    %%
    % ------------------------------------------
    % ------ LOAD CONFIGURATION INFO FROM ------
    % ------------- ONE DATA FILE --------------
    % ------------------------------------------

    %%

    % Just the Config structure variable
    dataread_aux = load(fullfile(dir_data_aux, list_dir_aux(1).name), 'Config');

    %
    sigL1.Config = dataread_aux.Config;


    %% Get height above the bottom (zhab) of the ADCP bins

    % In meters (based on the Solidworks drawing)
    sigL1.transducerHAB = (31.88)/100;

    % In meters
    sigL1.binsize = dataread_aux.Config.Burst_CellSize;
    
    % Height of the first cell center relative to transducer
    % (based on the Principles of Operation manual by Nortek, page 12)
    cellcenter_first_bin = dataread_aux.Config.Burst_BlankingDistance + ...
                           dataread_aux.Config.Burst_CellSize;

    %
    sigL1.cellcenter = cellcenter_first_bin + ...
                        (0:1:(double(dataread_aux.Config.Burst_NCells) - 1)) .* sigL1.binsize;

    %
    sigL1.zhab = sigL1.transducerHAB + sigL1.cellcenter;


% %     % Find all bins that are entirely below the maximum bottom depth
% % 
% % % %     %
% % % %     lin_verticalrange = ((sig1000.zhab + (sig1000.binsize/2)) < ...
% % % %                          max(sig1000.bottomdepthfrompres));
% %     % Or dummy for code developing purposes
% %     lin_verticalrange = true(1, length(sig1000.zhab));
% % 
% %     %
% %     sig1000.zhab = sig1000.zhab(lin_verticalrange);



    %%
    % ------------------------------------------
    % ------------- LOAD THE DATA --------------
    % ------------------------------------------


    %%

    %
    for i2 = 1:Ndatasegments

        %%

        % If the time segment is a time that is fully outside the
        % deployment trimming edges, than skip this time segment
        if (time_1 >= time_lims_proc(i2, 2)) || (time_2 < time_lims_proc(i2, 1))
            continue
        end 
       
        %
        time_lims_aux = datenum(time_lims_proc(i2, :));
        %
        [listfiles_perseg, struct_timelims] = Signature1000_filesintimelims(list_Signature{i1}(5:end), time_lims_aux);
        Nfilesperseg = length(listfiles_perseg);

% %         %
% %         disp(['--- Raw data that will processed is split into ' num2str(Nfilesperseg) ' data files. Loading data from ---'])
% %         listfiles_perseg

        %% Pre-allocate space for variables that will be read

        % Pre-allocate in cell arrays
        prealloc_aux = cell(Nfilesperseg, 1);

        %
        sigL1.timedatenum = prealloc_aux;
        sigL1.samplingrateHz = NaN;    % only added at the end
        %
% %         sigL1.dtime = prealloc_aux;    % variable not filled in the loop
        sigL1.pressure = prealloc_aux;
        sigL1.temperature = prealloc_aux;
        %
        sigL1.heading = prealloc_aux;
        sigL1.pitch = prealloc_aux;
        sigL1.roll = prealloc_aux;
    
        %
        sigL1.vel1 = prealloc_aux;
        sigL1.vel2 = prealloc_aux;
        sigL1.vel3 = prealloc_aux;
        sigL1.vel4 = prealloc_aux;
        %
        sigL1.amp1 = prealloc_aux;
        sigL1.amp2 = prealloc_aux;
        sigL1.amp3 = prealloc_aux;
        sigL1.amp4 = prealloc_aux;
        %
        sigL1.cor1 = prealloc_aux;
        sigL1.cor2 = prealloc_aux;
        sigL1.cor3 = prealloc_aux;
        sigL1.cor4 = prealloc_aux;

        %
        if sigL1.l5beams
            %
            sigL1.timedatenum5 = prealloc_aux;
            sigL1.dtime5 = prealloc_aux;
            %
            sigL1.vel5 = prealloc_aux;
            sigL1.amp5 = prealloc_aux;
            sigL1.cor5 = prealloc_aux;
        end
        

        %% Finally load the raw data

        %
        tic
        % Load over all files with data in the i2'th segment
        for i3 = 1:Nfilesperseg

            %
            dataread_aux = load(fullfile(dir_data_aux, ...
                                         listfiles_perseg(i3)));

            %
            lin_proclims_beam4time_aux = (dataread_aux.Data.Burst_Time >= datenum(time_lims_proc(i2, 1))) & ...
                                         (dataread_aux.Data.Burst_Time <  datenum(time_lims_proc(i2, 2)));

            % -------------------------------
            % First get all the scalars
            %
            sigL1.timedatenum{i3} = dataread_aux.Data.Burst_Time(lin_proclims_beam4time_aux);
            sigL1.pressure{i3} = dataread_aux.Data.Burst_Pressure(lin_proclims_beam4time_aux);
            sigL1.temperature{i3} = dataread_aux.Data.Burst_Temperature(lin_proclims_beam4time_aux);
            %
            sigL1.heading{i3} = dataread_aux.Data.Burst_Heading(lin_proclims_beam4time_aux);
            sigL1.pitch{i3} = dataread_aux.Data.Burst_Pitch(lin_proclims_beam4time_aux);
            sigL1.roll{i3} = dataread_aux.Data.Burst_Roll(lin_proclims_beam4time_aux);
            

            % -------------------------------
            % Then get velocity data
            
            % Dummy/for code development
            lin_verticalrange = true(1, size(dataread_aux.Data.Burst_VelBeam1, 2));

            % Get data from the 4 beams
            sigL1.vel1{i3} = dataread_aux.Data.Burst_VelBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel2{i3} = dataread_aux.Data.Burst_VelBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel3{i3} = dataread_aux.Data.Burst_VelBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel4{i3} = dataread_aux.Data.Burst_VelBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            %
            sigL1.amp1{i3} = dataread_aux.Data.Burst_AmpBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.amp2{i3} = dataread_aux.Data.Burst_AmpBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.amp3{i3} = dataread_aux.Data.Burst_AmpBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.amp4{i3} = dataread_aux.Data.Burst_AmpBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            %
            sigL1.cor1{i3} = dataread_aux.Data.Burst_CorBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.cor2{i3} = dataread_aux.Data.Burst_CorBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.cor3{i3} = dataread_aux.Data.Burst_CorBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.cor4{i3} = dataread_aux.Data.Burst_CorBeam4(lin_proclims_beam4time_aux, lin_verticalrange);

            %
            if sigL1.l5beams
                %
                lin_proclims_beam5time_aux = (dataread_aux.Data.IBurst_Time >= datenum(time_lims_proc(i2, 1))) & ...
                                             (dataread_aux.Data.IBurst_Time <  datenum(time_lims_proc(i2, 2)));
                %
                sigL1.vel5{i3} = dataread_aux.Data.IBurst_VelBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
                sigL1.amp5{i3} = dataread_aux.Data.IBurst_AmpBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
                sigL1.cor5{i3} = dataread_aux.Data.IBurst_CorBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
    
                %
                sigL1.timedatenum5{i3} = dataread_aux.Data.IBurst_Time(lin_proclims_beam5time_aux);
            end
            
            %
            disp(['--- Done loading data from file ' num2str(i3) ' out of ' num2str(Nfilesperseg) ' ---'])
        end


        %% Concatenate data in cell arrays into matrices

        %
        for i3 = 1:length(list_rawdata)
            %
            if isfield(sigL1, list_rawdata{i3})
                sigL1.(list_rawdata{i3}) = cat(1, sigL1.(list_rawdata{i3}){:});
            end
        end

% % % % %   TO BE DELETED!!!
% %         % ------------------------------------------------
% %         % Concatenate cell array into a long column vector
% %         
% %         sigL1.pressure = cat(1, sigL1.pressure{:});
% %         sigL1.temperature = cat(1, sigL1.temperature{:});
% %         %
% %         sigL1.heading = cat(1, sigL1.heading{:});
% %         sigL1.pitch = cat(1, sigL1.pitch{:});
% %         sigL1.roll = cat(1, sigL1.roll{:});
% %         % ------------------------------------------------
% % 
% %         % Concatenate cell arrays into matrices (loop over variables)
% %         for i3 = 1:length(list_beam_vars)
% %             sigL1.(list_beam_vars{i3}) = cat(1, sigL1.(list_beam_vars{i3}){:});
% %         end
% %         sigL1.timednum_fourbeams = cat(1, sigL1.timednum_fourbeams{:});
% % 
% %         % Concatenate 5th beam data
% %         if any(contains(list_5beams, list_Signature{i1}(1:3)))
% %             %
% %             sigL1.timednum_beam5 = cat(1, sigL1.timednum_beam5{:});
% %             %
% %             sigL1.vel5_raw = cat(1, sigL1.vel5_raw{:});
% %             sigL1.amp5_raw = cat(1, sigL1.amp5_raw{:});
% %             sigL1.corr5_raw = cat(1, sigL1.corr5_raw{:});
% %         end

        %
        disp('--- Done with loading all of the data for the current time segment in: ---')
        toc

    end

    %%
    % ------------------------------------------
    % ---------- BASIC DATA PROCESSING ---------
    % ------------------------------------------

    %% Clock drift correction

    %
    disp('----- Correcting for clock drift -----')

    % Now correct for clock drift
    time_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, ...
                                            list_Signature{i1}(5:end), ...
                                            sigL1.timedatenum);

    % Now correct for clock drift
    if sigL1.l5beam
        time5_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, ...
                                                 list_Signature{i1}(5:end), ...
                                                 sigL1.timedatenum5);
        keyboard
    end


    %% Convert time to date time

    %
    sigL1.dtime = datetime(time_aux, 'ConvertFrom', 'datenum', ...
                                     'TimeZone', 'America/Los_Angeles');

    %
    if sigL1.l5beam
        %
        sigL1.dtime5 = datetime(time5_aux, 'ConvertFrom', 'datenum', ...
                                           'TimeZone', 'America/Los_Angeles');
        
        keyboard
    end


    %% Apply (small) correction due to atmospheric pressure variability.
    % 
    % As opposed to Aquadopps, Signatures don't have a customizable
    % pressure offset (so pressure measurements don't include the
    % full pressure at the sensor)

    % Interpolate atmospheric pressure anomaly to timestamps
    % of the signature
    atmpresanomaly_aux = interp1(atmpresanomaly.dtime, ...
                                 atmpresanomaly.atm_anomaly, sigL1.dtime);

    %
    sigL1.pressure = sigL1.pressure - atmpresanomaly_aux;


    %% Check clock/gaps

    %% Interpolate 5th beam to the same timestamps as the other 4

% HOW IS THE W CALCULATION FROM 5TH BEAM?????

% %     %
% %     if any(contains(list_5beams, list_Signature{i1}(1:3)))
% %         tic
% %         %
% %         disp(['--- 5 beams will be used for velocity transformation. ' ...
% %               'First will interpolate 5th beam to time stamps of ' ...
% %               'other 4 beams ---'])
% % 
% %         %
% %         sig1000.vel5 = NaN(size(sig1000.vel1));
% %         sig1000.amp5 = sig1000.vel5;
% % 
% %         % Loop over bins of the 5th beam
% %         for i2 = 1:size(sig1000.vel5, 2)
% %             
% %             %
% %             sig1000.vel5(:, i2) = interp1(sig1000.timednum_beam5, ...
% %                                           sig1000.vel5_raw(:, i2), ...
% %                                           sig1000.timednum_fourbeams);
% %             %
% %             sig1000.amp5(:, i2) = interp1(sig1000.timednum_beam5, ...
% %                                           sig1000.amp5_raw(:, i2), ...
% %                                           sig1000.timednum_fourbeams);
% %         end
% %         disp('Done with 5th beam interpolation.') 
% %         toc
% %     end


    %% Interpolate variables to gridded time vector (after
    % making sure there are no major issues above)

    tic
    %
    disp('--- Gridding variables to time grid ---')

    % As time edges of the grid, round to the whole minute (done in a more
    % complicated way so that it doesn't rely on Matlab versions newer than
    % at least 2021b)
    dtime_edge_1 = datetime(sigL1.dtime(1).Year, sigL1.dtime(1).Month, sigL1.dtime(1).Day, ...
                            sigL1.dtime(1).Hour, (sigL1.dtime(1).Minute + 1), 00);
    dtime_edge_2 = datetime(sigL1.dtime(end).Year, sigL1.dtime(end).Month, sigL1.dtime(end).Day, ...
                            sigL1.dtime(end).Hour, sigL1.dtime(end).Minute, 00);
    %
    df_sampling = double(sigL1.Config.Burst_SamplingRate);    % in Hertz

    %
    dtime_grid = dtime_edge_1 : seconds(1/df_sampling) : dtime_edge_2;
    dtime_grid.TimeZone = sigL1.dtime.TimeZone;

    %
    Nlengthtimeseries = length(sigL1.dtime);
    list_time_vars = {'dtime', 'dtime5', 'timedatenum', 'timedatenum5'};

    %
    for i2 = 1:length(list_variables_aux)

        % Only interpolate fields that are NOT contained in list_time_vars
        % and are vectors of the correct length
        if ~any(contains(list_time_vars, list_variables_aux{i2})) && ...
           (length(sigL1.(list_variables_aux{i2})) == Nlengthtimeseries)
           
            %
            sigL1.(list_variables_aux{i2}) = ...
                            interp1(sigL1.dtime, ...              
                                    sigL1.(list_variables_aux{i2}), ...
                                    dtime_grid);

            % Turn to column vector
            sigL1.(list_variables_aux{i2}) = sigL1.(list_variables_aux{i2})(:);

        end
    end

    % Replace measured time stamps by time grid
    sigL1.dtime = dtime_grid(:);

    %
    disp('--- Done with time gridding ---')
    toc

    % Remove datenum time
    sigL1 = rmfield(sigL1, 'timedatenum');

    % Add sampling rate as a field
    sigL1.samplingrateHz = df_sampling;


    %% Compute bottom depth from pressure


    %%
    % ------------------------------------------
    % ---------- PROCESS VELOCITY DATA ---------
    % ------------------------------------------

    %% Transpose matrices so that row dimension is along bins

    %% Compute magnetic-ENU 3 components of velocity
    % (using library ADCPtools)
    %
    % Compute ENU velocities in different ways depending
    % whether 4 or 5 beams should be used (check the preamble
    % of this script. The preamble also makes sure that
    % all Signatures being processed are in the processing lists)
    %
    % Very long timeseries (e.g. a couple of weeks of 8 Hz
    % Signature1000 data (need to be broken into smaller chunks
    % otherwise janus5beam2earth will use all memory and crash Matlab).


    %% Rotate horizontal velocity from magnetic to true north

    %% QC velocity based on backscatter and/or correlation (????)

    %%
    % ------------------------------------------
    % ---- COMPUTE LOW-FREQUENCY QUANTITIES ----
    % ------------------------------------------

    %%


    %%
    % -----------------------------------------------------------
    % -------- FINAL ADJUSTMENTS, SAVE DATA AND QC PLOTS --------
    % -----------------------------------------------------------
    % ----------------------------------------------------
    % Turn all vectors into column vectors so that Matlab
    % can quickly displace the structure variable in the 
    % command window (Matlab displays first elements of
    % row vectors, and it takes longer)
%     list_fields_aux = fieldnames(aquadoppL1);
%     %
%     for i2 = 1:length(list_fields_aux)
%         %
%         if isvector(aquadoppL1.(list_fields_aux{i2})) && ...
%            ~isstruct(aquadoppL1.(list_fields_aux{i2})) && ...
%            ~ischar(aquadoppL1.(list_fields_aux{i2}))
%             %
%             aquadoppL1.(list_fields_aux{i2}) = aquadoppL1.(list_fields_aux{i2})(:);
%         end
%     end

    %% Add README

    % %     %
% %     time_dataproc = datetime('now', 'TimeZone', 'Local');
% %     time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');
% %     % Add a general README
% %     sig1000.README = ['Level 1 Signature1000 data from ROXSI 2022. The data is from Signature ' ...
% %                          ' with serial number SN and deployed at mooring site mooringID. ' ...
% %                          'Data processed by script ' mfilename() '.m on ' time_dataproc_char ' (TimeZone ' time_dataproc.TimeZone '). ' ...
% %                          'Horizontal velocity components are relative to ??????, where the magnetic ' ...
% %                          'decliation was taken from www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml. Data ' ...
% %                          'in the Level 1 structure has been trimmed for the deployment ' ...
% %                          'period, as defined in the table deploymentInfo_ROXSI2022.mat. ' ...
% %                          'Pressure is in dbar, where atmospheric pressure has been removed.'];


    %% Save stuff


    % ----------------------------------------------------
    % Save level 1 data with all variables

    % ----------------------------------------------------
    % Save level 1 data with scalar variables

    % ----------------------------------------------------
    % Save level 1 data along-beam/backscatter/correlation quantities

    % ----------------------------------------------------
    % Save QC figures


    %% Stuff before start processing the next Signature in the list

        % ----------------------------------------------------
    %
    disp(['----- DONE with RAW to L1 Signature1000 data processing: ' list_Signature{i1} ' -----'])
    toc(totalRunTime)
    
    %
    close all
    clear sigL1

end


%%

%
disp('###################### Done with processing from RAW to L1 for all Signature1000 ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');

