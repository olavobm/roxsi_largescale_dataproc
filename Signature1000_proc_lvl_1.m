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
% dir_output_L1 = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/';


%%
% -------------------------------------------
% --- DEFINE VARIABLES FOR DATA PROCESSING --
% -------------------------------------------

%% Add ADCPtools

roxsi_add_libraries()


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
% % list_Signature = {'B10_103045'};
% % list_Signature = {'B13_103046'};

% All Signatures
list_Signature = {'A01_103043', ...
                  'B10_103045', ...
                  'B13_103046', ...
                  'B15_103056', ...
                  'B17_101923', ...
                  'C01_102128', ...
                  'X11_101941'};
% % list_Signature = {'A01_103043'};
list_Signature = {'X05_100231'};

%
Nsignatures = length(list_Signature);


%% If you want to test the code/process a small
% amount of data, set this time limits accordingly.
% Otherwise, comment this out and all the data
% will be processed between deployment and 
% recovery

% % %
% % time_lims_proc = [datetime(2022, 06, 14, 00, 00, 00), ...
% %                   datetime(2022, 07, 06, 00, 00, 00)];
% % time_lims_proc.TimeZone = 'America/Los_Angeles';

%
time_lims_proc = [datetime(2022, 07, 05, 00, 00, 00), ...
                  datetime(2022, 07, 06, 00, 00, 00)];
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

% % %
% % list_rawdata = {'timedatenum', 'pressure', 'temperature', ...
% %                 'heading', 'pitch', 'roll', ...
% %                 'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
% %                 'amp1', 'amp2', 'amp3', 'amp4', 'amp5', ...
% %                 'cor1', 'cor2', 'cor3', 'cor4', 'cor5', ...
% %                 'timedatenum5'};

% Remove variables to free up memory
list_rawdata = {'timedatenum', 'pressure', 'temperature', ...
                'heading', 'pitch', 'roll', ...
                'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
                'timedatenum5'};


%% List of properties for converting along-beam to ENU
% velocity in Signature data

%
theta_beams = 25;

%
lGimbaled = true;
% lGimbaled = false;

BinmapType = 'none';
% BinmapType = 'linear';
% BinmapType = 'nearest';

%
luse3beams = false;


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
disp('------------------------------ Processing Signature1000 data from RAW to L1 ------------------------------')
disp('List of Signature1000s being processed:')
%
for i = 1:Nsignatures
    disp([num2str(i) ' - ' list_Signature{i}])
end


%% Loop over Signatures in the list and process data

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


    %
    if lpredeflimits
        disp('--- Time bounds for processing data are defined on top of the script. Processing data within these bounds ---')
    else
        disp('--- Predefined time bounds for processing data not found. Processing data for full deployment ---')
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
    % --- OPTIONAL PART: GET SCALAR DATA FILE --
    % --------- TO GET TRIMMING PERIOD ---------
    % ------------------------------------------
    %
    % If processing different segments separately
    % that will later be patched, then you should
    % first process the whole dataset of scalars,
    % and get the vertical trimming. Then this
    % will be applied to all segments, and the
    % matrices can be concatenated together.
    %
    % The other option is to comment out this part,
    % and then uncomment the definition of 
    % lin_verticalrange that appears in the next
    % code block. This bypass the requirement of
    % first processing the scalar data, but patching
    % data segments might fail.

    %
    L1scalar = load(fullfile(dir_output_L1, ...
                   ['roxsi_signature_L1_' char(sigL1.mooringID) '_' ...
                                          char(sigL1.SN) '_scalars.mat']));
    L1scalar = L1scalar.sigL1;

    %
    nzbinstrim = length(L1scalar.zhab);

    %
    lin_verticalrange = false(1, length(sigL1.zhab));
    lin_verticalrange(1:nzbinstrim) = true;


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
        sigL1.dtime = [];    % variable not filled in the loop. I just
                             % want to have its position here
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
% %         %
% %         sigL1.amp1 = prealloc_aux;
% %         sigL1.amp2 = prealloc_aux;
% %         sigL1.amp3 = prealloc_aux;
% %         sigL1.amp4 = prealloc_aux;
% %         %
% %         sigL1.cor1 = prealloc_aux;
% %         sigL1.cor2 = prealloc_aux;
% %         sigL1.cor3 = prealloc_aux;
% %         sigL1.cor4 = prealloc_aux;

        %
        if sigL1.l5beams
            %
            sigL1.timedatenum5 = prealloc_aux;
            sigL1.dtime5 = prealloc_aux;
            %
            sigL1.vel5 = prealloc_aux;
% %             sigL1.amp5 = prealloc_aux;
% %             sigL1.cor5 = prealloc_aux;
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
            
            % Get all bins if scalar data is not loaded above (uncomment this)
% %             lin_verticalrange = true(1, size(dataread_aux.Data.Burst_VelBeam1, 2));

            % Get data from the 4 beams
            sigL1.vel1{i3} = dataread_aux.Data.Burst_VelBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel2{i3} = dataread_aux.Data.Burst_VelBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel3{i3} = dataread_aux.Data.Burst_VelBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel4{i3} = dataread_aux.Data.Burst_VelBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            %
% %             sigL1.amp1{i3} = dataread_aux.Data.Burst_AmpBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
% %             sigL1.amp2{i3} = dataread_aux.Data.Burst_AmpBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
% %             sigL1.amp3{i3} = dataread_aux.Data.Burst_AmpBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
% %             sigL1.amp4{i3} = dataread_aux.Data.Burst_AmpBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
% %             %
% %             sigL1.cor1{i3} = dataread_aux.Data.Burst_CorBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
% %             sigL1.cor2{i3} = dataread_aux.Data.Burst_CorBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
% %             sigL1.cor3{i3} = dataread_aux.Data.Burst_CorBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
% %             sigL1.cor4{i3} = dataread_aux.Data.Burst_CorBeam4(lin_proclims_beam4time_aux, lin_verticalrange);

            %
            if sigL1.l5beams
                %
                lin_proclims_beam5time_aux = (dataread_aux.Data.IBurst_Time >= datenum(time_lims_proc(i2, 1))) & ...
                                             (dataread_aux.Data.IBurst_Time <  datenum(time_lims_proc(i2, 2)));
                %
                sigL1.vel5{i3} = dataread_aux.Data.IBurst_VelBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
% %                 sigL1.amp5{i3} = dataread_aux.Data.IBurst_AmpBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
% %                 sigL1.cor5{i3} = dataread_aux.Data.IBurst_CorBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
    
                %
                sigL1.timedatenum5{i3} = dataread_aux.Data.IBurst_Time(lin_proclims_beam5time_aux);
            end
            
            %
            disp(['--- Done loading data from ' num2str(i3) ' out of ' num2str(Nfilesperseg) ' files ---'])

            % Clear to free up some memory
            clear dataread_aux
        end


        %% Concatenate data in cell arrays into matrices

        %
        for i3 = 1:length(list_rawdata)
            %
            if isfield(sigL1, list_rawdata{i3})
                %
                sigL1.(list_rawdata{i3}) = cat(1, sigL1.(list_rawdata{i3}){:});

                % Turn integers into single precision variables
                % so they can be interpolated (correlations are
                % given as integers)
                if isinteger(sigL1.(list_rawdata{i3}))    
                    sigL1.(list_rawdata{i3}) = single(sigL1.(list_rawdata{i3}));
                end
            end
        end

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
    if sigL1.l5beams
        time5_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, ...
                                                 list_Signature{i1}(5:end), ...
                                                 sigL1.timedatenum5);
    end


    %% Convert time to date time

    %
    disp('--- Converting timestamps from datenum to datetime ---')

    %
    sigL1.dtime = datetime(time_aux, 'ConvertFrom', 'datenum', ...
                                     'TimeZone', 'America/Los_Angeles');

    %
    if sigL1.l5beams
        %
        sigL1.dtime5 = datetime(time5_aux, 'ConvertFrom', 'datenum', ...
                                           'TimeZone', 'America/Los_Angeles');
    end




    %% Check clock/gaps

    % Check for something that should never happen at this point
    if any(isnan(sigL1.timedatenum))
        warning(['###### Signature ' list_Signature{i1} ' has ' ...
                 'invalid (NaN) timestamps ######'])
    end

    % Now do a plot for cheching diff(time)
    disp('--- QC plot checking diff time ---')

    %
    inds_time = 1:length(sigL1.dtime);
    inds_difftime = (inds_time(1:end-1) + inds_time(2:end))./2;

    %
    fig_L1_QC_clock = figure;
    set(fig_L1_QC_clock, 'units', 'normalized')
    set(fig_L1_QC_clock, 'Position', [0.2, 0.2, 0.4, 0.6])
        %
        haxs_1 = axes(fig_L1_QC_clock, 'Position', [0.2, 0.575, 0.6 , 0.325]);
        haxs_2 = axes(fig_L1_QC_clock, 'Position', [0.2, 0.150, 0.6, 0.325]);
        %
        haxs_all = [haxs_1, haxs_2];
        hold(haxs_all, 'on')
        %
        plot(haxs_1, inds_time, sigL1.dtime, '-k')
        plot(haxs_2, inds_difftime, seconds(diff(sigL1.dtime)), '-k')

    %
    set(haxs_all, 'FontSize', 12, 'Box', 'on', ...
                  'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_all, 'XLim', [0, (inds_time(end) + 1)])
    %
    ylim(haxs_1, sigL1.dtime([1, end]))

    %
    xlabel(haxs_2, 'Indices', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    ylabel(haxs_1, 'Time', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, 'seconds', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, ['ROXSI 2022: Signature  ' char(sigL1.mooringID) ' - SN ' ...
                   char(sigL1.SN) ': time and diff(time) (in seconds)'], ...
                  'Interpreter', 'Latex', 'FontSize', 12)
    %
    linkaxes(haxs_all, 'x')


    % Plot horizontal lines for trimming edges
    xlims_aux = xlim(haxs_all(1));
    %
    plot(haxs_all(1), xlims_aux, [time_1, time_1], '--r')
    plot(haxs_all(1), xlims_aux, [time_2, time_2], '--r')
    %
    xlim(haxs_all(1), xlims_aux)
    

    %% Apply clock correction for Signatures that had weird issue

    %
    if strcmp(list_Signature{i1}, 'X05_100231')
        %
        [lok_4beams, Nclockinvs, timelims_clockinvs] = correct_Sig1000_clock(sigL1.dtime);

        %
        [lok_5thbeam, Nclockinvs, timelims_clockinvs] = correct_Sig1000_clock(sigL1.dtime5);

        %
        sigL1.dtime = sigL1.dtime(lok_4beams);
        %
        sigL1.pressure = sigL1.pressure(lok_4beams, :);
        sigL1.temperature = sigL1.temperature(lok_4beams, :);
        sigL1.heading = sigL1.heading(lok_4beams, :);
        sigL1.pitch = sigL1.pitch(lok_4beams, :);
        sigL1.roll = sigL1.roll(lok_4beams, :);
        %
        sigL1.vel1 = sigL1.vel1(lok_4beams, :);
        sigL1.vel2 = sigL1.vel3(lok_4beams, :);
        sigL1.vel3 = sigL1.vel4(lok_4beams, :);
        sigL1.vel4 = sigL1.vel4(lok_4beams, :);

        %
        sigL1.dtime5 = sigL1.dtime5(lok_5thbeam);
        sigL1.vel5 = sigL1.vel5(lok_5thbeam, :);

    end


    %% Remove time variables that won't be used anymore

    sigL1 = rmfield(sigL1, 'timedatenum');
    if sigL1.l5beams
        sigL1 = rmfield(sigL1, 'timedatenum5');
    end
    clear time_aux time5_aux
    

    %% Interpolate 5th beam to the same timestamps as the other 4

    %
    if sigL1.l5beams

        tic
        %
        disp(['--- 5 beams will be used for velocity transformation. ' ...
              'Interpolating 5th along-beam velocity (and backscatter ' ...
              'and correlation) to time stamps of other 4 beams ---'])

        %
        vel5_aux = single(NaN(size(sigL1.vel1)));
% %         amp5_aux = vel5_aux;
% %         cor5_aux = vel5_aux;

        % Loop over bins of the 5th beam
        for i2 = 1:size(sigL1.vel5, 2)
            
            %
            vel5_aux(:, i2) = interp1(sigL1.dtime5, ...
                                      sigL1.vel5(:, i2), ...
                                      sigL1.dtime);
% %             %
% %             amp5_aux(:, i2) = interp1(sigL1.dtime5, ...
% %                                       sigL1.amp5(:, i2), ...
% %                                       sigL1.dtime);
% %             %
% %             cor5_aux(:, i2) = interp1(sigL1.dtime5, ...
% %                                       single(sigL1.cor5(:, i2)), ...
% %                                       sigL1.dtime);
        end

        % Replace
        sigL1.vel5 = vel5_aux;
% %         sigL1.amp5 = amp5_aux;
% %         sigL1.cor5 = cor5_aux;
        
        %
        disp('Done with 5th beam interpolation.') 
        toc

        %
        sigL1 = rmfield(sigL1, 'dtime5');
    end


    %% Interpolate variables to gridded time (after
    % making sure there are no major issues above)
    %
    % At this point, time is always along the row dimension

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
    dtime_grid = dtime_grid(:);

    %
    Nlengthtimeseries = length(sigL1.dtime);
    %
    list_time_vars = {'dtime', 'dtime5', 'timedatenum', 'timedatenum5'};
    list_all_fields = fieldnames(sigL1);

    %
    for i2 = 1:length(list_all_fields)

        % Only interpolate fields that are NOT contained in list_time_vars
        % and are arrays with expected number of points along time
        if ~any(contains(list_time_vars, list_all_fields{i2})) && ...
           (size(sigL1.(list_all_fields{i2}), 1) == Nlengthtimeseries)
           
            %
            sigL1.(list_all_fields{i2}) = ...
                            interp1(sigL1.dtime, ...              
                                    sigL1.(list_all_fields{i2}), ...
                                    dtime_grid);
            % PS: matrices are interpolated column-wise

% %             % Turn to column vector
% %             sigL1.(list_variables_aux{i2}) = sigL1.(list_variables_aux{i2})(:);

        end
    end

    % Replace measured time stamps by time grid
    sigL1.dtime = dtime_grid;

    %
    disp('--- Done with time gridding ---')
    toc

    % Add sampling rate as a field
    sigL1.samplingrateHz = df_sampling;


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


    %% Compute bottom depth from pressure

    % Simplest estimate: assume pressure is hydrostatic
    rho0 = 1030;
    g = 9.8;
    sigL1.bottomdepthfrompres = 1e4*sigL1.pressure ./ (rho0*g);


    %%
    % ------------------------------------------
    % ---------- PROCESS VELOCITY DATA ---------
    % ------------------------------------------


    %% First step is to trim out all bins that never give good data

    % Find all bins that are below the maximum bottom depth
    lin_verticalrange = ((sigL1.zhab + (sigL1.binsize)) < max(sigL1.bottomdepthfrompres));
    %
    % PS: it's more intuitive to use binsize/2, however I don't
    % know what is the precise definition of a "full cell" in
    % the Signature1000 data. From one part of the manual, it seemed
    % that measurements 2*binsize might be taken into 1 cell, but the
    % windowing applied by the preprocessing inside the instrument
    % may remove the edges.

    %
    Nptstime = length(sigL1.dtime);
    Nptsbins = length(sigL1.zhab);
    %
    list_all_fields = fieldnames(sigL1);
    %
    for i2 = 1:length(list_all_fields)
        %
        if (size(sigL1.(list_all_fields{i2}), 1)==Nptstime) && ...
           (size(sigL1.(list_all_fields{i2}), 2)==Nptsbins)
            %
            sigL1.(list_all_fields{i2}) = sigL1.(list_all_fields{i2})(:, lin_verticalrange); 

            % Transpose matrices so that row dimension is along bins
            sigL1.(list_all_fields{i2}) = sigL1.(list_all_fields{i2}).';
        end
    end

    % And trim zhab after trimming the matrices above
    sigL1.zhab = sigL1.zhab(lin_verticalrange);


    %% NaN data at/above the instantaneous ocean surface

    %
    tic
    disp('--- Removing data above the instantaneous ocean surface ---')

    %
%     zhab_halfstep = sigL1.zhab + (sigL1.binsize/2);   % same question as above
    zhab_halfstep = sigL1.zhab + (sigL1.binsize);
% %     Nbins = length(sigL1.zhab);   % not used
    %
    ind_abovesurface = 1:1:(size(sigL1.vel1, 1) * size(sigL1.vel1, 2));
    ind_abovesurface = reshape(ind_abovesurface, size(sigL1.vel1, 1), size(sigL1.vel2, 2));
    %
    for i2 = 1:length(sigL1.dtime)
        %
        ind_withinocean_aux = find((zhab_halfstep < sigL1.bottomdepthfrompres(i2)), 1, 'last');
        %
        ind_abovesurface(1:ind_withinocean_aux, i2) = NaN;
    end
    %
    ind_abovesurface = ind_abovesurface(:);
    ind_abovesurface = ind_abovesurface(~isnan(ind_abovesurface));
    %
    sigL1.vel1(ind_abovesurface) = NaN;
    sigL1.vel2(ind_abovesurface) = NaN;
    sigL1.vel3(ind_abovesurface) = NaN;
    sigL1.vel4(ind_abovesurface) = NaN;
    %
    if sigL1.l5beams
        sigL1.vel5(ind_abovesurface) = NaN;
    end

    %
    disp('--- Done with trimming data ---')
    toc


    %% Compute magnetic-ENU 3 components of velocity
    % (using library ADCPtools)
    %
    % Compute ENU velocities in different ways depending
    % whether 4 or 5 beams should be used (check the preamble
    % of this script. The preamble also makes sure that
    % all Signatures being processed are in the processing lists)
    %

    % ------------------------
    % Very long timeseries (e.g. a couple of weeks of 8 Hz
    % Signature1000 data (need to be broken into smaller chunks
    % otherwise janus5beam2earth will use all memory and crash Matlab).

    %
    npts_rot_TH = 1000000;
    %
    if Nptstime<=npts_rot_TH
        indbreak_rot = [1; Nptstime];
    else
        %
        inds_edges_breakrot = 1 : npts_rot_TH : Nptstime;
        if (inds_edges_breakrot(end)~=Nptstime)
            inds_edges_breakrot = [inds_edges_breakrot, Nptstime];
        end
        %
        indbreak_rot = [inds_edges_breakrot(1:end-1); ...
                        (inds_edges_breakrot(2:end) - 1)];
        % The last one shouldn't have 1 subtracted. Add again
        indbreak_rot(2, end) = indbreak_rot(2, end) + 1;
    end

    % ------------------------
    tic
    %
    disp('--- Converting along-beam velocities to magnetic ENU ---')
    %
    disp(['--- Coordinate transformation will be computed ' ...
          'for ' num2str(size(indbreak_rot, 2)) ' separate chunks ' ...
          'of the timeseries (to avoid crashing Matlab) ---'])

    %
    sigL1.u = NaN(size(sigL1.vel1));
    sigL1.v = sigL1.u;
    sigL1.w = sigL1.u;

% % In the future I might want to add the 2 Janus w estimates
% %     sigL1.Wbeam5 = NaN(size(sigL1.vel1));


    % If the Signature is in the 5-beam list
    if sigL1.l5beams

        % Compute velocity using 5 beams
        disp('Using 5 beams')

        %
        for i2 = 1:size(indbreak_rot, 2)

            %
            disp(['--- Transformation in chunk ' num2str(i2) ' out of ' num2str(size(indbreak_rot, 2)) ' ---'])
            %
            ind_sub_aux = indbreak_rot(1, i2) : indbreak_rot(2, i2);

            % Use w from the 5th beam (it does incorporate ux
            % and uy from the Janus measurements)
            [sigL1.u(:, ind_sub_aux), ...
             sigL1.v(:, ind_sub_aux), ...
             ~, ...
             sigL1.w(:, ind_sub_aux)] = ...    % this is the w from 5th beam
                    janus5beam2earth((sigL1.heading(ind_sub_aux).' - 90), ...
                                     sigL1.roll(ind_sub_aux).', -sigL1.pitch(ind_sub_aux).', ...
                                     25, ...
                                     -sigL1.vel1(:, ind_sub_aux), -sigL1.vel3(:, ind_sub_aux), ...
                                     -sigL1.vel4(:, ind_sub_aux), -sigL1.vel2(:, ind_sub_aux), ...
                                     -sigL1.vel5(:, ind_sub_aux), ...
                                     sigL1.cellcenter, lGimbaled, BinmapType, true, luse3beams);
    
            % PS: The column dimension should be the time dimension
            % for all variables (that's why vectores (e.g. pitch)
            % are transposed)
        end

    % If the Signature is not in the 5-beam list
    else
        % Use 4 beams instead
        disp('Using 4 beams')

        %
        for i2 = 1:size(indbreak_rot, 2)

            %
            disp(['--- Transformation in chunk ' num2str(i2) ' out of ' num2str(size(indbreak_rot, 2)) ' ---'])
            %
            ind_sub_aux = indbreak_rot(1, i2) : indbreak_rot(2, i2);
            %
            [sigL1.u(:, ind_sub_aux), ...
             sigL1.v(:, ind_sub_aux), ...
             sigL1.w(:, ind_sub_aux)] = ...
                        janus2earth(sigL1.heading(ind_sub_aux).' - 90, ...
                                    sigL1.roll(ind_sub_aux).', -sigL1.pitch(ind_sub_aux).', ...
                                    25, ...
                                    -sigL1.vel1(:, ind_sub_aux), -sigL1.vel3(:, ind_sub_aux), ...
                                    -sigL1.vel4(:, ind_sub_aux), -sigL1.vel2(:, ind_sub_aux), ...
                                    sigL1.cellcenter, lGimbaled, BinmapType, luse3beams);
        end
    end

    %
    disp('--- Done with coordinate transformation to magnetic ENU. It took: ---')
    toc

    
    %% Rotate horizontal velocity from magnetic to true north

    %
    tic
    disp('--- Rotating from magnetic ENU to local XY coordinate system ---')
    
    %
    disp(['----- Rotating horizontal velocity from magnetic ENU ' ...
          'to local XY coordinate system. Magnetic declination is ' ...
          num2str(sigL1.magdec, '%.2f') ' degrees -----'])

    %
    [sigL1.u, sigL1.v] = ROXSI_uv_ENtoXY(sigL1.u, sigL1.v, sigL1.site, true);

    %
    disp('--- Done with rotation to local XY. It took: ---')
    toc


    %% QC velocity based on backscatter and/or correlation (????)


    %%
    % ------------------------------------------
    % ---- COMPUTE LOW-FREQUENCY QUANTITIES ----
    % ------------------------------------------

    %
    tic
    disp('--- Computing lower frequency fields ---')

    %
    sigL1.averaged.dt = 10*60;    % in seconds
    %
    npts_avg = sigL1.averaged.dt / (1/sigL1.samplingrateHz);

    %
    timesmooth_lims_1 = datetime(sigL1.dtime(1).Year, sigL1.dtime(1).Month, sigL1.dtime(1).Day, ...
                                 sigL1.dtime(1).Hour, (sigL1.dtime(1).Minute + 1 + (sigL1.averaged.dt/60)), 00);
    timesmooth_lims_2 = datetime(sigL1.dtime(end).Year, sigL1.dtime(end).Month, sigL1.dtime(end).Day, ...
                                 sigL1.dtime(end).Hour, (sigL1.dtime(end).Minute - (sigL1.averaged.dt/60)), 00);
    %
    sigL1.averaged.dtime = timesmooth_lims_1 : seconds(sigL1.averaged.dt) : timesmooth_lims_2;
    sigL1.averaged.dtime.TimeZone = sigL1.dtime.TimeZone;

    %
    sigL1.averaged.pressure = time_smooth_reg(sigL1.dtime, sigL1.pressure, ...
                                              sigL1.averaged.dt, ...
                                              sigL1.averaged.dtime([1, end]));
    sigL1.averaged.bottomdepthfrompres = (1e4 * sigL1.averaged.pressure) ./ (1030*9.8);
    

    %
    sigL1.averaged.dtime = sigL1.averaged.dtime(:);
    sigL1.averaged.pressure = sigL1.averaged.pressure(:);
    sigL1.averaged.bottomdepthfrompres = sigL1.averaged.bottomdepthfrompres(:);


    % Find bins below the surface
    ind_avg_abovesurface = 1:1:(length(sigL1.zhab) * length(sigL1.averaged.dtime));
    ind_avg_abovesurface = reshape(ind_avg_abovesurface, length(sigL1.zhab), length(sigL1.averaged.dtime));
    %
    for i2 = 1:length(sigL1.averaged.dtime)
        %
        ind_withinocean_aux = find(((sigL1.zhab + sigL1.binsize) < sigL1.averaged.bottomdepthfrompres(i2)), 1, 'last');
        %
        ind_avg_abovesurface(1:(ind_withinocean_aux-1), i2) = NaN;    % and remove another bin because of waves
    end
    %
    ind_avg_abovesurface = ind_avg_abovesurface(:);
    ind_avg_abovesurface = ind_avg_abovesurface(~isnan(ind_avg_abovesurface));

    % Now average the beam data
    prealloc_aux = NaN(length(sigL1.zhab), length(sigL1.averaged.dtime));
    %
    sigL1.averaged.u = prealloc_aux;
    sigL1.averaged.v = prealloc_aux;
    sigL1.averaged.w = prealloc_aux;
    %
% %     sigL1.averaged.amp1 = prealloc_aux;
% %     sigL1.averaged.amp2 = prealloc_aux;
% %     sigL1.averaged.amp3 = prealloc_aux;
% %     sigL1.averaged.amp4 = prealloc_aux;

    %
% %     list_fields_aux = {'u', 'v', 'w', 'amp1', 'amp2', 'amp3', 'amp4'};
    list_fields_aux = {'u', 'v', 'w'};

    %
    for i2 = 1:length(list_fields_aux)
        for i3 = 1:length(sigL1.zhab)
            %
            sigL1.averaged.(list_fields_aux{i2})(i3, :) = ...
                             time_smooth_reg(sigL1.dtime, ...
                                             sigL1.(list_fields_aux{i2})(i3, :), ...
                                             sigL1.averaged.dt, ...
                                             sigL1.averaged.dtime([1, end]));
        end

        % Remove the above the surface
        sigL1.averaged.(list_fields_aux{i2})(ind_avg_abovesurface) = NaN;
        
    end

    %
    disp('--- Done with computing lower frequency quantities ---')
    toc



    %%
    % -----------------------------------------------------------
    % -------- FINAL ADJUSTMENTS, SAVE DATA AND QC PLOTS --------
    % -----------------------------------------------------------
    % ----------------------------------------------------


    %% Add README

    %
    time_dataproc = datetime('now', 'TimeZone', 'Local');
    time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');
    % Add a general README
    sigL1.README = ['Level 1 Signature1000 data from ROXSI 2022. The data is from Signature ' ...
                    'with serial number SN and deployed at mooring site mooringID. ' ...
                    'Data processed by script ' mfilename() '.m on ' time_dataproc_char ' (TimeZone ' time_dataproc.TimeZone '). ' ...
                    'Horizontal velocity components are relative to local XY coordinate system, and the magnetic ' ...
                    'decliation was taken from www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml. Data ' ...
                    'in the Level 1 structure has been trimmed for the deployment ' ...
                    'period, as defined in the table deploymentInfo_ROXSI2022.mat. ' ...
                    'Pressure is in dbar, where atmospheric pressure has been removed.'];


    %% Move along-beam, backscatter, and correlation to a separate structure

    %
    sigL1beamdata.SN = sigL1.SN;
    sigL1beamdata.mooringID = sigL1.mooringID;
    %
    sigL1beamdata.dtime = sigL1.dtime;
    sigL1beamdata.zhab = sigL1.zhab;

    %
    list_vars_move = {'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
                      'amp1', 'amp2', 'amp3', 'amp4', 'amp5', ...
                      'cor1', 'cor2', 'cor3', 'cor4', 'cor5'};

% %     %
% %     for i2 = 1:length(list_vars_move)
% %         %
% %         sigL1beamdata.SN = sigL1.SN;
% %         sigL1beamdata.mooringID = sigL1.mooringID;
% %         %
% %         sigL1beamdata.dtime = sigL1.dtime;
% %         sigL1beamdata.zhab = sigL1.zhab;
% %         %
% %         if isfield(sigL1, list_vars_move{i2})
% %             %
% %             sigL1beamdata.(list_vars_move{i2}) = sigL1.(list_vars_move{i2});
% %             %
% %             sigL1 = rmfield(sigL1, list_vars_move{i2});
% %         end
% %     end


    %% Save stuff

    % ----------------------------------------------------
    % Save QC figures

    % ----------------------------------------------------
    % Save level 1 data along-beam/backscatter/correlation quantities

    %
    disp('----- Saving beam data in level 1 data file -----')
    str_filename = ['roxsi_signature_L1_' char(sigL1.mooringID) '_' char(sigL1.SN) '_alongbeamvel'];
    %
    save(fullfile(dir_output_L1, [str_filename '.mat']), 'sigL1beamdata', '-v7.3')
    %
    if sigL1.l5beams
        save(fullfile(dir_output_L1, [str_filename '.mat']), '-struct', 'sigL1', 'pressure', 'vel1', 'vel2', 'vel3', 'vel4', 'vel5', '-append')
    else
        save(fullfile(dir_output_L1, [str_filename '.mat']), '-struct', 'sigL1', 'pressure', 'vel1', 'vel2', 'vel3', 'vel4', '-append')
    end

    % ----------------------------------------------------
    % Save level 1 data with all variables
    %
    for i2 = 1:length(list_vars_move)
        if isfield(sigL1, list_vars_move{i2})
            sigL1 = rmfield(sigL1, list_vars_move{i2});
        end
    end
    %
    disp('----- Saving primary level 1 data structure -----')
    str_filename = ['roxsi_signature_L1_' char(sigL1.mooringID) '_' char(sigL1.SN)];
    %
    save(fullfile(dir_output_L1, [str_filename '.mat']), 'sigL1', '-v7.3')


    % ----------------------------------------------------
    % Save level 1 data with scalar variables only
    %
    sigL1.cellcenter = sigL1.cellcenter(1);
    sigL1.zhab = sigL1.zhab(1);
    %
    sigL1.u = sigL1.u(1, :);    sigL1.u = sigL1.u(:);
    sigL1.v = sigL1.v(1, :);    sigL1.v = sigL1.v(:);
    sigL1.w = sigL1.w(1, :);    sigL1.w = sigL1.w(:);
    %
    sigL1.averaged.u = sigL1.averaged.u(1, :);    sigL1.averaged.u = sigL1.averaged.u(:);
    sigL1.averaged.v = sigL1.averaged.v(1, :);    sigL1.averaged.v = sigL1.averaged.v(:);
    sigL1.averaged.w = sigL1.averaged.w(1, :);    sigL1.averaged.w = sigL1.averaged.w(:);

    %
    disp('----- Saving level 1 data with scalars only -----')
    str_filename = ['roxsi_signature_L1_' char(sigL1.mooringID) '_' char(sigL1.SN) '_scalars'];
    %
    save(fullfile(dir_output_L1, [str_filename '.mat']), 'sigL1', '-v7.3')



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

