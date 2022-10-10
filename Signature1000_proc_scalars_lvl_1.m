%% Script that processes (from RAW to level 1) the scalar variales
% from Signatures 1000 deployed in ROXSI's large-scale array. This
% is a separate script from velocity processing because pressure
% measurements may be useful at times when velocity measurements
% are bad (e.g. due to large tilt).
%
% While ADCP heading/tilt is not necessarily useful for analysis
% of pressure data, it gives an idea of the ADCP deployment.
% 
%
% All Signatures (ADCPs) were programmed to start
% sampling at 2022/06/21 18:00:00 (local). Signature
% SN 103045, at B10, was only deployed later due to
% weather conditions.

clear
close all


%%
% --------------------------------------
% -------------- PREAMBLE --------------
% --------------------------------------

%%

% % roxsi_add_libraries()

%%

%
% % dirparent_data = data_dirpath();
dirparent_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
%
dir_rawdata_parent = fullfile(dirparent_data, 'RAW', 'Signature1000');


%%

%
% % dir_output_parent = data_dirpath();
% % dir_output_parent = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/';
dir_output_parent = '/home/omarques/Documents/MATLAB/ROXSIproc_output/';
%
% dir_output_data_L1 = fullfile(dir_output_parent, 'Level1_Data', 'Signature1000_Level1');
% dir_output_figs_L1 = fullfile(dir_output_parent, 'Level1_Data', 'Signature1000_Level1', 'qc_plots');
dir_output_data_L1 = dir_output_parent;

% Logical switches to save or not save data and figures
lsave_file = true;
lsave_fig = true;


% % %% Define time limits and time step for processing velocity data. Since
% % % there are a lot of Signature1000, it seems appropriate to break
% % % it apart in smaller chunks of time. Roughly speaking,
% % % files with 1/4 of a day (6 hours) will be around 300 MB.
% % %
% % % Note that the pressure (and heading/tilt/pitch/temperature) are
% % % still loaded for the full timeseries (for general QC) and only
% % % the processing of velocity will be within these time limits.
% % %
% % % Note that in the definition I will use, each file will have data
% % % starting and including the first time limit of each segment
% % % (given by time_lims_proc) and ending just before (not including)
% % % the next time limit.
% % %
% % % Note these time limits are "real times" -- i.e. they should be
% % % compared to clock-drift-corrected timestamps.
% % 
% % % -----------------------------------
% % % % %
% % % % time_beginend_proc = [datetime(2022, 06, 29, 00, 00, 00), ...
% % % %                       datetime(2022, 06, 29, 06, 00, 00)];
% % % % time_beginend_proc.TimeZone = 'America/Los_Angeles';
% % % % 
% % % % %
% % % % dtstep_proc = hours(6);
% % % % 
% % % % %
% % % % time_lims_proc = time_beginend_proc(1) : dtstep_proc : time_beginend_proc(2);
% % 
% % % -----------------------------------
% % %
% % % % time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 06, 00, 00); ...
% % % %                   datetime(2022, 07, 03, 18, 00, 00), datetime(2022, 07, 04, 00, 00, 00)];
% % %
% % % % time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 06, 00, 00)];    % 6 hour chunk
% % %
% % time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 07, 08, 00, 00, 00)];   % around 10 days
% % %
% % % time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 00, 05, 00)];   % 5 minutes
% % % time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 00, 01, 00)];   % 1 minute
% % time_lims_proc.TimeZone = 'America/Los_Angeles';
% % 
% % %
% % Ndatasegments = size(time_lims_proc, 1);


%% Load ADCP deployment information

% Just to be clear: file and variable have
% the same name (though not a requirement)
%
% dir_coderepo = repo_dirpath();
dir_coderepo = '/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc';
load(fullfile(dir_coderepo, 'deploymentInfo_ROXSI2022.mat'), 'deploymentInfo_ROXSI2022')


%% List of Signatures that will be processed

% All Signatures
list_Signature = {'A01_103043', ...
                  'B10_103045', ...
                  'B13_103046', ...
                  'B15_103056', ...
                  'B17_101923', ...
                  'C01_102128', ...
                  'X05_100231', ...
                  'X11_101941'};

% % % Just a test
list_Signature = {'A01_103043'};

%
Nsignatures = length(list_Signature);


% % %% The list of ADCPs processed with either 4 or 5 beams
% % 
% % %
% % % % list_4beams = {'B10', 'B13', 'B15', 'B17', 'X11'};    % B10 has 5th beam, but on HR mode (up to 4 meters above the transducer). Others did the echosounder.
% % % % list_5beams = {'A01', 'C01', 'X05'};
% % %
% % list_4beams = {'B10'};    % B10 has 5th beam, but on HR mode (up to 4 meters above the transducer). Others did the echosounder.
% % list_5beams = {'A01', 'B13', 'B15', 'B17', 'C01', 'X05', 'X11'};
% % 
% % % Check all of the list of Signatures is in either of the
% % % lists for processing with 4 or 5 beams
% % for i = 1:length(list_Signature)
% % 
% %     %
% %     linanylist = any(contains(list_4beams, list_Signature{i}(1:3))) || ...
% %                  any(contains(list_5beams, list_Signature{i}(1:3)));
% % 
% %     %
% %     if ~linanylist
% % 
% %         error(['Signature ' list_Signature{i} ' is not present in ' ...
% %                'any of the lists specifying velocity processing.'])
% % 
% %     end
% % 
% % end

    

%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------------- DO DATA PROCESSING ----------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------


%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_Signature_procL1_scalars_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_output_data_L1, log_file_name))

%
totalRunTime = tic;

% % 
% % %% Fields to be deleters

%%
% % 
% % list_Data_fields = {'Burst_Time', ...
% %                     'Burst_VelBeam1', 'Burst_VelBeam2', 'Burst_VelBeam3', 'Burst_VelBeam4', ...
% %                     'Burst_AmpBeam1', 'Burst_AmpBeam2', 'Burst_AmpBeam3', 'Burst_AmpBeam4', ...
% %                     'Burst_Heading', 'Burst_Pitch', 'Burst_Roll', ...
% %                     'Burst_Temperature', 'Burst_Pressure'};


%% Display on the screen:

%
disp(' '), disp(' ')
disp('------------------------------ Processing scalar variables from Signature1000s ------------------------------')
disp('List of Signature1000s being processed:')
%
for i = 1:Nsignatures
    disp([num2str(i) ' - ' list_Signature{i}])
end



% Loop over Signature1000's in the list
for i1 = 1:Nsignatures


    %%
    %
    disp(' '), disp(' ')
    disp(['----- Start processing raw Signature1000 data from: ' list_Signature{i1} ' -----'])


    %% Get the Signature1000 files in the correct order

    %
    dir_data_aux = fullfile(dir_rawdata_parent, list_Signature{i1}, 'converted');
    %
    list_dir_aux = Signature_orderfiles(dir_data_aux, list_Signature{i1}(1:3));

    %
    Nfiles = length(list_dir_aux);


    %% First add basic metadata to structure variable

    %
    sig1000.SN = convertCharsToStrings(list_Signature{i1}(5:end));
    sig1000.mooringID = convertCharsToStrings(list_Signature{i1}(1:3));

    % Latitude/longitude
    info_mooringtable = ROXSI_mooringlocation(sig1000.mooringID, "ADCP");
    %
    sig1000.latitude = info_mooringtable.latitude;
    sig1000.longitude = info_mooringtable.longitude;


    %% Get trimming times for the deployment period of this Signature

    %
    ind_row_match = find(strcmp(deploymentInfo_ROXSI2022.SN, list_Signature{i1}(5:end)));
    
    %
    time_1 = deploymentInfo_ROXSI2022.time_begin_trim(ind_row_match);
    time_2 = deploymentInfo_ROXSI2022.time_end_trim(ind_row_match);

    %
    time_1 = datetime(datenum(time_1, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    time_2 = datetime(datenum(time_2, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');

% %     %
% %     lin_deployment = (senAQDP_aux.time >= time_1) & (senAQDP_aux.time <= time_2);


    %%
    % ---------------------------------------------------------------
    % -------- NOW GET HEADING AND TILT AND DO FIRST QC PLOT --------
    % ---------------------------------------------------------------

    %% Load Configuration from first file

    % Just the Config structure variable
    dataread_aux = load(fullfile(dir_data_aux, list_dir_aux(1).name), 'Config');

    %
    sig1000.Config = dataread_aux.Config;
    

    %% First get the timeseries of scalar variables (these are the
    % "Burst" variables, which are given at the same time as the
    % data from the 4 beams. Signature1000 also gives these variables
    % for the 5th beam "IBurst" variables)

    % Pre-allocate in cell arrays
    prealloc_aux = cell(1, Nfiles);
    %
    sig1000.timedatenum = prealloc_aux;
    sig1000.dtime = prealloc_aux;    % variable not filled  in the loop
    sig1000.pressure = prealloc_aux;
    sig1000.temperature = prealloc_aux;
    %
    sig1000.heading = prealloc_aux;
    sig1000.pitch = prealloc_aux;
    sig1000.roll = prealloc_aux;

    %
    disp('--- Load all scalar variables ---')

    tic
    % Loop over files
    for i2 = 1:Nfiles

        %% Load data in file i2'th
        %
        dataread_aux = load(fullfile(dir_data_aux, list_dir_aux(i2).name));

% %         %% Delete variables...
% % 
% %         %
% %         dataread_aux = rmfield(dataread_aux, {'Units', 'Descriptions'});
% %         %
% %         list_fields_in_Data = fieldnames(dataread_aux.Data);
% %         %
% %         for i3 = 1:length(list_fields_in_Data)
% %             %
% %             if ~any(strcmp(list_fields_in_Data{i3}, list_Data_fields))
% %                 dataread_aux.Data = rmfield(dataread_aux.Data, list_fields_in_Data{i3});
% %             end
% %         end


        %% Put time, pressure (which is already in dbar),
        % heading, pitch, and roll in output data structure

        %
        sig1000.timedatenum{i2} = dataread_aux.Data.Burst_Time;
        sig1000.pressure{i2} = dataread_aux.Data.Burst_Pressure;
        sig1000.temperature{i2} = dataread_aux.Data.Burst_Temperature;
        %
        sig1000.heading{i2} = dataread_aux.Data.Burst_Heading;
        sig1000.pitch{i2} = dataread_aux.Data.Burst_Pitch;
        sig1000.roll{i2} = dataread_aux.Data.Burst_Roll;

        %%
        clear dataread_aux

    end
    disp('--- Done with loading all scalar variables ---')
    toc

    % Concatenate cell array into a long column vector
    sig1000.timedatenum = cat(1, sig1000.timedatenum{:});
    sig1000.pressure = cat(1, sig1000.pressure{:});
    sig1000.temperature = cat(1, sig1000.temperature{:});
    %
    sig1000.heading = cat(1, sig1000.heading{:});
    sig1000.pitch = cat(1, sig1000.pitch{:});
    sig1000.roll = cat(1, sig1000.roll{:});


    disp('--- Done getting timeseries of scalar variables ---')


    %% Do clock drift correction

    %
    disp('----- Correcting for clock drift -----')

    % Now correct for clock drift
    time_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, ...
                                            list_Signature{i1}(5:end), ...
                                            sig1000.timedatenum);


    %% Convert time to date time

    %
    sig1000.dtime = datetime(time_aux, 'ConvertFrom', 'datenum', ...
                                       'TimeZone', 'America/Los_Angeles');


    %% Apply (small) correction due to atmospheric pressure variability.
    % 
    % As opposed to Aquadopps, Signatures don't have
    % a customizable pressure offset.




    %% Find all bins that are entirely below the maximum bottom depth

% % % %     %
% % % %     lin_verticalrange = ((sig1000.zhab + (sig1000.binsize/2)) < ...
% % % %                          max(sig1000.bottomdepthfrompres));
% %     % Or dummy for code developing purposes
% %     lin_verticalrange = true(1, length(sig1000.zhab));
% % 
% %     %
% %     sig1000.zhab = sig1000.zhab(lin_verticalrange);


    %%
    % ---------------------------------------------------------------
    % --------------- PLOT PRESSURE, HEADING, AND TILT --------------
    % ---------------------------------------------------------------

    %%

% %     %
% %     disp('--- QC plot with timeseries of scalar variables ---')
% % 
% %     %
% %     fig_L1_QC_tilt = figure;
% %     set(fig_L1_QC_tilt, 'units', 'normalized')
% %     set(fig_L1_QC_tilt, 'Position', [0.2, 0.2, 0.4, 0.6])
% %         %
% %         haxs_1 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.73, 0.8, 0.17]);
% %         haxs_2 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.52, 0.8, 0.17]);
% %         haxs_3 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.31, 0.8, 0.17]);
% %         haxs_4 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.10, 0.8, 0.17]);
% %         %
% %         haxs_all = [haxs_1, haxs_2, haxs_3, haxs_4];
% %         hold(haxs_all, 'on')
% %         %
% %         plot(haxs_1, sig1000.dtime, sig1000.pressure, '-k')
% %         plot(haxs_2, sig1000.dtime, sig1000.heading, '-k')
% %         plot(haxs_3, sig1000.dtime, sig1000.pitch, '-k')
% %         plot(haxs_4, sig1000.dtime, sig1000.roll, '-k')
% % 
% %     %
% %     set(haxs_all, 'FontSize', 16, 'Box', 'on', ...
% %                   'XGrid', 'on', 'YGrid', 'on')
% %     %
% %     set(haxs_all, 'XLim', sig1000.dtime([1, end]) + [-hours(12); hours(12)])
% %     %
% %     lin_deployment = (sig1000.dtime >= time_1) & (sig1000.dtime <= time_2);
% %     %
% %     pres_sub_aux = sig1000.pressure(lin_deployment);
% %     ylim(haxs_1, [min(pres_sub_aux), max(pres_sub_aux)])
% % % % %     ylim(haxs_2, [0, 360])
% % % %     ylim(haxs_3, [-6, 6])
% % % %     ylim(haxs_4, [-6, 6])
% %     %
% %     heading_sub_aux = sig1000.heading(lin_deployment);
% %     pitch_sub_aux = sig1000.pitch(lin_deployment);
% %     roll_sub_aux = sig1000.roll(lin_deployment);
% %     ylim(haxs_2, [min(heading_sub_aux), max(heading_sub_aux)])
% %     ylim(haxs_3, [min(pitch_sub_aux), max(pitch_sub_aux)])
% %     ylim(haxs_4, [min(roll_sub_aux), max(roll_sub_aux)])
% % 
% %     %
% %     ylabel(haxs_1, '[dbar]', 'Interpreter', 'Latex', 'FontSize', 16)
% %     ylabel(haxs_2, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
% %     ylabel(haxs_3, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
% %     ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
% %     %
% %     title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
% %                    char(sig1000.SN) ': pressure, heading, pitch, and roll'], ...
% %                   'Interpreter', 'Latex', 'FontSize', 16)
% %     %
% %     linkaxes([haxs_1, haxs_2, haxs_3, haxs_4], 'x')
% % 
% % 
% %     %
% %     for i2 = 1:length(haxs_all)
% %         ylims_aux = ylim(haxs_all(i2));
% %         %
% %         plot(haxs_all(i2), [time_1, time_1], ylims_aux, '--r')
% %         plot(haxs_all(i2), [time_2, time_2], ylims_aux, '--r')
% %         %
% %         ylim(haxs_all(i2), ylims_aux)
% %     end
% % 
% % 
% %     %
% %     exportgraphics(fig_L1_QC_tilt, fullfile(pwd, ['sig_pres_tilt_' list_Signature{i1} '.png']), 'Resolution', 300)


    %% Includes temperature and plot pitch and roll on same subplot

    %
    disp('--- QC plot with timeseries of scalar variables ---')

    %
    fig_L1_QC_tilt = figure;
    set(fig_L1_QC_tilt, 'units', 'normalized')
    set(fig_L1_QC_tilt, 'Position', [0.2, 0.2, 0.4, 0.6])
        %
        haxs_1 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.73, 0.8, 0.17]);
        haxs_2 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.52, 0.8, 0.17]);
        haxs_3 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.31, 0.8, 0.17]);
        haxs_4 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.10, 0.8, 0.17]);
        %
        haxs_all = [haxs_1, haxs_2, haxs_3, haxs_4];
        hold(haxs_all, 'on')
        %
        plot(haxs_1, sig1000.dtime, sig1000.pressure, '-k')
        plot(haxs_2, sig1000.dtime, sig1000.temperature, '-k')
        plot(haxs_3, sig1000.dtime, sig1000.heading, '-k')
        %
        plot(haxs_4, sig1000.dtime, sig1000.pitch, '-b')
        plot(haxs_4, sig1000.dtime, sig1000.roll, '-r')

    %
    set(haxs_all, 'FontSize', 16, 'Box', 'on', ...
                  'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_all, 'XLim', sig1000.dtime([1, end]) + [-hours(12); hours(12)])
    %
    lin_deployment = (sig1000.dtime >= time_1) & (sig1000.dtime <= time_2);
    %
    pres_sub_aux = sig1000.pressure(lin_deployment);
    ylim(haxs_1, [min(pres_sub_aux), max(pres_sub_aux)])
% % %     ylim(haxs_2, [0, 360])
% %     ylim(haxs_3, [-6, 6])
% %     ylim(haxs_4, [-6, 6])
    %
    temperature_sub_aux = sig1000.temperature(lin_deployment);
    heading_sub_aux = sig1000.heading(lin_deployment);
    pitch_sub_aux = sig1000.pitch(lin_deployment);
    roll_sub_aux = sig1000.roll(lin_deployment);
    %
    ylim(haxs_2, [min(temperature_sub_aux), max(temperature_sub_aux)])
    ylim(haxs_3, [min(heading_sub_aux), max(heading_sub_aux)])
    ylim(haxs_4, [min([min(pitch_sub_aux), min(roll_sub_aux)]), ...
                  max([max(pitch_sub_aux), max(roll_sub_aux)])])

    %
    ylabel(haxs_1, '[dbar]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, '[$^\circ$C]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_3, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
                   char(sig1000.SN) ': pressure, temperature, heading, ' ...
                   'pitch (blue), and roll (red)'], ...
                  'Interpreter', 'Latex', 'FontSize', 16)
    %
    linkaxes([haxs_1, haxs_2, haxs_3, haxs_4], 'x')


    %
    for i2 = 1:length(haxs_all)
        ylims_aux = ylim(haxs_all(i2));
        %
        plot(haxs_all(i2), [time_1, time_1], ylims_aux, '--r')
        plot(haxs_all(i2), [time_2, time_2], ylims_aux, '--r')
        %
        ylim(haxs_all(i2), ylims_aux)
    end


    %
    exportgraphics(fig_L1_QC_tilt, fullfile(pwd, ['sig_scalars_' list_Signature{i1} '.png']), 'Resolution', 300)


    %% Trim variables for the desired period where analysis will be done
    %
    % The logical variable lin_deployment is created in the code block
    % above.

    %
    list_variables_aux = fieldnames(sig1000);

    %
    for i2 = 1:length(list_variables_aux)
        %
        if length(sig1000.(list_variables_aux{i2}))==length(lin_deployment)
            %
            sig1000.(list_variables_aux{i2}) = sig1000.(list_variables_aux{i2})(lin_deployment);
        end
    end


    %% Check clock

    %
    disp('--- QC plot checking diff time ---')

    %
    inds_time = 1:length(sig1000.dtime);
    inds_difftime = (inds_time(1:end-1) + inds_time(2:end))./2;

    %
    fig_L1_QC_clock = figure;
    set(fig_L1_QC_clock, 'units', 'normalized')
    set(fig_L1_QC_clock, 'Position', [0.2, 0.2, 0.4, 0.6])
        %
        haxs_1 = axes(fig_L1_QC_clock, 'Position', [0.1, 0.575, 0.8, 0.325]);
        haxs_2 = axes(fig_L1_QC_clock, 'Position', [0.1, 0.150, 0.8, 0.325]);
        %
        haxs_all = [haxs_1, haxs_2];
        hold(haxs_all, 'on')
        %
        plot(haxs_1, inds_time, sig1000.dtime, '-k')
        plot(haxs_2, inds_difftime, seconds(sig1000.dtime), '-k')

    %
    set(haxs_all, 'FontSize', 16, 'Box', 'on', ...
                  'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_all, 'XLim', sig1000.dtime([1, end]) + [-hours(12); hours(12)])
    %
    ylim(haxs_1, sig1000.dtime([1, end]))

    %
    ylabel(haxs_1, 'Time', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, 'seconds', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
                   char(sig1000.SN) ': time and diff(time) (in seconds)'], ...
                  'Interpreter', 'Latex', 'FontSize', 16)
    %
    linkaxes([haxs_1, haxs_2, haxs_3, haxs_4], 'x')


    %
    for i2 = 1:length(haxs_all)
        ylims_aux = ylim(haxs_all(i2));
        %
        plot(haxs_all(i2), [time_1, time_1], ylims_aux, '--r')
        plot(haxs_all(i2), [time_2, time_2], ylims_aux, '--r')
        %
        ylim(haxs_all(i2), ylims_aux)
    end


    %
    exportgraphics(fig_L1_QC_tilt, fullfile(pwd, ['sig_clock_' list_Signature{i1} '.png']), 'Resolution', 300)


    %% Interpolate variables to gridded time vector (after
    % making sure there are no major issues above)

% %     %
% % %     dtime_grid = 
% %     %
% %     Nlengthtimeseries = length(sig1000.dtime);
% %     %
% % %     list_variables_interp = {};
% %     %
% %     for i2 = 1:length(list_variables_aux)
% %         %
% %         if length(sig1000.(list_variables_aux{i2})) == Nlengthtimeseries
% %             %
% % %             interp1(sig1000.dtime, sig1000.(list_variables_aux{i2}), dtime_grid)
% % % %              = sig1000.(list_variables_aux{i2})(lin_deployment);
% %         end
% %     end
% % % 
% %   
% %     % Now replace the original time by gridded time
% % % %     sig1000.dtime = dtime_grid;
% % 
% % % %  Remove other variables

%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ------------------------- FINAL STEPS --------------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------

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
return

    %%


    % ----------------------------------------------------
    % Save level 1 data and QC figures

    %
    if lsave_file
        %
        disp('----- Saving level 1 data -----')

        %
        str_filename = ['roxsi_signature_L1_' char(sig1000.mooringID) '_' char(sig1000.SN)];
        %
        save(fullfile(dir_output_data_L1, [str_filename '.mat']), 'sig1000', '-v7.3')
    end

    %
    if lsave_fig
% %         %
% %         disp('----- Saving level 1 QC plot figures -----')
% % 
% %         %
% %         str_filename = ['roxsi_aquadopp_L1_' list_Signature{i1}(1:3) '_' char(sig1000.SN) '_QC_1'];
% %         % Save figure as *.png
% %         exportgraphics(fig_L1_QC_tilt, fullfile(dir_output_figs_L1, [str_filename '.png']), 'Resolution', 300)
% % 
% %         %
% %         str_filename = ['roxsi_aquadopp_L1_' list_Signature{i1}(1:3) '_' char(sig1000.SN) '_QC_2'];
% %         % Save figure as *.png
% %         exportgraphics(fig_L1_QC, fullfile(dir_output_figs_L1, [str_filename '.png']), 'Resolution', 300)
    
    end


    %%

% %     % Clear some variables to avoid issues in the next loop iteration
% %     close(fig_L1_QC_tilt), close(fig_L1_QC)
% %     clear aquadoppL1 header_aux beamAQDP_aux senAQDP_aux
    %
    clear sig1000


    %%
    
    %
    disp(['----- DONE with Signature data proc: ' list_Signature{i1} ' -----'])
    toc(totalRunTime)


end

%
disp('###################### Done with processing scalar variables for all Signature1000 ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');
