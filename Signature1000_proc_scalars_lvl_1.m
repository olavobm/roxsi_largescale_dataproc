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

%
% % dirparent_data = data_dirpath();
dirparent_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
%
dir_rawdata_parent = fullfile(dirparent_data, 'RAW', 'Signature1000');


%%

%
% % dir_output_parent = data_dirpath();
% % dir_output_parent = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/';
% % dir_output_parent = '/home/omarques/Documents/MATLAB/ROXSIproc_output/';
dir_output_parent = pwd;
%
% dir_output_data_L1 = fullfile(dir_output_parent, 'Level1_Data', 'Signature1000_Level1');
% dir_output_figs_L1 = fullfile(dir_output_parent, 'Level1_Data', 'Signature1000_Level1', 'qc_plots');
dir_output_data_L1 = dir_output_parent;

% Logical switches to save or not save data and figures
lsave_file = true;
lsave_fig = true;


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
                  'X11_101941', ...
                  'X05_100231'};

% Just a test
list_Signature = {'A01_103043'};

%
Nsignatures = length(list_Signature);


%% Set whether the first velocity bin will be included
% (no QC on this velocity)

%
% lvelbin1 = false;
lvelbin1 = true;

%
if lvelbin1
    roxsi_add_libraries()    % add ADCPtools
end


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

    %
    sig1000.site = info_mooringtable.roxsiarray;


    %% Get trimming times for the deployment period of this Signature

    %
    ind_row_match = find(strcmp(deploymentInfo_ROXSI2022.SN, list_Signature{i1}(5:end)));
    
    %
    time_1 = deploymentInfo_ROXSI2022.time_begin_trim(ind_row_match);
    time_2 = deploymentInfo_ROXSI2022.time_end_trim(ind_row_match);

    %
    time_1 = datetime(datenum(time_1, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    time_2 = datetime(datenum(time_2, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');


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
    sig1000.samplingrateHz = NaN;    % only added at the end
    sig1000.dtime = prealloc_aux;    % variable not filled  in the loop
    sig1000.pressure = prealloc_aux;
    sig1000.temperature = prealloc_aux;
    %
    sig1000.heading = prealloc_aux;
    sig1000.pitch = prealloc_aux;
    sig1000.roll = prealloc_aux;

    %
    if lvelbin1
        %
        disp('--- Adding velocity from 1st bin ---')
        %
        sig1000.vel1 = prealloc_aux;
        sig1000.vel2 = prealloc_aux;
        sig1000.vel3 = prealloc_aux;
        sig1000.vel4 = prealloc_aux;
        %
        sig1000.cor1 = prealloc_aux;
        sig1000.cor2 = prealloc_aux;
        sig1000.cor3 = prealloc_aux;
        sig1000.cor4 = prealloc_aux;
    end

    %
    disp('--- Load all scalar variables ---')

    tic
    % Loop over files
    for i2 = 1:Nfiles

        %% Load data in file i2'th
        %
        dataread_aux = load(fullfile(dir_data_aux, list_dir_aux(i2).name));


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

        %
        if lvelbin1
            % Take only the first bin (or any other, but just one bin)
            sig1000.vel1{i2} = dataread_aux.Data.Burst_VelBeam1(:, 1);
            sig1000.vel2{i2} = dataread_aux.Data.Burst_VelBeam2(:, 1);
            sig1000.vel3{i2} = dataread_aux.Data.Burst_VelBeam3(:, 1);
            sig1000.vel4{i2} = dataread_aux.Data.Burst_VelBeam4(:, 1);
            %
            sig1000.cor1{i2} = dataread_aux.Data.Burst_CorBeam1(:, 1);
            sig1000.cor2{i2} = dataread_aux.Data.Burst_CorBeam2(:, 1);
            sig1000.cor3{i2} = dataread_aux.Data.Burst_CorBeam3(:, 1);
            sig1000.cor4{i2} = dataread_aux.Data.Burst_CorBeam4(:, 1);
        end


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

    %
    if lvelbin1
        sig1000.vel1 = cat(1, sig1000.vel1{:});
        sig1000.vel2 = cat(1, sig1000.vel2{:});
        sig1000.vel3 = cat(1, sig1000.vel3{:});
        sig1000.vel4 = cat(1, sig1000.vel4{:});
        %
        sig1000.cor1 = double(cat(1, sig1000.cor1{:}));
        sig1000.cor2 = double(cat(1, sig1000.cor2{:}));
        sig1000.cor3 = double(cat(1, sig1000.cor3{:}));
        sig1000.cor4 = double(cat(1, sig1000.cor4{:}));
    end


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
    % As opposed to Aquadopps, Signatures don't have a customizable
    % pressure offset (so pressure measurements don't include the
    % full pressure at the sensor)

    % Interpolate atmospheric pressure anomaly to timestamps
    % of the signature
    atmpresanomaly_aux = interp1(atmpresanomaly.dtime, ...
                                 atmpresanomaly.atm_anomaly, sig1000.dtime);

    %
    sig1000.pressure = sig1000.pressure - atmpresanomaly_aux;


    %%
    % ---------------------------------------------------------------
    % ----- PLOT PRESSURE, TEMPERATURE, HEADING, PITCH, AND ROLL ----
    % ---------------------------------------------------------------

    %% Older version (without temperature)

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


    %% Similar as above, but includes temperature
    % and plot pitch and roll on same subplot

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
    ylim(haxs_3, [min(heading_sub_aux), max(heading_sub_aux)] + [-1, +1])
    ylim(haxs_4, [min([min(pitch_sub_aux), min(roll_sub_aux)]), ...
                  max([max(pitch_sub_aux), max(roll_sub_aux)])] + [-1, +1])

    %
    ylabel(haxs_1, '[dbar]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, '[$^\circ$C]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_3, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
                   char(sig1000.SN) ': pressure, temperature, heading, ' ...
                   'pitch (blue), and roll (red)'], ...
                  'Interpreter', 'Latex', 'FontSize', 12)
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

% %     % Saved at the end of the for loop
% %     exportgraphics(fig_L1_QC_tilt, fullfile(pwd, ['sig_scalars_' list_Signature{i1} '.png']), 'Resolution', 300)


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

    % Check for something that should never happen at this point
    if any(isnan(sig1000.timedatenum))
        warning(['###### Signature ' list_Signature{i1} ' has ' ...
                 'invalid (NaN) timestamps ######'])
    end

    % Now do a plot for cheching diff(time)
    disp('--- QC plot checking diff time ---')

    %
    inds_time = 1:length(sig1000.dtime);
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
        plot(haxs_1, inds_time, sig1000.dtime, '-k')
        plot(haxs_2, inds_difftime, seconds(diff(sig1000.dtime)), '-k')

    %
    set(haxs_all, 'FontSize', 12, 'Box', 'on', ...
                  'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_all, 'XLim', [0, (inds_time(end) + 1)])
    %
    ylim(haxs_1, sig1000.dtime([1, end]))

    %
    xlabel(haxs_2, 'Indices', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    ylabel(haxs_1, 'Time', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, 'seconds', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
                   char(sig1000.SN) ': time and diff(time) (in seconds)'], ...
                  'Interpreter', 'Latex', 'FontSize', 16)
    %
    linkaxes(haxs_all, 'x')


    % Plot horizontal lines for trimming edges
    xlims_aux = xlim(haxs_all(1));
    %
    plot(haxs_all(1), xlims_aux, [time_1, time_1], '--r')
    plot(haxs_all(1), xlims_aux, [time_2, time_2], '--r')
    %
    xlim(haxs_all(1), xlims_aux)
    

% %     % Save at the end of the loop
% %     exportgraphics(fig_L1_QC_clock, fullfile(pwd, ['sig_clock_' list_Signature{i1} '.png']), 'Resolution', 300)


    %% Interpolate variables to gridded time vector (after
    % making sure there are no major issues above)

    %
    disp('--- Gridding scalar variables to a time grid ---')

    % As time edges of the grid, round to the whole minute (done in a more
    % complicated way so that it doesn't rely on Matlab versions newer than
    % at least 2021b)
    dtime_edge_1 = datetime(sig1000.dtime(1).Year, sig1000.dtime(1).Month, sig1000.dtime(1).Day, ...
                            sig1000.dtime(1).Hour, (sig1000.dtime(1).Minute + 1), 00);
    dtime_edge_2 = datetime(sig1000.dtime(end).Year, sig1000.dtime(end).Month, sig1000.dtime(end).Day, ...
                            sig1000.dtime(end).Hour, sig1000.dtime(end).Minute, 00);
    %
    df_sampling = double(sig1000.Config.Burst_SamplingRate);    % in Hertz

    %
    dtime_grid = dtime_edge_1 : seconds(1/df_sampling) : dtime_edge_2;
    dtime_grid.TimeZone = sig1000.dtime.TimeZone;

    %
    Nlengthtimeseries = length(sig1000.dtime);
    list_time_vars = {'timedatenum', 'dtime'};

    %
    for i2 = 1:length(list_variables_aux)

        % Only interpolate fields that are NOT contained in list_time_vars
        % and are vectors of the correct length
        if ~any(contains(list_time_vars, list_variables_aux{i2})) && ...
           (length(sig1000.(list_variables_aux{i2})) == Nlengthtimeseries)
           
            %
            sig1000.(list_variables_aux{i2}) = ...
                            interp1(sig1000.dtime, ...              
                                    sig1000.(list_variables_aux{i2}), ...
                                    dtime_grid);

            % Turn to column vector
            sig1000.(list_variables_aux{i2}) = sig1000.(list_variables_aux{i2})(:);

        end
    end

    % Replace measured time stamps by time grid
    sig1000.dtime = dtime_grid(:);

    %
    disp('--- Done with time gridding ---')

    % Remove datenum time
    sig1000 = rmfield(sig1000, 'timedatenum');

    % Add sampling rate as a field
    sig1000.samplingrateHz = df_sampling;


    %% Repeat the plot of scalar variables again

    %
    disp('--- Plot of the timeseries in the output structure ---')

    %
    fig_L1_procdata = figure;
    set(fig_L1_procdata, 'units', 'normalized')
    set(fig_L1_procdata, 'Position', [0.2, 0.2, 0.4, 0.6])
        %
        haxs_1 = axes(fig_L1_procdata, 'Position', [0.1, 0.73, 0.8, 0.17]);
        haxs_2 = axes(fig_L1_procdata, 'Position', [0.1, 0.52, 0.8, 0.17]);
        haxs_3 = axes(fig_L1_procdata, 'Position', [0.1, 0.31, 0.8, 0.17]);
        haxs_4 = axes(fig_L1_procdata, 'Position', [0.1, 0.10, 0.8, 0.17]);
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
    ylim(haxs_1, [min(sig1000.pressure), max(sig1000.pressure)])
    %
    ylim(haxs_2, [min(sig1000.temperature), max(sig1000.temperature)])
    ylim(haxs_3, [min(sig1000.heading), max(sig1000.heading)] + [-1, +1])
    ylim(haxs_4, [min([min(sig1000.pitch), min(sig1000.roll)]), ...
                  max([max(sig1000.pitch), max(sig1000.roll)])] + [-1, +1])

    %
    ylabel(haxs_1, '[dbar]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, '[$^\circ$C]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_3, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
                   char(sig1000.SN) ': pressure, temperature, heading, ' ...
                   'pitch (blue), and roll (red)'], ...
                  'Interpreter', 'Latex', 'FontSize', 12)
    %
    linkaxes([haxs_1, haxs_2, haxs_3, haxs_4], 'x')

% %     % Saved at the end of the for loop
% %     exportgraphics(fig_L1_procdata, fullfile(pwd, ['sig_scalars_proc_' list_Signature{i1} '.png']), 'Resolution', 300)


    %% Compute rotate velocity if lvelbin1 is true

    %
    if lvelbin1

        % -------------------------------------------------------
        % Rotate velocity to magnetic ENU

        % If data is too long, break it apart (maybe not necessary when
        % there is only 1 velocity bin)
        npts_rot_TH = 4000000;
        Npts_alldata = length(sig1000.dtime);
        %
        if Npts_alldata<=npts_rot_TH
            indbreak_rot = [1; Npts_alldata];
        else
            %
            inds_edges_breakrot = 1 : npts_rot_TH : Npts_alldata;
            if (inds_edges_breakrot(end)~=Npts_alldata)
                inds_edges_breakrot = [inds_edges_breakrot, Npts_alldata];
            end
            %
            indbreak_rot = [inds_edges_breakrot(1:end-1); ...
                            (inds_edges_breakrot(2:end) - 1)];
            % The last one shouldn't have 1 subtracted. Add again
            indbreak_rot(2, end) = indbreak_rot(2, end) + 1;
        end

        %
        sig1000.Ue = NaN(size(sig1000.vel1));
        sig1000.Vn = NaN(size(sig1000.vel1));
        sig1000.Wup = NaN(size(sig1000.vel1));
        
        %
        disp('--- Using 4 beams ---')

        %
        disp(['--- Coordinate transformation will be computed ' ...
              'for ' num2str(size(indbreak_rot, 2)) ' separate chunks ' ...
              'of the timeseries (to avoid crashing Matlab) ---'])

        %
        for i2 = 1:size(indbreak_rot, 2)

            %
            ind_sub_aux = indbreak_rot(1, i2) : indbreak_rot(2, i2);
            %
            [sig1000.Ue(ind_sub_aux), ...
             sig1000.Vn(ind_sub_aux), ...
             sig1000.Wup(ind_sub_aux)] = ...
                        janus2earth(sig1000.heading(ind_sub_aux).' - 90, ...
                                    sig1000.roll(ind_sub_aux).', -sig1000.pitch(ind_sub_aux).', ...
                                    25, ...
                                    -sig1000.vel1(ind_sub_aux).', -sig1000.vel3(ind_sub_aux).', ...
                                    -sig1000.vel4(ind_sub_aux).', -sig1000.vel2(ind_sub_aux).');

            %
            disp(['----- Done with chunk ' num2str(i2) ' out of ' num2str(size(indbreak_rot, 2)) ' -----'])
        end


        % -------------------------------------------------------
        % Rotate velocity to geographic ENU
        sig1000.magdec = 12.86;

        %
        disp(['----- Rotating horizontal velocity from magnetic ENU ' ...
              'to local XY coordinate system. Magnetic declination is ' ...
              num2str(sig1000.magdec, '%.2f') ' degrees -----'])
    
        %
        [sig1000.Ux, sig1000.Uy] = ROXSI_uv_ENtoXY(sig1000.Ue, sig1000.Vn, sig1000.site, true);

% %         %
% %         rotMatrix = [cosd(sig1000.magdec), sind(sig1000.magdec); ...
% %                      -sind(sig1000.magdec), cosd(sig1000.magdec)];
% %     
% %     % %     % Check the rotation (i.e. velocity aligned with magnetic north
% %     % %     % should have a small zonal component and large meridional
% %     % %     % component in a geographical north coordinate system)
% %     % %     rotMatrix * [0; 1]
% %     
% %         %
% %         uv_aux = [sig1000.Ue.'; sig1000.Vn.'];
% %         %
% %         uv_rot_aux = rotMatrix * uv_aux;
% % 
% %         %
% %         sig1000.Ue = uv_aux(1, :).';
% %         sig1000.Vn = uv_aux(2, :).';
% % 
% %         % -------------------------------------------------------
% %         %  Rotate velocity to XY ROXSI coordinate system(s)
% % 
% %         %
% %         [sig1000.Ux, sig1000.Uy] = ROXSI_uv_ENtoXY(sig1000.Ue, sig1000.Vn, sig1000.site);
        

        % -------------------------------------------------------
        % Remove other variables
        sig1000 = rmfield(sig1000, {'vel1', 'vel2', 'vel3', 'vel4', 'Ue', 'Vn'});

    end


    


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ------------------------- FINAL STEPS --------------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------

    %% Add README

    %
    time_dataproc = datetime('now', 'TimeZone', 'Local');
    time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');
    % Add a general README
    sig1000.README = ['Level 1 scalar data form Signature1000 in ROXSI ' ...
                      '2022. The data is from Signature with serial ' ...
                      'number SN and deployed at mooring site mooringID. ' ...
                      'Data processed by script ' mfilename() '.m ' ...
                      'on ' time_dataproc_char ' (TimeZone ' time_dataproc.TimeZone '). ' ...
                      'This data structure may be useful if there are ' ...
                      'problems with the velocity data but pressure ' ...
                      'is fine. Deployment period is defined in the ' ...
                      'table deploymentInfo_ROXSI2022.mat. Pressure ' ...
                      'is in dbar, where atmospheric pressure has been removed.'];


    %%


    % ----------------------------------------------------
    % Save level 1 data and QC figures

    %
    if lsave_file
        %
        disp('----- Saving level 1 data -----')

        %
        str_filename = ['roxsi_signature_scalars_L1_' char(sig1000.mooringID) '_' char(sig1000.SN)];
        %
        save(fullfile(dir_output_data_L1, [str_filename '.mat']), 'sig1000', '-v7.3')
    end

    %
    if lsave_fig
        %
        disp('----- Saving level 1 QC plot figures -----')

        % Save figures as *.png
        %
        exportgraphics(fig_L1_QC_tilt, fullfile(pwd, ['sig_scalars_' list_Signature{i1} '.png']), 'Resolution', 300)
        %
        exportgraphics(fig_L1_QC_clock, fullfile(pwd, ['sig_clock_' list_Signature{i1} '.png']), 'Resolution', 300)
        %
        exportgraphics(fig_L1_procdata, fullfile(pwd, ['sig_scalars_proc_' list_Signature{i1} '.png']), 'Resolution', 300)

    end


    %%

    % Clear some variables to avoid issues in the next loop iteration
    close(fig_L1_QC_tilt)
    close(fig_L1_QC_clock)
    close(fig_L1_procdata)
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

