%% Plot initial QC plots for Aquadopp data:
%   1) timeseries of pressure, heading, pitch, and roll
%   2) timeseries zoomed around recovery time to get the
%      recovery time and look at pressure when instrument
%      was OUT of the water.
%
% This is a good script to run to (semi-)mannually get 
% recovery times for each instrument so that these can
% be added to deploymentInfo_ROXSI2022.mat.

clear
close all


%%
% --------------------------------------
% --------- PRELIMINARY STUFF ----------
% --------------------------------------

%%

%
dir_rawdata_parent = fullfile(data_dirpath(), 'RAW', 'Aquadopp');


%%

%
% % dir_output_figs = fullfile(data_dirpath(), 'Level1_Data', 'Aquadopp_Level1', 'qc_plots');
dir_output_figs = pwd;

% Logical switch to save or not save data and figuress
lsave_fig = true;


%% Load ADCP deployment information

% Just to be clear: file and variable have
% the same name (though not a requirement)
load(fullfile(repo_dirpath(), 'deploymentInfo_ROXSI2022.mat'), 'deploymentInfo_ROXSI2022')


%% List of Aquadopps that will be processed

% % % All Aquadopps
% % list_Aquadopp = {'A03_5380', ...
% %                  'B02_12507', 'B04_2147', 'B07_2141', 'B08_13288', 'B11_12280', ...
% %                  'C03_0709', ...
% %                  'D01_12346', 'D02_0653', ...
% %                  'E03_13300', 'E04_13172', 'E06_9736', 'E12_11150', ...
% %                  'F01_9995', 'F02_5838', 'F03_5384', 'F04_5401', 'F05_14032', ...
% %                  'X06_13290', 'X13_9945'};


% A few tests
list_Aquadopp = {'B02_12507'};   % a little data (averaging)
% % list_Aquadopp = {'E12_11150'};
% list_Aquadopp = {'B07_2141'};
% % list_Aquadopp = {'F02_5838'};    % lot of data (1 Hz)

%
Naquadopps = length(list_Aquadopp);


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------------- NOW DO THE QC PLOTS ---------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------


%% Define time limits when ALL instruments were in the water
% to format plots and make them more useful

time_limits_format = [datetime(2022, 06, 25), datetime(2022, 07, 18)];
time_limits_format.TimeZone = 'America/Los_Angeles';


%% Load NOAA's atmospheric pressure

%
atm_pressure = load(['/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch' ...
                     '/figures_bydate/2022_08_17/obm_edited_noaa_mry_barometric_pressure/' ...
                     'atm_pressure.mat']);


%% Variables that will be read from the *.sen file

%
list_senfile_vars = ["time", "heading", "pitch", "roll", "pressure", "temperature"];


%% Start the QC plots

% Loop over Aquadopps
tic
for i = 1:Naquadopps

    % ----------------------------------------------------
    %
    disp(' ')
    disp(' ')
    disp(['----- Start *.sen QC for Aquadopp data: ' list_Aquadopp{i} ' -----'])


    % ----------------------------------------------------
    % Load variables in *.sen file
    file_sen_aux = dir(fullfile(dir_rawdata_parent, list_Aquadopp{i}, '*.sen'));

    %
    senAQDP = Aquadopp_read_senfile(fullfile(file_sen_aux.folder, file_sen_aux.name), list_senfile_vars);


    % ----------------------------------------------------
    % Correct for clock drift
    time_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, list_Aquadopp{i}(5:end), datenum(senAQDP.time));
    %
    senAQDP.time = datetime(time_aux, ...
                                    'ConvertFrom', 'datenum', ...
                                    'TimeZone', 'America/Los_Angeles');

    % Get which row the i'th Aquadopp is in deploymentInfo_ROXSI2022
    ind_row_match = find(strncmp(deploymentInfo_ROXSI2022.mooringID, list_Aquadopp{i}(1:3), 3));

    % ----------------------------------------------------
    % ------- Make simple QC plot with timeseries --------
    % ------ of pressure, heading, pitch, and roll -------
    % ----------------------------------------------------
    
    %
    figure
        %
        hfig_timeseries = gcf;
        set(hfig_timeseries, 'units', 'normalized', 'Position', [0.3, 0.3, 0.4, 0.4])

        %
        haxs_1 = axes('Position', [0.1, 0.7375, 0.8, 0.1625]);
        haxs_2 = axes('Position', [0.1, 0.5250, 0.8, 0.1625]);
        haxs_3 = axes('Position', [0.1, 0.3125, 0.8000 0.1625]);
        haxs_4 = axes('Position', [0.1, 0.1, 0.8, 0.1625]);
        %
        hold(haxs_1, 'on')
        hold(haxs_2, 'on')
        hold(haxs_3, 'on')
        hold(haxs_4, 'on')

        %
        plot(haxs_1, senAQDP.time, senAQDP.pressure, 'k')
        plot(haxs_2, senAQDP.time, senAQDP.heading, 'k')
        plot(haxs_3, senAQDP.time, senAQDP.pitch, 'k')
        plot(haxs_4, senAQDP.time, senAQDP.roll, 'k')

    %
    lin_timelims = (senAQDP.time >= time_limits_format(1)) & ...
                   (senAQDP.time <= time_limits_format(2));
    %
    min_max_pressure = [min(senAQDP.pressure(lin_timelims)), max(senAQDP.pressure(lin_timelims))];
    min_max_heading = [min(senAQDP.heading(lin_timelims)), max(senAQDP.heading(lin_timelims))];
    min_max_pitch = [min(senAQDP.pitch(lin_timelims)), max(senAQDP.pitch(lin_timelims))];
    min_max_roll = [min(senAQDP.roll(lin_timelims)), max(senAQDP.roll(lin_timelims))];
    %
    ylim(haxs_1, min_max_pressure + [-1.5, +0.5]);
    ylim(haxs_2, min_max_heading)
    ylim(haxs_3, min_max_pitch + [-1, 1])
    ylim(haxs_4, min_max_roll + [-1, 1])

    % Time limits starting before ALL instrument deployments and after ALL
    % instrument recoveries
    xaxs_lims = [datetime(2022, 06, 21, 12, 00, 0, 'TimeZone', 'America/Los_Angeles'), ...
                 datetime(2022, 07, 25, 10, 40, 0, 'TimeZone', 'America/Los_Angeles')];
    %
    set([haxs_1, haxs_2, haxs_3, haxs_4], 'XLim', xaxs_lims)
    %
    set([haxs_1, haxs_2, haxs_3, haxs_4], ...
                    'FontSize', 16, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
                    'XTick', xaxs_lims(1) : days(5) : xaxs_lims(2))
    set([haxs_1, haxs_2, haxs_3], 'XTickLabel', [])

    %
    ylabel(haxs_1, '[dbar]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel(haxs_2, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel(haxs_3, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 18)
    %
    title(haxs_1, ['ROXSI 2022: Aquadopp ' list_Aquadopp{i}(1:3) ' - SN ' ...
                   list_Aquadopp{i}(5:end) ': pressure, heading, ' ...
                   'pitch, and roll (clock drift = ' ...
                   num2str(deploymentInfo_ROXSI2022.clockdrift(ind_row_match), '%.2f') ' s)'], ...
                  'Interpreter', 'Latex', 'FontSize', 18)

    %
    % ----------------------------------------------------
    % ------- MAKE PLOT ZOOMING IN AROUND RECOVERY -------
    % ----------------------------------------------------
    

    % ---------------------------------
    % Just select the data later in time
    % to deal with less data
    llater_record = (senAQDP.time >= ...
                         datetime(2022, 07, 19, 'TimeZone', 'America/Los_Angeles'));
    % B10 was recovered later, so use a difference time reference
    if strcmp(list_Aquadopp{i}(1:3), 'B10')
        llater_record = (senAQDP.time >= ...
                          datetime(2022, 07, 24, 'TimeZone', 'America/Los_Angeles'));
    end

    %
    time_aux = senAQDP.time(llater_record);
    pressure_aux = senAQDP.pressure(llater_record);
    pitch_aux = senAQDP.pitch(llater_record);    % just for plotting
    %
    time_mid_aux = datetime((datenum(time_aux(1:end-1)) + ...
                             datenum(time_aux(2:end)))./2, ...
                             'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    diff_pres_aux = diff(pressure_aux);
    

    % ---------------------------------
    % Do bin averaging to get automatically
    % identify a pressure jump that is
    % likely associated with recovery
    %
    timediff_grid = 20;    % in minutes
    %
    time_grid_aux = time_aux(1) : minutes(timediff_grid) : time_aux(end);
    pres_grid_aux = NaN(1, length(time_grid_aux));

    % Do a very inneficient (but general) bin averaging
    tic
    % This can take O(1) min to average
    for i2 = 1:length(time_grid_aux)
        %
        lin_gridlims_aux = (time_aux >= (time_grid_aux(i2) - minutes(1))) & ...
                           (time_aux <= (time_grid_aux(i2) + minutes(1)));
        %
        pres_grid_aux(i2) = mean(pressure_aux(lin_gridlims_aux));
    end
    toc


    % Define zoom limits around recovery
    [~, ind_recovery] = max(abs(diff(pres_grid_aux)));
    %
    time_mid_grid_aux = datetime((datenum(time_grid_aux(1:end-1)) + ...
                                  datenum(time_grid_aux(2:end)))./2, ...
                                'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    %
    time_zoom_lims = [(time_mid_grid_aux(ind_recovery) - minutes(timediff_grid)), ...
                      (time_mid_grid_aux(ind_recovery) + minutes(timediff_grid))];

    %
    linzoomlim_time_aux = (time_aux >= time_zoom_lims(1)) & (time_aux <= time_zoom_lims(2));
    linzoomlim_midtime_aux = (time_mid_aux >= time_zoom_lims(1)) & (time_mid_aux <= time_zoom_lims(2));


    % ---------------------------------
    % Make figure
    figure
        %
        hfig_recovery = gcf;
        set(hfig_recovery, 'units', 'normalized', 'Position', [0.3, 0.3, 0.4, 0.4])

        %
        haxs_fullpres = axes('Position', [0.075, 0.575, 0.35, 0.325]);
        haxs_zoompres = axes('Position', [0.5, 0.7, 0.475, 0.25]);
        haxs_zoomdiffpres = axes('Position', [0.5, 0.4, 0.475, 0.25]);
        %
        haxs_fullpitch = axes('Position', [0.075, 0.1, 0.35, 0.325]);
        haxs_zoompitch = axes('Position', [0.5, 0.1, 0.475, 0.25]);
        %
        hold(haxs_fullpres, 'on')
        hold(haxs_zoompres, 'on')
        hold(haxs_zoomdiffpres, 'on')
        hold(haxs_fullpitch, 'on')
        hold(haxs_zoompitch, 'on')

            % ------------------------------------
            %
            plot(haxs_fullpres, time_aux, pressure_aux, '.-')
            %
            plot(haxs_zoompres, time_aux(linzoomlim_time_aux), pressure_aux(linzoomlim_time_aux), '.-')
            %
            plot(haxs_zoomdiffpres, time_mid_aux(linzoomlim_midtime_aux), diff_pres_aux(linzoomlim_midtime_aux), '.-')

            % Plot atmospheric pressure corrected for what it should be
            % (atmospheric pressure was around
            % 1017 mbar when pressure was offset before deployment)
            plot(haxs_fullpres, atm_pressure.time_vec, (atm_pressure.atm_pres - 1017 + 100)./100)
            plot(haxs_zoompres, atm_pressure.time_vec, (atm_pressure.atm_pres - 1017 + 100)./100)

            % ------------------------------------
            %
            plot(haxs_fullpitch, time_aux, pitch_aux, '.-')
            plot(haxs_zoompitch, time_aux(linzoomlim_time_aux), pitch_aux(linzoomlim_time_aux), '.-')

    % ------------------------------------
    %
    set([haxs_fullpres, haxs_zoompres, haxs_zoomdiffpres, ...
         haxs_fullpitch, haxs_zoompitch], ...
                      'FontSize', 16, 'Box', 'on', ...
                      'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_zoomdiffpres, 'YLim', max(abs(get(haxs_zoomdiffpres, 'YLim'))) .* [-1, 1]);

    % Plot vertical dashed lines for limits
    for i2 = [haxs_fullpres, haxs_zoompres, haxs_zoomdiffpres, haxs_zoompitch]
        %
        y_axs_aux = ylim(i2);
        %
        plot((i2), [time_zoom_lims(1), time_zoom_lims(1)], y_axs_aux, '--k')
        plot((i2), [time_zoom_lims(2), time_zoom_lims(2)], y_axs_aux, '--k')
    end

    %
    edge_extra = seconds(30 * timediff_grid/2);
    %
    set([haxs_fullpres, haxs_fullpitch], 'XLim', time_aux([1, end]))
    set([haxs_zoompres, haxs_zoomdiffpres, haxs_zoompitch], 'XLim', time_zoom_lims + [-edge_extra, edge_extra])

    %
    xlabel(haxs_fullpres, 'Time (drift-corrected)', 'Interpreter', 'latex', 'FontSize', 16)
    xlabel(haxs_fullpitch, 'Time (drift-corrected)', 'Interpreter', 'latex', 'FontSize', 16)
    xlabel(haxs_zoompitch, 'Time (drift-corrected)', 'Interpreter', 'latex', 'FontSize', 16)
    %
    ylabel(haxs_fullpres, 'Pressure [dbar]', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel(haxs_zoompres, 'Pressure [dbar]', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel(haxs_zoomdiffpres, 'diff(pres) [dbar]', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel(haxs_fullpitch, 'pitch [degrees]', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel(haxs_zoompitch, 'pitch [degrees]', 'Interpreter', 'latex', 'FontSize', 16)

    %
    title(haxs_fullpres, {['ROXSI 2022: Aquadopp ' list_Aquadopp{i}(1:3) ' - SN ' list_Aquadopp{i}(5:end) ':']; 'pressure around recovery'}, 'Interpreter', 'latex', 'FontSize', 20)
    title(haxs_fullpitch, 'pitch around recovery', 'Interpreter', 'latex', 'FontSize', 20)
    %
    title(haxs_zoompres, 'Zoom: pressure (and atm pres), diff(pres), pitch', 'Interpreter', 'latex', 'FontSize', 20)
    
    % ----------------------------------------------------
    % Save *.sen QC plots

    % Save figure
    if lsave_fig
        %
        disp('Saving *.sen QC 1 QC plot figure.')

        %
        str_filename = ['roxsi_aquadopp_senQC_' list_Aquadopp{i}(5:end) '_' list_Aquadopp{i}(1:3)];

        % Save figure as *.png
        exportgraphics(hfig_timeseries, fullfile(dir_output_figs, [str_filename '_timeseries.png']), 'Resolution', 300)
        exportgraphics(hfig_recovery, fullfile(dir_output_figs, [str_filename '_recovery.png']), 'Resolution', 300)

    end

    %
    disp(['----- DONE *.sen QC for Aquadopp data proc: ' list_Aquadopp{i} ' -----'])
    %
    toc

end

    
