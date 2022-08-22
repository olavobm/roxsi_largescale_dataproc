%% Script that process Smart Mooring/Spotter pressure data only
% from RAW to level 1.

clear
close all


%%
% --------------------------------------
% --------- PRELIMINARY STUFF ----------
% --------------------------------------

%%

%
dir_rawdata_parent = fullfile(data_dirpath(), 'RAW', 'Spotters_Smart', 'SDcards');


%%

%
dir_output_data_L1 = fullfile(data_dirpath(), 'Level1_Data', 'Spotters_Smart_Level1');
dir_output_figs_L1 = fullfile(data_dirpath(), 'Level1_Data', 'Spotters_Smart_Level1', 'qc_plots');
%
% dir_output_data_L1 = '/Volumes/OBM-HD/docs/researchPostdoc/datasets/ROXSI/fieldworks/experiment_2022/Aquadopp/';
% dir_output_figs_L1 = fullfile(dir_output_data_L1, 'qc_p lots');

% Logical switches to save or not save data and figures
lsave_file = false;
lsave_fig = false;


%% List of Spotters/Smart moorings that will be processed

% All Spotters/Smart moorings
list_SmartMoorings = {'E01_spot1851', 'E02_spot1859', ...
                     'E05_spot1853', 'E07_spot1855', 'E07_spot1857', ...
                     'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
                     'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

%
Nspotters = length(list_SmartMoorings);


%% Display on the screen:

% Skip a few lines for clarity
disp(' ')
disp(' ')
disp(' ')
% Print message
disp(['Processing ' num2str(Nspotters) ' Smart Mooring(s) from RAW to level 1. ' ...
      'Here are all Aquadopps:'])

list_SmartMoorings


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------------- DO DATA PROCESSING ----------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------


%% Load atmospheric pressure

% % %
% % atm_pressure = load(['/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch' ...
% %                      '/figures_bydate/2022_08_17/obm_edited_noaa_mry_barometric_pressure/' ...
% %                      'atm_pressure.mat']);


%% Datenum limits of the full deployment (from before the first
% went into the water to after the last came out of the water)

%
deployment_timelimits = [datenum(2022, 06, 16, 0, 0, 0), ...
                         datenum(2022, 07, 23, 0, 0, 0)];


%%

%
for i = 1:length(list_SmartMoorings)

    %
    bla = Spotter_read_SMD_allSDcard(fullfile(dir_rawdata_parent, list_SmartMoorings{i}));

    % ------------------------------------------
    % Convert to datenumb in local time (PDT)
    bla.allfiles.dtime = bla.allfiles.dtime - (7/24);

    % ------------------------------------------
    % Plot timeseries of the clock and pressure
    hfig_basic = figure;
        %
        set(hfig_basic, 'Units', 'normalized')
        set(hfig_basic, 'Position', [0.2, 0.1722, 0.3387, 0.3368])
        %
        haxs_1 = axes('Position', [0.1, 0.75, 0.35, 0.2]);
        haxs_2 = axes('Position', [0.1, 0.425, 0.35, 0.2]);
        haxs_3 = axes('Position', [0.1, 0.1, 0.35, 0.2]);
        %
        haxs_4 = axes('Position', [0.525, 0.2, 0.45, 0.5]);
        %
        hold(haxs_1, 'on')
        hold(haxs_2, 'on')
        hold(haxs_3, 'on')
        hold(haxs_4, 'on')

            %
            plot(haxs_1, bla.allfiles.dtime, '.-k')
            plot(haxs_2, 24*3600*diff(bla.allfiles.dtime), '.-k')
            plot(haxs_3, bla.allfiles.pressure, '.-k')
            %
            plot(haxs_4, bla.allfiles.dtime, bla.allfiles.pressure, '.-k')


        %
        set([haxs_1, haxs_2, haxs_3, haxs_4], ...
                            'FontSize', 16, 'Box', 'on', ...
                            'XGrid', 'on', 'YGrid', 'on')
        %
        set(haxs_1, 'YLim', deployment_timelimits)
        set(haxs_2, 'YLim', [-2, 4])
        set(haxs_3, 'YLim', [0, 3e6])
        set([haxs_1, haxs_2, haxs_3], 'XLim', [-0.05*length(bla.allfiles.dtime), ...
                                              1.05.*length(bla.allfiles.dtime)])
        %
        set(haxs_4, 'XLim', deployment_timelimits, 'YLim', [0.7e6, 2.8e6])
        %
% %         set([haxs_1, haxs_2], 'XTickLabel', [])

        %
        xlabel(haxs_3, 'Indices of each Spotter', 'Interpreter', 'Latex', 'FontSize', 16)
        xlabel(haxs_4, 'Time [PDT]', 'Interpreter', 'Latex', 'FontSize', 16)
        %
        ylabel(haxs_1, 'Time [PDT]', 'Interpreter', 'Latex', 'FontSize', 16)
        ylabel(haxs_2, 'diff(time) [s]', 'Interpreter', 'Latex', 'FontSize', 16)
        ylabel(haxs_3, 'pressure', 'Interpreter', 'Latex', 'FontSize', 16)
        ylabel(haxs_4, 'pressure', 'Interpreter', 'Latex', 'FontSize', 16)
        %
        title(haxs_1, 'Timestamps', 'Interpreter', 'Latex', 'FontSize', 16)
        title(haxs_2, 'diff(time)', 'Interpreter', 'Latex', 'FontSize', 16)
        title(haxs_3, 'pressure', 'Interpreter', 'Latex', 'FontSize', 16)
        %
        title(haxs_4, {['ROXSI 2022: ' list_SmartMoorings{i}(1:3) ' - SN ' ...
                       list_SmartMoorings{10}(end-3:end)];'pressure timeseries'}, ...
                       'Interpreter', 'Latex', 'FontSize', 24)
end
