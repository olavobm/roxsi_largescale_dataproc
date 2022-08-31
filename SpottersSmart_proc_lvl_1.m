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

% % % All Spotters/Smart moorings
% % list_SmartMoorings = {'E01_spot1851', 'E02_spot1859', ...
% %                      'E05_spot1853', 'E07_spot1855', 'E07_spot1857', ...
% %                      'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
% %                      'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

% Look at just one that we didn't have problems with
list_SmartMoorings = {'E02_spot1859'};

%
Nspotters = length(list_SmartMoorings);


%% Display on the screen:

% Skip a few lines for clarity
disp(' ')
disp(' ')
disp(' ')
% Print message
disp(['Processing ' num2str(Nspotters) ' Smart Mooring(s) from RAW to level 1. ' ...
      'Processing Aquadopps:'])

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
    raw_readdata = Spotter_readmulti_SMD(fullfile(dir_rawdata_parent, list_SmartMoorings{i}));

    % ------------------------------------------
    % In datenum, convert from UTC to local time (PDT)
    raw_readdata.allfiles.dtime = raw_readdata.allfiles.dtime - (7/24);


    % ------------------------------------------
    % Variable to only get data that has the correct
    % flag (link) from the Spotter saying that it
    % is indeed data (it seems more appropriate
    % than NaNing data)
    l_gooddata = (raw_readdata.allfiles.link==1);

    %
    spotterSmartdata.dtime = raw_readdata.allfiles.dtime(l_gooddata);
    spotterSmartdata.pressure = raw_readdata.allfiles.pressure(l_gooddata);
    %
    spotterSmartdata.unixEpoch = raw_readdata.allfiles.unixEpoch(l_gooddata);


    % ------------------------------------------
    % Trim between times when instrument was in
    % the water (which also removes timestamps at
    % 31-Dec-1969 17:00:00, thrown by the GPS when????)

    % JUST A PLACEHOLDER!!!!
    time_lims_aux = [datenum(2022, 06, 17, 07, 0, 0), datenum(2022, 07, 19, 0, 0, 0)];

    %
    lintime_lims_aux = (spotterSmartdata.dtime >= time_lims_aux(1)) & ...
                       (spotterSmartdata.dtime <= time_lims_aux(2));

    %
    spotterSmartdata.dtime = spotterSmartdata.dtime(lintime_lims_aux);
    spotterSmartdata.pressure = spotterSmartdata.pressure(lintime_lims_aux);
    spotterSmartdata.unixEpoch = spotterSmartdata.unixEpoch(lintime_lims_aux);


    % ------------------------------------------
    % Fix times when clock goes backwards -- UNLESS THERE IS AN EXTRA
    % PIECE OF INFORMATION ON HOW THIS ERROR HAPPENS, THEN I HAVE TO
    % MAKE AN ASSUMPTION ABOUT THE TIME DIFFERENCE BETWEEN THE NORMAL
    % TIMESERIES AND THE SEGMENT WITH THE PROBLEM

    %
    timediff_aux = 24*3600*diff(spotterSmartdata.dtime);
    % THESE TIME DIFFERENCES ARE NOT AS EXACT AS IN UNIXEPOCH!!!
    % (e.g. a -1 becomes -0.999994575977325).
    % Maybe because of chainging time zone through 7/24????

    %
    inds_gobacks = find(timediff_aux < -0.1);

    % Threshold of the maximum length (in number of points) of
    % the segment to be fixed -- MAYBE NOT THE BEST APPROACH, BUT
    % MAYBE GOOD ENOUGH
    NsegTH = 20;

    %
    for i2 = 1:length(inds_gobacks)

        %
        [~, ind_back_normaltime_aux] = max(timediff_aux(inds_gobacks(i2):(inds_gobacks(i2) + NsegTH)));

        % Indices of the segment when the timestamps have gone back in time
        ind_seg_tofix = inds_gobacks(i2) + (1 : (ind_back_normaltime_aux - 1));

        % Make a simple correction in terms of multiples
        % of the sampling period (0.5 s)
        %
        integer_division_aux = floor(10*abs(timediff_aux(inds_gobacks(i2))))/5;
        %
        if integer_division_aux < 1
            integer_division_aux = 0;
        end
        %
        time_factor_fix_aux = 0.5*(integer_division_aux + 1);    % in seconds


        % Add the time factor to correct the timestamps of the segment
        spotterSmartdata.dtime(ind_seg_tofix) = spotterSmartdata.dtime(ind_seg_tofix) + (time_factor_fix_aux ./ (24*3600));
        
       
        % ------------------------
        % Make diagnostic plot -- there are hundreds of
        % instances per spotter, so it's unlikely you
        % want to save all of them

        %
        lmakeplot = false;
        if integer_division_aux < 1
            lmakeplot = true;
        end

        %
        if lmakeplot

            % Indices to plot
            indsplt = (ind_seg_tofix(1)-10):(ind_seg_tofix(end)+10);
    % %         indsplt = (ind_seg_tofix(1)-200):(ind_seg_tofix(end)+400);
    
            %
            hfig = figure;
                set(hfig, 'units', 'normalized')
                set(hfig, 'Position', [0.39, 0.49, 0.26, 0.43])

            %
            newFigDims([9.25, 8.6389])
                %
                haxs_1 = axes('Position', [0.15, 0.7, 0.75, 0.18]);
                haxs_2 = axes('Position', [0.15, 0.4, 0.75, 0.18]);
                haxs_3 = axes('Position', [0.15, 0.1, 0.75, 0.18]);
                %
                hold(haxs_1, 'on')
                hold(haxs_2, 'on')
                hold(haxs_3, 'on')
    
                    % Uncorrected
                    plot(haxs_1, indsplt, datetime((719529 + (spotterSmartdata.unixEpoch(indsplt)./86400) - (7/24)), 'ConvertFrom', 'datenum'), '.-', 'MarkerSize', 20)
                    % Corrected
                    plot(haxs_1, indsplt, datetime(spotterSmartdata.dtime(indsplt), 'ConvertFrom', 'datenum'), '.-', 'MarkerSize', 20)
    
                    % Plot the time difference in subplots 2 and 3
                    for indhaxs = [haxs_2, haxs_3]
                        plot(indhaxs, (indsplt(1:end-1)+indsplt(2:end))./2, 24*3600*diff((719529 + (spotterSmartdata.unixEpoch(indsplt)./86400) - (7/24))), '.-', 'MarkerSize', 20)
                        plot(indhaxs, (indsplt(1:end-1)+indsplt(2:end))./2, 24*3600*diff(spotterSmartdata.dtime(indsplt)), '.-', 'MarkerSize', 20)
                    end    

            %
            hleg = legend(haxs_1, 'uncorrected', 'corrected', 'Location', 'NorthWest', 'FontSize', 14);
    
    
            %
            set([haxs_1, haxs_2, haxs_3], 'FontSize', 16, 'Box', 'on', ...
                                          'XGrid', 'on', 'YGrid', 'on')
            set([haxs_1, haxs_2, haxs_3], 'XLim', indsplt([1, end]) + [-2, 2])
            set(haxs_3, 'YLim', [0.42, 0.58])
    
            %
            plot(haxs_2, get(haxs_2, 'XLim'), [0, 0], 'k')
            plot(haxs_3, get(haxs_2, 'XLim'), [0, 0], 'k')
    
            %
            xlabel(haxs_3, 'Indices of data points', 'Interpreter', 'Latex', 'FontSize', 22)
            %
            ylabel(haxs_1, 'Local time', 'Interpreter', 'Latex', 'FontSize', 22)
            ylabel(haxs_2, '[s]', 'Interpreter', 'Latex', 'FontSize', 22)
            ylabel(haxs_3, '[s]', 'Interpreter', 'Latex', 'FontSize', 22)
            
            %
            title(haxs_1, {['ROXSI 2022: ' list_SmartMoorings{i}(1:3) ' - SN ' ...
                           list_SmartMoorings{i}(end-3:end) ':'];'sequential timestamps'}, ...
                           'Interpreter', 'Latex', 'FontSize', 18)
            title(haxs_2, 'Time difference between timestamps', 'Interpreter', 'Latex', 'FontSize', 18)
            title(haxs_3, 'Same as above, but zoomed in the $y$ axis', 'Interpreter', 'Latex', 'FontSize', 18)
    
            %
            linkaxes([haxs_1, haxs_2, haxs_3], 'x')
    
            %
            keyboard
        end

    end




    keyboard
    %
    inds_break_periods = find(inds_gobacks > 6000);   % ????????
    % Remove those earlier on
    inds_break_periods = inds_break_periods(inds_break_periods > 40);

    % Select about 10 of them (about because depends on index stepping):
    inds_break_selec = 1 : round(length(inds_break_periods)/10) : length(inds_break_periods);

    % Now select corresponding indices in the data
    inds_plt_segments = inds_gobacks(inds_break_periods(inds_break_selec));




    
    %
%     diff(spotterSmartdata.dtime) < 0

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
            plot(haxs_1, raw_readdata.allfiles.dtime, '.-k')
            plot(haxs_2, 24*3600*diff(raw_readdata.allfiles.dtime), '.-k')
            plot(haxs_3, raw_readdata.allfiles.pressure, '.-k')
            %
            plot(haxs_4, raw_readdata.allfiles.dtime, raw_readdata.allfiles.pressure, '.-k')


        %
        set([haxs_1, haxs_2, haxs_3, haxs_4], ...
                            'FontSize', 16, 'Box', 'on', ...
                            'XGrid', 'on', 'YGrid', 'on')
        %
        set(haxs_1, 'YLim', deployment_timelimits)
        set(haxs_2, 'YLim', [-2, 4])
        set(haxs_3, 'YLim', [0, 3e6])
        set([haxs_1, haxs_2, haxs_3], 'XLim', [-0.05*length(raw_readdata.allfiles.dtime), ...
                                              1.05.*length(raw_readdata.allfiles.dtime)])
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
                       list_SmartMoorings{i}(end-3:end)];'pressure timeseries'}, ...
                       'Interpreter', 'Latex', 'FontSize', 24)

        %
        linkaxes([haxs_1, haxs_2, haxs_3], 'x')

    % ------------------------------------------
    % Similar as above, but plot loaded variables and processed data
    % on top of each other.



return
    % ------------------------------------------
    % Plot gap associated with chaging Spotter mode on the dashboard
    
    % For Spotter E02, the gap is between 17-Jun-2022 16:48:53 and 17-Jun-2022 17:11:51
    %
    time_lims_gap = [datenum(2022, 06, 17, 17, 1, 0), datenum(2022, 06, 17, 17, 04, 30)];
    %
    ind_first_aux = find(raw_readdata.allfiles.dtime > time_lims_gap(1), 1, 'first');
    ind_last_aux = find(raw_readdata.allfiles.dtime < time_lims_gap(2), 1, 'last');
    %
    hfig_gap = figure;
        %
        set(hfig_gap, 'Units', 'normalized')
        set(hfig_gap, 'Position', [0.2, 0.1722, 0.3387, 0.3368])
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
            plot(haxs_1, raw_readdata.allfiles.dtime(ind_first_aux:ind_last_aux), '.-k')
            plot(haxs_2, 24*3600*diff(raw_readdata.allfiles.dtime(ind_first_aux:ind_last_aux)), '.-k')
            plot(haxs_3, raw_readdata.allfiles.pressure(ind_first_aux:ind_last_aux), '.-k')
            %
            plot(haxs_4, raw_readdata.allfiles.dtime(ind_first_aux:ind_last_aux), ...
                         raw_readdata.allfiles.pressure(ind_first_aux:ind_last_aux), '.-k')


        %
        set([haxs_1, haxs_2, haxs_3, haxs_4], ...
                            'FontSize', 16, 'Box', 'on', ...
                            'XGrid', 'on', 'YGrid', 'on')

        %
        pres_aux = raw_readdata.allfiles.pressure(ind_first_aux:ind_last_aux);

% %         %
        set(haxs_1, 'YLim', time_lims_gap)
        set(haxs_2, 'YLim', [-2, 4])
        set(haxs_3, 'YLim', [min(pres_aux), max(pres_aux)])
        set([haxs_1, haxs_2, haxs_3], 'XLim', [0, length(ind_first_aux:ind_last_aux)])
        %
        set(haxs_4, 'XLim', time_lims_gap, 'YLim', [min(pres_aux), max(pres_aux)])

        %
        xticks_aux = linspace(raw_readdata.allfiles.dtime(ind_first_aux), ...
                              raw_readdata.allfiles.dtime(ind_last_aux), ...
                              7);
        set(haxs_4, 'XTick', xticks_aux)
        %
        cell_xticklabel = cell(1, 7);
        %
        for i3 = 1:7
            cell_xticklabel{i3} = num2str(24*3600*(xticks_aux(i3) - xticks_aux(1)), '%.1f');
        end
        set(haxs_4, 'XTickLabel', cell_xticklabel)


        %
        xlabel(haxs_3, 'Indices of each Spotter', 'Interpreter', 'Latex', 'FontSize', 16)
        xlabel(haxs_4, ['Time, seconds since ' datestr(xticks_aux(1)) ' [PDT]'], 'Interpreter', 'Latex', 'FontSize', 16)
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
                       list_SmartMoorings{i}(end-3:end)];'pressure timeseries'}, ...
                       'Interpreter', 'Latex', 'FontSize', 24)

        %
        linkaxes([haxs_1, haxs_2, haxs_3], 'x')

    %
    exportgraphics(gcf, 'smartmooring_E02_gapmode.png', 'Resolution', 300);

    %
    return

    % ------------------------------------------
    % Plot unique values of clock finite difference
% % %     hclockfig = figure;
% % %         %
% % %         time_diff_unique = unique(24*3600*diff(bla.allfiles.dtime));
% % %         %
% % %         lskipstodealwith = abs(time_diff_unique) < 100; 
% % %         
% % %         %
% % %         plot(time_diff_unique(lskipstodealwith), '.-')
% % % 
% % %         %
% % %         grid on
% % %         set(gca, 'FontSize', 16)
% % % 
% % %     %
% % %     time_diff_ok = time_diff_unique(lskipstodealwith);
% % %     %
% % %     time_diff_ok_wrapped = time_diff_ok;
% % %     time_diff_ok_wrapped(time_diff_ok<0) = abs(time_diff_ok_wrapped(time_diff_ok<0));
% % %     time_diff_ok_wrapped(time_diff_ok<0) = time_diff_ok_wrapped(time_diff_ok<0) + 0;
% % % 
% % %     %
% % %     figure
% % %         %
% % %         plot(time_diff_ok_wrapped, '.-')
% % %         %
% % %         grid on
% % %         set(gca, 'FontSize', 16)


    % ------------------------------------------
% %     % Plot diagnostic plot for several examples when clock stopped
% %     %
    timediff_vec = 24*3600*diff(raw_readdata.allfiles.dtime);
% %     %
% %     inds_stopped = find((timediff_vec > -0.1) & (timediff_vec < 0.1));
% %     %
% %     inds_break_periods = find(inds_stopped > 6000);    % what's this 6000?????
% %     % Remove those earlier on
% %     inds_break_periods = inds_break_periods(inds_break_periods > 1500);
% % 
% %     % Select about 10 of them (about because depends on index stepping):
% %     inds_break_selec = 1 : round(length(inds_break_periods)/10) : length(inds_break_periods);
% % 
% %     % Now select corresponding indices in the data
% %     inds_plt_segments = inds_stopped(inds_break_periods(inds_break_selec));

    % ---------------------
    % For clock going backwards
    %
    inds_stopped = find(timediff_vec < -0.1);
    %
    inds_break_periods = find(inds_stopped > 6000);
    % Remove those earlier on
    inds_break_periods = inds_break_periods(inds_break_periods > 40);

    % Select about 10 of them (about because depends on index stepping):
    inds_break_selec = 1 : round(length(inds_break_periods)/10) : length(inds_break_periods);

    % Now select corresponding indices in the data
%     inds_plt_segments = inds_stopped(inds_break_periods(inds_break_selec));
    inds_plt_segments = inds_stopped(inds_break_periods(1));



    %
    for i2 = 1:length(inds_plt_segments)
        %
        inds_segment_aux = (inds_plt_segments(i2) - 70):1:(inds_plt_segments(i2) + 70);
        %
% %         inds_segment_aux = (inds_plt_segments(i2) - 30):1:(inds_plt_segments(i2) + 60);
    
        %
        hfig_diag_aux = figure;
            %
            set(hfig_diag_aux, 'Units', 'normalized')
            set(hfig_diag_aux, 'Position', [0.37, 0.17, 0.47, 0.38])
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
                plot(haxs_1, raw_readdata.allfiles.dtime(inds_segment_aux), '.-k')
                plot(haxs_2, 24*3600*diff(raw_readdata.allfiles.dtime(inds_segment_aux)), '.-k')
                plot(haxs_3, raw_readdata.allfiles.pressure(inds_segment_aux), '.-k')
                %
                plot(haxs_4, raw_readdata.allfiles.dtime(inds_segment_aux), raw_readdata.allfiles.pressure(inds_segment_aux), '.-k')
    
    
    
            %
            set([haxs_1, haxs_2, haxs_3, haxs_4], ...
                                'FontSize', 16, 'Box', 'on', ...
                                'XGrid', 'on', 'YGrid', 'on')

            %
            minmax_p_aux = [min(raw_readdata.allfiles.pressure(inds_segment_aux)), ...
                            max(raw_readdata.allfiles.pressure(inds_segment_aux))];

            %
            set(haxs_1, 'YLim', raw_readdata.allfiles.dtime(inds_segment_aux([1, end])))
%             set(haxs_2, 'YLim', [-0.2, 0.7])
            %
            ytimediff_aux = 24*3600*diff(raw_readdata.allfiles.dtime(inds_segment_aux));
            %
            set(haxs_2, 'YLim', [min(ytimediff_aux), max(ytimediff_aux)])
            set(haxs_3, 'YLim', minmax_p_aux)
            set([haxs_1, haxs_2, haxs_3], 'XLim', [0, length(inds_segment_aux)])
            %
            set(haxs_4, 'XLim', raw_readdata.allfiles.dtime(inds_segment_aux([1, end])), 'YLim', minmax_p_aux)
            %
            xticks_aux = linspace(raw_readdata.allfiles.dtime(inds_segment_aux(1)), ...
                                  raw_readdata.allfiles.dtime(inds_segment_aux(end)), ...
                                  7);
            set(haxs_4, 'XTick', xticks_aux)
            %
            cell_xticklabel = cell(1, 7);
            %
            for i3 = 1:7
                cell_xticklabel{i3} = num2str(24*3600*(xticks_aux(i3) - xticks_aux(1)), '%.1f');
            end
            set(haxs_4, 'XTickLabel', cell_xticklabel)
    
            %
            xlabel(haxs_3, 'Indices along this data segment', 'Interpreter', 'Latex', 'FontSize', 16)
            xlabel(haxs_4, ['Time, seconds since ' datestr(xticks_aux(1)) ' [PDT]'], 'Interpreter', 'Latex', 'FontSize', 16)
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
                           list_SmartMoorings{i}(end-3:end)];'pressure timeseries'}, ...
                           'Interpreter', 'Latex', 'FontSize', 24)
    
            %
            linkaxes([haxs_1, haxs_2, haxs_3], 'x')

            %
% %             exportgraphics(gcf, ['smartmooring_clockQC_timestops_' num2str(i2, '%.02d') '.png'], 'Resolution', 300);
% % % %             exportgraphics(gcf, ['smartmooring_clockQC_timeinversion_' num2str(i2, '%.02d') '.png'], 'Resolution', 300);
    end
    
        %
        return
end
