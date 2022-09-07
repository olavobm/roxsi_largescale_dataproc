%% Script that process Smart Mooring/Spotter pressure data only
% from RAW to level 1.
%
% This script can make ??? figures per Spotter that:
%       -
%       -
%       -
%       -
%       -
%
% Some of the figures are already commented out. Feel
% free to comment/comment out the figures you do or
% don't want to see.
% 
%
%


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

% % Look at just one that we didn't have problems with
% list_SmartMoorings = {'E02_spot1859'};

% All good Smart Moorings
% list_SmartMoorings = {'E01_spot1851', 'E02_spot1859', 'E08_spot1852', 'E10_spot1848'};

% All seemingly great, apart from not extending the whole deployment
% list_SmartMoorings = {'E05_spot1853', 'E07_spot1857', 'E09_spot1856'};
% list_SmartMoorings = {'E07_spot1857', 'E09_spot1856'};

%
list_SmartMoorings = {'E08_spot1852'};


% %
% list_SmartMoorings = {'E05_spot1853'};    % ONE ANNOYING PROBLEM

% % list_SmartMoorings = {'E09_spot1850'};    % one annoying problem

%
list_SmartMoorings = {'E11_spot1860'};

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


%% Load mooring locations
%
% Here we load locations of all moorings. The deployed locations of the
% Spotter moorings are calculated in Spotters_get_all_locations.m. For
% the Spotters, these locations are calculated by fitting a watch
% circle to the locations reported by the Spotter. For the Smart Moorings,
% this location should correspond to the location of bottom mount/pressure
% sensor.

%
mooringtable = load(fullfile(repo_dirpath(), 'ROXSI2022_mooringtable.mat'));
mooringtable = mooringtable.mooringtable;


%% Load deploymente info to trim the data during deployment

%
dplySpotters = load(fullfile(repo_dirpath(), 'deploymentInfo_Spotters_ROXSI2022'));
dplySpotters = dplySpotters.dployInfo_Spotters;


%% Load atmospheric pressure

%
% % atmpres_NOAA = load(fullfile(data_dirpath(), 'noaa_mry_barometric_pressure', 'atm_pressure.mat'));
atmpres_NOAA = load('atm_pressure.mat');


%% Datenum limits of the full deployment (from before the first
% went into the water to after the last came out of the water)

%
deployment_timelimits = [datenum(2022, 06, 16, 0, 0, 0), ...
                         datenum(2022, 07, 23, 0, 0, 0)];


%% Do data processing: first with a section of QC
% plots from the data that is read, and then a
% second section where the data is processed to level 1

% Loop over Smart Moorings
for i1 = 1:length(list_SmartMoorings)

    % ------------------------------------------
    % Load raw data
    disp(' ')
    disp(' ')
    %
    disp(['-------------- Starting level 1 data processing of smart mooring ' list_SmartMoorings{i1}(1:3) ' - SN ' list_SmartMoorings{i1}(end-3:end) ' --------------'])
    %
    disp('---- Reading data ---- ')

    %
    raw_readdata = Spotter_readmulti_SMD(fullfile(dir_rawdata_parent, list_SmartMoorings{i1}));

    %
    disp('---- Done with reading data ---- ')

    % ------------------------------------------
    % In datenum, convert from UTC to local time (PDT)
    raw_readdata.allfiles.dtime = raw_readdata.allfiles.dtime - (7/24);

    % ------------------------------------------
    % Plot full timeseries QC
    SpotterSmart_QC_plot_pressure(raw_readdata.allfiles.dtime, ...
                                  raw_readdata.allfiles.pressure, ...
                                  deployment_timelimits, ...
                                  list_SmartMoorings{i1}(1:3), list_SmartMoorings{i1}(end-3:end));

    % ------------------------------------------
    % Plot QC of gap due to changing Spotter mode
    
    % For Spotter E02, the gap is between 17-Jun-2022 16:48:53 and 17-Jun-2022 17:11:51
    if strcmp(list_SmartMoorings{i1}(1:3), 'E02')
        %
        time_lims_gap = [datenum(2022, 06, 17, 17, 1, 0), datenum(2022, 06, 17, 17, 04, 30)];
        %
        SpotterSmart_QC_plot_pressure(raw_readdata.allfiles.dtime, ...
                                      raw_readdata.allfiles.pressure, ...
                                      time_lims_gap, ...
                                      list_SmartMoorings{i1}(1:3), list_SmartMoorings{i1}(end-3:end));
    end

    % ------------------------------------------
    % Plot QC of clock stopping -- this is associated
    % with times when there is actually no pressure data
    % and the flag (link) is 0. So this "clock stopping"
    % thing is not really part of the data processing.
    %
    timediff_vec = 24*3600*diff(raw_readdata.allfiles.dtime);
    %
    inds_stopped = find((timediff_vec > -0.1) & (timediff_vec < 0.1));
    %
    inds_break_periods = find(inds_stopped > 6000);    % what's this 6000?????
    % Remove those earlier on
    inds_break_periods = inds_break_periods(inds_break_periods > 1500);

    % Select about 10 of them (about because depends on index stepping):
    inds_break_selec = 1 : round(length(inds_break_periods)/10) : length(inds_break_periods);

    % Now select corresponding indices in the data
    inds_plt_segments = inds_stopped(inds_break_periods(inds_break_selec));
%     %
%     for i2 = 1:length(inds_plt_segments)
%         %
%         inds_segment_aux = (inds_plt_segments(i2) - 70):1:(inds_plt_segments(i2) + 70);
%         %
%         SpotterSmart_QC_plot_pressure(raw_readdata.allfiles.dtime(inds_segment_aux), ...
%                                       raw_readdata.allfiles.pressure(inds_segment_aux), ...
%                                       raw_readdata.allfiles.dtime(inds_segment_aux([1, end])), ...
%                                       list_SmartMoorings{i}(1:3), list_SmartMoorings{i}(end-3:end));
%     end


    % ------------------------------------------
    % Plot QC of clock going back in time
    %
    timediff_vec = 24*3600*diff(raw_readdata.allfiles.dtime);
    %
    inds_timereversal = find(timediff_vec < 0);
    % 
    inds_break_periods = find(inds_timereversal > 6000);    % what's this 6000?????
    % Remove those earlier on
    inds_break_periods = inds_break_periods(inds_break_periods > 40);

    % Select about 10 of them (about because depends on index stepping):
    inds_break_selec = 1 : round(length(inds_break_periods)/10) : length(inds_break_periods);

    % Now select corresponding indices in the data
    inds_plt_segments = inds_timereversal(inds_break_periods(inds_break_selec));

%     %
%     for i2 = 1:length(inds_plt_segments)
%         %
%         inds_segment_aux = (inds_plt_segments(i2) - 70):1:(inds_plt_segments(i2) + 70);
%         %
%         SpotterSmart_QC_plot_pressure(raw_readdata.allfiles.dtime(inds_segment_aux), ...
%                                       raw_readdata.allfiles.pressure(inds_segment_aux), ...
%                                       raw_readdata.allfiles.dtime(inds_segment_aux([1, end])), ...
%                                       list_SmartMoorings{i}(1:3), list_SmartMoorings{i}(end-3:end));
%     end

    % ------------------------------------------
    % Plot the unique values of clock finite difference
    hclockfig = figure;
        %
        time_diff_unique = unique(24*3600*diff(raw_readdata.allfiles.dtime));
        %
        lskipstodealwith = abs(time_diff_unique) < 100; 
        
        %
        plot(time_diff_unique(lskipstodealwith), '.-')

        %
        grid on
        set(gca, 'FontSize', 16)

        %
        xlabel('Indice of unique time difference', 'Interpreter', 'Latex', 'FontSize', 16)
        ylabel('Time difference [s]', 'Interpreter', 'Latex', 'FontSize', 16)
        %
        title(['ROXSI 2022: ' list_SmartMoorings{i1}(1:3) ' SN ' ...
               list_SmartMoorings{i1}(end-3:end) ': unique ' ...
               '(short) time differences'], 'Interpreter', 'Latex', 'FontSize', 16)

    
    %%
    % ----------------------------------------------------------
    % --------------- DO LEVEL 1 DATA PROCESSING ---------------
    % ----------------------------------------------------------

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
    % *****
    % AT LEAST IN 2022, THE SMART MOORINGS
    % HAD A BUG/FEATURE WHERE ALL PRESSURE VALUES
    % ARE RECORDED WITH AN EXTRA ZERO AS THE LAST
    % DIGIT. THAT IS, THE PRESSURE IN PASCAL IS
    % THE VALUE STORED IN THE SD CARD DIVIDED BY 10.
    % *****

    %
    spotterSmartdata.pressure = spotterSmartdata.pressure./10;


    % ------------------------------------------
    % Trim between times when instrument was in
    % the water (which also removes timestamps at
    % 31-Dec-1969 17:00:00, thrown by the GPS when????)

    % Find the Smart Mooring in the deployment table
    lmatch = strcmp(dplySpotters.SN, list_SmartMoorings{i1}(end-3:end));

    %
    time_lims_aux = [datenum(dplySpotters(lmatch, :).time_begin_trim, "yyyy/mm/dd HH:MM:SS"), ...
                     datenum(dplySpotters(lmatch, :).time_end_trim, "yyyy/mm/dd HH:MM:SS")];

% %     % JUST A PLACEHOLDER!!!!
% %     time_lims_aux = [datenum(2022, 06, 17, 07, 0, 0), datenum(2022, 07, 19, 0, 0, 0)];

    %
    lintime_lims_aux = (spotterSmartdata.dtime >= time_lims_aux(1)) & ...
                       (spotterSmartdata.dtime <= time_lims_aux(2));

    %
    spotterSmartdata.dtime = spotterSmartdata.dtime(lintime_lims_aux);
    spotterSmartdata.pressure = spotterSmartdata.pressure(lintime_lims_aux);
    spotterSmartdata.unixEpoch = spotterSmartdata.unixEpoch(lintime_lims_aux);    % this unixEpoch is only used below to have a copy of
                                                                                  % uncorrected time, so that it's compared to corrected time.

    %%
    % -------------------------------------------------
    % -------------------------------------------------
    % -------------------------------------------------
    % Fix times when clock goes backwards -- UNLESS THERE IS AN EXTRA
    % PIECE OF INFORMATION ON HOW THIS ERROR HAPPENS, THEN I HAVE TO
    % MAKE AN ASSUMPTION ABOUT THE TIME DIFFERENCE BETWEEN THE NORMAL
    % PART OF THE TIMESERIES AND THE SEGMENT WITH THE PROBLEM
    %
    % Actually, the same thing also happens with the clock slowing
    % down (i.e. the diff(time) is positive, but closer to 0 than
    % 0.5, and plots of examples convincinly show the problem is
    % in the clock). But in this case, I really need to check if
    % this conclusion from a few examples is valid for all cases
    % of diff(time) below a positive threshold.

    %
    timediff_aux = 24*3600*diff(spotterSmartdata.dtime);
    % THESE TIME DIFFERENCES ARE NOT AS EXACT AS IN UNIXEPOCH.
    % (e.g. a -1 becomes -0.999994575977325).
    % Maybe because of chainging time zone through 7/24????

    %
%     inds_gobacks = find(timediff_aux < 0);
    inds_gobacks = find(timediff_aux < 0.2);

    % Threshold of the maximum length (in number of points) of
    % the segment to be fixed -- MAYBE NOT THE BEST APPROACH, BUT
    % MAYBE GOOD ENOUGH
    NsegTH = 20;

    %
    for i2 = 1:length(inds_gobacks)

        %
        lfix_segment_aux = true;

        % Finds what is likely the corresponding time difference when
        % the clock resets (HOWEVER, there are some underlying assumptions
        % here).
        [~, ind_back_normaltime_aux] = max(timediff_aux(inds_gobacks(i2):(inds_gobacks(i2) + NsegTH)));

        % Indices of the segment when the timestamps have gone back in time
        ind_seg_tofix = inds_gobacks(i2) + (1 : (ind_back_normaltime_aux - 1));

% %         if (spotterSmartdata.dtime(ind_seg_tofix(1))>datenum(2022, 07, 01, 18, 12, 00)) && ...
% %            (spotterSmartdata.dtime(ind_seg_tofix(1))<datenum(2022, 07, 01, 18, 13, 30))
% %             keyboard
% %         end

        % Check the time differences before and after the segment
        % are DEFINITELY CONSISTENT with an error in the clock.
        % (one counter example, which I don't if ever happens would
        % be a real clock delay in the middle of the segment)
        time_diff_before_aux = spotterSmartdata.dtime(ind_seg_tofix(1)) - spotterSmartdata.dtime((ind_seg_tofix(1)-1));
        time_diff_after_aux = spotterSmartdata.dtime((ind_seg_tofix(end)+1)) - spotterSmartdata.dtime(ind_seg_tofix(end));
        %
        time_diff_before_aux = time_diff_before_aux*24*3600;
        time_diff_after_aux = time_diff_after_aux*24*3600;
        
        %
        time_diff_inseg_aux = 24*3600*diff(spotterSmartdata.dtime(ind_seg_tofix));

% %         if (spotterSmartdata.dtime(ind_seg_tofix(1))>datenum(2022, 07, 01, 18, 12, 00)) && ...
% %            (spotterSmartdata.dtime(ind_seg_tofix(1))<datenum(2022, 07, 01, 18, 13, 30))
% %             keyboard
% %         end

        %
        if time_diff_before_aux > 0
            %
            if abs(time_diff_after_aux - 0.5) < 0.8*abs(time_diff_before_aux - 0.5)
% %                 %
% %                 warning('error 1!!!!!!!!')
% %                 
% %                 % From a couple of Spotters error 2 always happens
% %                 % when error 1 happens
% %                 if any(abs(time_diff_inseg_aux - 0.5) > abs(time_diff_after_aux - 0.5))
% %                     warning('error 2!!!!!!!!')
% % 
% %                 %
% %                 figure
% %                     plot((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+54), ...
% %                          24*3600*diff(spotterSmartdata.dtime((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+55))), '.-')
% %                     hold on
% %                     plot(ind_seg_tofix(1:end-1), ...
% %                          24*3600*diff(spotterSmartdata.dtime(ind_seg_tofix)), '.-')
% %                     grid on
% %                     %
% %                     set(gca, 'FontSize', 16)
% %                 
% %                 end

                % Remove the segment in this iteration from
                % those that need to be corrected
                lfix_segment_aux = false;

            end

        %
        else 
            
        % Throw a warning if there is a clock problem within the segment
        % that it's getting fixed -- nice to check, but so far I haven't
        % seen an example of where this goes wrong

% %             %
% %             if abs(time_diff_inseg_aux - 0.5) >= 0.1
% %                 warning('error 3!!!!!!!!')
% %             end
% % 
% %             %
% %             vec_diffaux = timediff_aux(inds_gobacks(i2):(inds_gobacks(i2) + NsegTH));
% %             %
% %             max_timediff = max(vec_diffaux);
% %             %
% %             new_vec_diffaux = setdiff(vec_diffaux, max_timediff, 'stable');
% %             
% %             %
% %             second_max_timediff = max(new_vec_diffaux);
% % 
% %             %
% %             if abs(second_max_timediff - 0.5) >= 0.1
% %                 warning('error 4!!!!!!!!')
% %                 %
% %                 figure
% %                     plot((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+24), ...
% %                          24*3600*diff(spotterSmartdata.dtime((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+25))), '.-')
% %                     hold on
% %                     plot(ind_seg_tofix(1:end-1), ...
% %                          24*3600*diff(spotterSmartdata.dtime(ind_seg_tofix)), '.-')
% %                     grid on
% %                     %
% %                     set(gca, 'FontSize', 16)
% % 
% %             end
        end

        % Only apply correction if it's not identified as
        % an exception by the if block just above
        if lfix_segment_aux

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
        
% %             if (spotterSmartdata.dtime(ind_seg_tofix(1))>datenum(2022, 07, 01, 18, 12, 00)) && ...
% %                (spotterSmartdata.dtime(ind_seg_tofix(1))<datenum(2022, 07, 01, 18, 13, 30))
% %                 keyboard
% %             end
       
            % ------------------------
            % Make diagnostic plot -- there are hundreds of
            % instances per spotter, so it's unlikely you
            % want to save all of them
    
            %
            lmakeplot = false;
    % %         if integer_division_aux < 1
    % %             lmakeplot = true;
    % %         end
    
            %
            if lmakeplot
    
                % Indices to plot -- you may want to plot a shorter
                % ir longer range around the segment that it's being
                % fixed
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
                title(haxs_1, {['ROXSI 2022: ' list_SmartMoorings{i1}(1:3) ' - SN ' ...
                               list_SmartMoorings{i1}(end-3:end) ':'];'sequential timestamps'}, ...
                               'Interpreter', 'Latex', 'FontSize', 18)
                title(haxs_2, 'Time difference between timestamps', 'Interpreter', 'Latex', 'FontSize', 18)
                title(haxs_3, 'Same as above, but zoomed in the $y$ axis', 'Interpreter', 'Latex', 'FontSize', 18)
        
                %
                linkaxes([haxs_1, haxs_2, haxs_3], 'x')
    
            end


        end
    end
    % -------------------------------------------------
    % -------------------------------------------------
    % -------------------------------------------------

    %%
    keyboard
    % ------------------------------------------
    % Plot the unique values of clock finite difference
    % (to check my processing has removed all
    % instances of clock going backwards)
    hclockfig = figure;
        %
        time_diff_unique = unique(24*3600*diff(spotterSmartdata.dtime));
        %
        lskipstodealwith = abs(time_diff_unique) < 100; 
        
        %
        plot(time_diff_unique(lskipstodealwith), '.-')

        %
        grid on
        set(gca, 'FontSize', 16)

        %
        xlabel('Indice of unique time difference', 'Interpreter', 'Latex', 'FontSize', 16)
        ylabel('Time difference [s]', 'Interpreter', 'Latex', 'FontSize', 16)
        %
        title({['ROXSI 2022: ' list_SmartMoorings{i1}(1:3) ' SN ' ...
               list_SmartMoorings{i1}(end-3:end) ': unique ' ...
               '(short) time'];'differences in processed data'}, 'Interpreter', 'Latex', 'FontSize', 20)

    %% Convert pressure from Pascal to decibar

    %
    spotterSmartdata.pressure = (spotterSmartdata.pressure ./ 1e4);


    %% Remove atmospheric pressure

    % From milibar to decibar
    atm_pressure = (atmpres_NOAA.atm_pres./100);

    % Based on SoloD at B18p (at China Rock, SN 77272) that was exposed
    % at low tide, there is a spatial variability in atmospheric pressure,
    % where the pressure at China Rock is slightly smaller than at the
    % Monterey harbor. Apply this correction for the spatial variability:
    atm_pressure = atm_pressure - 0.032;    % all in dbar

    % Interpolate atmospheric pressure to Smart mooring timestamps
    atm_pressure_interp = interp1(datenum(atmpres_NOAA.time_vec), atm_pressure, spotterSmartdata.dtime);

    % Remove atmospheric pressure from observations
    spotterSmartdata.pressure = spotterSmartdata.pressure - atm_pressure_interp;


    %% Compute average pressure in 30 min intervals
    % so that it gives a diagnostic view of the data
    
    %
    disp('---- Computing low-passed pressure data ---- ')

    %
    dtavg = 30/(24*60);

    %
    timeavg_vec = spotterSmartdata.dtime(1) : dtavg : spotterSmartdata.dtime(end);
    timeavg_vec = timeavg_vec(2:end-1);    % remove edges
    Pavg_vec = NaN(1, length(timeavg_vec));

    %
    tic
    for i2 = 1:length(timeavg_vec)
        %
        lin_time_avg = (spotterSmartdata.dtime >= (timeavg_vec(i2) - (dtavg/2))) & ...
                       (spotterSmartdata.dtime <= (timeavg_vec(i2) + (dtavg/2)));
        %
        pressure_inlims = spotterSmartdata.pressure(lin_time_avg);

        % Only compute average if there is enough data
        if length(pressure_inlims) > (0.9 * 30*60*2)
            %
            Pavg_vec(i2) = mean(pressure_inlims, 'omitnan');
        end
    end
    toc

    %
    disp('---- Done with computing low-passed pressure data ---- ')

% %     %
% %     figure
% %         plot(datetime(timeavg_vec, 'ConvertFrom', 'datenum'), Pavg_vec, '.-')


    %%
    % ----------------------------------------------------------
    % ------- ORGANIZE LEVEL 1 DATA STRUCTURE AND SAVE IT ------
    % ----------------------------------------------------------

    %
    spotsmart.mooringID = [list_SmartMoorings{i1}(1:3) 'sp'];
    %
    spotsmart.SN = list_SmartMoorings{i1}(end-3:end);

    %
    lmatch = strcmp(dplySpotters.SN, list_SmartMoorings{i1}(end-3:end));
    %
    spotsmart.latitude = mooringtable(lmatch, :).latitude;
    spotsmart.longitude = mooringtable(lmatch, :).longitude;

    % Height of the sensor above the bottom
    % (this was measured in inches, = 5 inches).
    spotsmart.zhab = 12.7 * 1e-2;    % in meters

    %
    spotsmart.dtime = datetime(spotterSmartdata.dtime, 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    %
    spotsmart.pressure = spotterSmartdata.pressure;

    %
    spotsmart.dtime_dt_avg = dtavg;
    %
    spotsmart.dtime_avg = datetime(timeavg_vec, 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    %
    spotsmart.pressure_avg = Pavg_vec;

    %
    spotsmart.time_zone = 'PDT';

    %
    time_dataproc = datestr(datetime('now', 'TimeZone', 'America/Los_Angeles'), 'yyyy/mm/dd HH:MM:SS');
    %
    spotsmart.REAMDE = ['Level 1 smart mooring data from ROXSI 2022. Data processed by script ' ...
                        mfilename() '.m on ' time_dataproc ' (PDT). ' ...
                        'Pressure is in decibar and atmospheric pressure was subtracted from ' ...
                        'the data. Longitude and latitude (for the pressure measurement) ' ...
                        'were computed from the mean watch circle obtained ' ...
                        'from Spotter coordinates over the whole deployment. ' ...
                        'zhab is the the height in meters of the pressure sensor above the bottom.'];

    %
    disp('---- Saving smart mooring level 1 data ---- ')

    % Save Level 1 data
    save(fullfile(repo_dirpath(), ['smart_mooring_' spotsmart.mooringID '_' spotsmart.SN '_L1.mat']), 'spotsmart');

    %
    disp(['-------------- Done with level 1 data processing of smart mooring ' spotsmart.mooringID ' - SN ' spotsmart.SN ' --------------'])


    %% Clear variables/close figures as needed before
    % going to the next loop iteration

% %     %
% %     close all


end


% If you want to save a figure, do something like this
% right after the code section that plots the figure:
% 
% exportgraphics(gcf, ['smartmooring_clockQC_timestops_' num2str(i2, '%.02d') '.png'], 'Resolution', 300);


