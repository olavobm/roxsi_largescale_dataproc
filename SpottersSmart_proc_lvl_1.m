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

% All good Smart Moorings
list_SmartMoorings = {'E01_spot1851', 'E02_spot1859', 'E08_spot1852', 'E10_spot1848'};
% Same, but separately
% % list_SmartMoorings = {'E01_spot1851'};
% % list_SmartMoorings = {'E02_spot1859'};
% % list_SmartMoorings = {'E08_spot1852'};
% % list_SmartMoorings = {'E10_spot1848'};

% All seemingly great, apart from not extending the whole deployment
% list_SmartMoorings = {'E05_spot1853', 'E07_spot1857', 'E09_spot1856'};
%
% list_SmartMoorings = {'E05_spot1853'};
% list_SmartMoorings = {'E07_spot1857'};
% list_SmartMoorings = {'E09_spot1856'};


% % 
% list_SmartMoorings = {'E09_spot1850'};    % one annoying problem
% list_SmartMoorings = {'E11_spot1860'};

% % % Three that broke (though only E07 and E13 will have sparse segments)
% % list_SmartMoorings = {'E07_spot1855, 'E09_spot1856', 'E13_spot1849'};

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
%     NsegTH = 20;    % not good in one instance in E11 - SN 1860, with 2-step correction that is a bit longer than 20 pts 
    NsegTH = 30;


    %
    if any(diff(inds_gobacks) <= NsegTH)
        %
        warning([list_SmartMoorings{i1} ': THERE ARE CLOCK INVERSIONS THAT ARE TOO CLOSE TOGETHER!!!!'])
    end

    %
    for i2 = 1:length(inds_gobacks)

        %
        lfix_segment_aux = true;
% % % 
% % %         % Finds what is likely the corresponding time difference when
% % %         % the clock resets (HOWEVER, there are some underlying assumptions
% % %         % here).
% % %         [~, ind_back_normaltime_aux] = max(timediff_aux(inds_gobacks(i2):(inds_gobacks(i2) + NsegTH)));
% % % 
% % %         % Indices of the segment when the timestamps have gone back in time
% % %         ind_seg_tofix = inds_gobacks(i2) + (1 : (ind_back_normaltime_aux - 1));



        % ------------------------------------------
        %
        ind_seg_lookatdiff = inds_gobacks(i2): 1 : (inds_gobacks(i2) + NsegTH);
        %
        time_diff_lookat = 24*3600*diff(spotterSmartdata.dtime(ind_seg_lookatdiff));

        %
        time_diff_anomaly = time_diff_lookat - 0.5;

        %
        big_timediff_anomaly = time_diff_anomaly(abs(time_diff_anomaly) > 0.2);

        %
        big_timediff_anomaly_relative = time_diff_anomaly(abs(time_diff_anomaly) > 0.2);

        %
        if any(big_timediff_anomaly(2:end) < 0)
            warning('!!!!! CLOCK GOES BACKWARD TWICE IN TWO INSTANCES VERY CLOSE TO EACH OTHER !!!!!')
        end

        % This assumes that big_timediff_anomaly(2:end) are all positive
        lbig_enough = big_timediff_anomaly(2:end) > (0.3*abs(big_timediff_anomaly(1)));
        lbig_enough = [true; lbig_enough];

        %
        big_enough_timediff_anomaly_relative = big_timediff_anomaly_relative(lbig_enough);

% %         %
% %         if (inds_gobacks(i2) > 2492600) && (inds_gobacks(i2) < 2492680)
% %             keyboard
% %         end
% % 
% %         if inds_gobacks(i2)==2372744
% %             keyboard
% %         end
% %         if inds_gobacks(i2)==2372755
% %             keyboard
% %         end

% %         if i2 >= 291
% %             keyboard
% %         end


        %
        if (length(big_enough_timediff_anomaly_relative) == 1) && (big_enough_timediff_anomaly_relative > -0.5)
% %             disp('------- CLOCK SLOWED DOWN, BUT IT IS LIKELY REAL AND NO CORRECTION IS NECESSARY -------')
            lfix_segment_aux = false;

        %
        elseif length(big_enough_timediff_anomaly_relative) == 2
% %             disp('------- SIMPLE CLOCK INVERSION -------')

            %
            [~, ind_backtonormal_relative] = min(abs(time_diff_anomaly(2:end) + time_diff_anomaly(1)));
            ind_backtonormal = inds_gobacks(i2) + (ind_backtonormal_relative);
            %
            ind_seg_tofix = (inds_gobacks(i2)+1) : 1 : ind_backtonormal;

% %             % Plot to check that the identification of data points that
% %             % need to be adjusted is correct
% %             figure
% %                 grid on
% %                 hold on
% %                 inds_plt = (inds_gobacks(i2)-10):(ind_seg_tofix(end)+10);
% %                 plot(inds_plt, spotterSmartdata.dtime(inds_plt), '.-', 'MarkerSize', 20)
% %                 plot(ind_seg_tofix, spotterSmartdata.dtime(ind_seg_tofix), '.-', 'MarkerSize', 20)

        %
        elseif length(big_enough_timediff_anomaly_relative) == 3

            % Check if a 2-point correction may be good enough
            balance_between_anomalies = big_enough_timediff_anomaly_relative(2:end) + big_enough_timediff_anomaly_relative(1);
            %
            [min_balance, ind_balance] = min(abs(balance_between_anomalies));

            % Test whether 2 points make a reasonable (ad-hoc) correction,
            % i.e. balance each other/add up to something close to 0 (such
            % that we only have a simple clock inversion to correct)
            if min_balance < 0.2
                % Then only make a 2-point correction
% %                 disp('------- 3 TIMESTAMP ANOMALIES, BUT ONLY A SIMPLE CLOCK INVERSION -------')

                % Get the indices of data points that need to be fixed.
                % This is the same as for the simple inversion in the
                % elseif above, and it does (or should do) something
                % analogous to the min_balance calculated above
                [~, ind_backtonormal_relative] = min(abs(time_diff_anomaly(2:end) + time_diff_anomaly(1)));
                ind_backtonormal = inds_gobacks(i2) + (ind_backtonormal_relative);
                %
                ind_seg_tofix = (inds_gobacks(i2)+1) : 1 : ind_backtonormal;
    
% %                 % Plot to check that the identification of data points that
% %                 % need to be adjusted is correct
% %                 figure
% %                     grid on
% %                     hold on
% %                     inds_plt = (inds_gobacks(i2)-10):(ind_seg_tofix(end)+10);
% %                     plot(inds_plt, spotterSmartdata.dtime(inds_plt), '.-', 'MarkerSize', 20)
% %                     plot(ind_seg_tofix, spotterSmartdata.dtime(ind_seg_tofix), '.-', 'MarkerSize', 20)
    
            else

                % Use a criterion to infer if a 2-step clock correction
                % is reasonable. Otherwise do something different

                %
                if abs(sum(big_enough_timediff_anomaly_relative))<0.25
                    % If there is a reasonable balance between 3 anomalies,
                    % then do the 2-step clock correction

% %                     disp('------- NEEDS 2-step CLOCK CORRECTION -------')

                    % Here I'll just do the second part of the segment
                    % that needs to be corrected, turning into a
                    % simple clock inversion that will be addressed
                    % below just like the other cases in the if statement
                    
                    %
                    ind_relative_secondhalf_1 = find(time_diff_anomaly == big_enough_timediff_anomaly_relative(2));
                    ind_relative_secondhalf_2 = find(time_diff_anomaly == big_enough_timediff_anomaly_relative(3));
                    %
                    ind_secondhalf_1 = inds_gobacks(i2) + ind_relative_secondhalf_1;
                    ind_secondhalf_2 = inds_gobacks(i2) + ind_relative_secondhalf_2 - 1;
                    %
                    ind_secondhalf = ind_secondhalf_1 : ind_secondhalf_2;
    
                    % Make partial correction
                    spotterSmartdata.dtime(ind_secondhalf) = spotterSmartdata.dtime(ind_secondhalf) - big_enough_timediff_anomaly_relative(2)/(24*3600);
    
                    % Now create the indices of the full segment to
                    % be adjusted later like all the other cases
                    % in the if statement
                    ind_seg_tofix = (inds_gobacks(i2)+1) : 1 : ind_secondhalf(end);

                else
                    % Otherwise we have something different. Correct
                    % the inversion ONLY if there is a subsequent
                    % skip that is bigger than the inversion (thus
                    % the final correction is OK). Otherwise, throw
                    % an error.
                    %
                    max_anomaly = max(big_enough_timediff_anomaly_relative);
                    %
                    if max_anomaly > abs(big_enough_timediff_anomaly_relative(1))
                        % If the maximum anomaly is larger than the
                        % negative anomaly, then we can make a correction
% % keyboard
                    else
                        % Otherwise, the correction above will still have a
                        % clock reversal. Then doesn't do anything

                        % %                       lfix_segment_aux
                    end
                    
                    
                end



                % Or instead make a 3-point correction


% %                 % Plot to check that the identification of data points that
% %                 % need to be adjusted is correct
% %                 figure
% %                     grid on
% %                     hold on
% %                     inds_plt = (inds_gobacks(i2)-10):(ind_seg_tofix(end)+10);
% %                     plot(inds_plt, spotterSmartdata.dtime(inds_plt), '.-', 'MarkerSize', 20)
% %                     plot(ind_seg_tofix, spotterSmartdata.dtime(ind_seg_tofix), '.-', 'MarkerSize', 20)

            end

        %
        else
            % Anything can happen here. I'll will only implement
            % a correction for simple clock inversion. More complex
            % cases may need to be removed altogether.
            length(big_enough_timediff_anomaly_relative)
            warning('****** WEIRD STUFF ******')

            % Check if a 2-point correction may be good enough
            balance_between_anomalies = big_enough_timediff_anomaly_relative(2:end) + big_enough_timediff_anomaly_relative(1);
            ind_goodbalance_relative = find(abs(balance_between_anomalies) < 0.1);
            %
            if length(ind_goodbalance_relative) == 1
                % Then there is a good balance between 2 anomalies,
                % and it seems reasonable to apply the simple clock
                % inversion correction
                ind_insegment_balance = find(time_diff_anomaly == ...
                                             big_enough_timediff_anomaly_relative(1 + ind_goodbalance_relative));

                %
                ind_seg_tofix = (inds_gobacks(i2)+1) : 1 : (inds_gobacks(i2) + ind_insegment_balance - 1);

            else
                % Otherwise just throw a warning (maybe should be error)
                length(big_enough_timediff_anomaly_relative)
                warning('****** --------------- A VERY COMPLEX/WEIRD PROBLEM --------------- ******')

% %                 % Doesn't try to fix
% %                 lfix_segment_aux = false;

                % Select indice for a simple clock inversion --
                % this could be OK or a problem. Similar as the
                % first instance in the if statement, but take
                % indice with the best balance. I think this should be
                % generally fine though (as long as the later
                % timestamp anomalies are not related to the
                % clock inversion).

% %                 % First one
% %                 ind_insegment_balance = find(time_diff_anomaly == ...
% %                                              big_enough_timediff_anomaly_relative(1 + ind_goodbalance_relative(1)));

                %
                [~, ind_bestbalance] = min(abs(balance_between_anomalies));
                
                % Also get the second best balance
                ind_others = setdiff(1:length(balance_between_anomalies), ind_bestbalance);
                ind_secondbestbalance = find(balance_between_anomalies == min(balance_between_anomalies(ind_others)));

                %
                if abs(balance_between_anomalies(ind_bestbalance) - ...
                       balance_between_anomalies(ind_secondbestbalance)) < 1e-3
                    % If difference is too small, just take
                    % the one that appears earlier
                    ind_choice_aux = min([ind_bestbalance, ind_secondbestbalance]);
                else
                    % Otherwise takes the best balance
                    ind_choice_aux = ind_bestbalance;
                end
              
                %
                ind_insegment_balance = find(time_diff_anomaly == ...
                                             big_enough_timediff_anomaly_relative(1 + ind_choice_aux));
                %
                ind_seg_tofix = (inds_gobacks(i2)+1) : 1 : (inds_gobacks(i2) + ind_insegment_balance - 1);

                

% %                 % Plot to check that the identification of data points that
% %                 % need to be adjusted is correct
% %                 figure
% %                     grid on
% %                     hold on
% %                     inds_plt = (inds_gobacks(i2)-10):(ind_seg_tofix(end)+10);
% %                     plot(inds_plt, spotterSmartdata.dtime(inds_plt), '.-', 'MarkerSize', 20)
% %                     plot(ind_seg_tofix, spotterSmartdata.dtime(ind_seg_tofix), '.-', 'MarkerSize', 20)

            end
        end
        

% %         % A very useful diagnostic plot for the potential clock issues
% %         % identified above -- just make sure this is commented when
% %         % running the code as a whole and just make the plot for
% %         % the desired iterations in the loop
% % % %         ind_seg_plt = ind_seg_lookatdiff;
% %         ind_seg_plt = (ind_seg_lookatdiff(1) - 10):(ind_seg_lookatdiff(end) + 10);
% %         %
% %         figure
% %             %
% %             set(gcf, 'Units', 'normalized')
% %             set(gcf, 'Position', [0.47, 0.23, 0.19, 0.51])
% %             %
% %             haxs_1 = axes('Position', [0.1, 0.7, 0.8, 0.2]);
% %             haxs_2 = axes('Position', [0.1, 0.4, 0.8, 0.2]);
% %             haxs_3 = axes('Position', [0.1, 0.1, 0.8, 0.2]);
% %             %
% %             plot(haxs_1, ind_seg_plt(1:end-1) + 0.5, 24*3600*diff(spotterSmartdata.dtime(ind_seg_plt)) - 0.5, '.-')
% %             plot(haxs_2, ind_seg_plt, spotterSmartdata.pressure(ind_seg_plt), '.-')
% %             plot(haxs_3, spotterSmartdata.dtime(ind_seg_plt), spotterSmartdata.pressure(ind_seg_plt), '.-')
% % 
% %         %
% %         set([haxs_1, haxs_2, haxs_3], 'FontSize', 16, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
% %         %
% %         axis(haxs_1, 'tight')
% %         axis(haxs_2, 'tight')
% %         axis(haxs_3, 'tight')


        % ------------------------------------------


%         % Check the time differences before and after the segment
%         % are DEFINITELY CONSISTENT with an error in the clock.
%         % (one counter example, which I don't if ever happens would
%         % be a real clock delay in the middle of the segment)
%         time_diff_before_aux = spotterSmartdata.dtime(ind_seg_tofix(1)) - spotterSmartdata.dtime((ind_seg_tofix(1)-1));
%         time_diff_after_aux = spotterSmartdata.dtime((ind_seg_tofix(end)+1)) - spotterSmartdata.dtime(ind_seg_tofix(end));
%         %
%         time_diff_before_aux = time_diff_before_aux*24*3600;
%         time_diff_after_aux = time_diff_after_aux*24*3600;
%         
%         %
%         time_diff_inseg_aux = 24*3600*diff(spotterSmartdata.dtime(ind_seg_tofix));



% % % % % % % % % %  I THINK THIS WAS USEFUL ONLY FOR WHEN CLOCK SLOWS DOWN
% % % % % % % % % %  BUT THEN IT'S NOT A PROBLEM. I SHOULDN'T NEED THIS
% % % % % % % % % %  ANYMORE
% % % % 
% % % %         %
% % % %         if time_diff_before_aux > 0
% % % %             %
% % % %             if abs(time_diff_after_aux - 0.5) < 0.8*abs(time_diff_before_aux - 0.5)
% % % % % %                 %
% % % % % %                 warning('error 1!!!!!!!!')
% % % % % %                 
% % % % % %                 % From a couple of Spotters error 2 always happens
% % % % % %                 % when error 1 happens
% % % % % %                 if any(abs(time_diff_inseg_aux - 0.5) > abs(time_diff_after_aux - 0.5))
% % % % % %                     warning('error 2!!!!!!!!')
% % % % % % 
% % % % % %                 %
% % % % % %                 figure
% % % % % %                     plot((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+54), ...
% % % % % %                          24*3600*diff(spotterSmartdata.dtime((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+55))), '.-')
% % % % % %                     hold on
% % % % % %                     plot(ind_seg_tofix(1:end-1), ...
% % % % % %                          24*3600*diff(spotterSmartdata.dtime(ind_seg_tofix)), '.-')
% % % % % %                     grid on
% % % % % %                     %
% % % % % %                     set(gca, 'FontSize', 16)
% % % % % %                 
% % % % % %                 end
% % % % 
% % % %                 % Remove the segment in this iteration from
% % % %                 % those that need to be corrected
% % % %                 lfix_segment_aux = false;
% % % % 
% % % %             end
% % % % 
% % % %         %
% % % %         else 
% % % %             
% % % %         % Throw a warning if there is a clock problem within the segment
% % % %         % that it's getting fixed -- nice to check, but so far I haven't
% % % %         % seen an example of where this goes wrong
% % % % 
% % % % % %             %
% % % % % %             if abs(time_diff_inseg_aux - 0.5) >= 0.1
% % % % % %                 warning('error 3!!!!!!!!')
% % % % % %             end
% % % % % % 
% % % % % %             %
% % % % % %             vec_diffaux = timediff_aux(inds_gobacks(i2):(inds_gobacks(i2) + NsegTH));
% % % % % %             %
% % % % % %             max_timediff = max(vec_diffaux);
% % % % % %             %
% % % % % %             new_vec_diffaux = setdiff(vec_diffaux, max_timediff, 'stable');
% % % % % %             
% % % % % %             %
% % % % % %             second_max_timediff = max(new_vec_diffaux);
% % % % % % 
% % % % % %             %
% % % % % %             if abs(second_max_timediff - 0.5) >= 0.1
% % % % % %                 warning('error 4!!!!!!!!')
% % % % % %                 %
% % % % % %                 figure
% % % % % %                     plot((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+24), ...
% % % % % %                          24*3600*diff(spotterSmartdata.dtime((ind_seg_tofix(1)-8):(ind_seg_tofix(1)+25))), '.-')
% % % % % %                     hold on
% % % % % %                     plot(ind_seg_tofix(1:end-1), ...
% % % % % %                          24*3600*diff(spotterSmartdata.dtime(ind_seg_tofix)), '.-')
% % % % % %                     grid on
% % % % % %                     %
% % % % % %                     set(gca, 'FontSize', 16)
% % % % % % 
% % % % % %             end
% % % %         end

        % Applies correction if it's not an exception
        if lfix_segment_aux

            % Make a simple correction in terms of multiples
            % of the sampling period (0.5 s)
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

    %% Plot the unique values of clock finite difference
    % (to check my processing has removed all
    % instances of clock going backwards) -- not a very
    % insightful/useful plot though
    
% %     % ------------------------------------------

% %     hclockfig = figure;
% %         %
% %         time_diff_unique = unique(24*3600*diff(spotterSmartdata.dtime));
% %         %
% %         lskipstodealwith = abs(time_diff_unique) < 100; 
% %         
% %         %
% %         plot(time_diff_unique(lskipstodealwith), '.-')
% % 
% %         %
% %         grid on
% %         set(gca, 'FontSize', 16)
% % 
% %         %
% %         xlabel('Indice of unique time difference', 'Interpreter', 'Latex', 'FontSize', 16)
% %         ylabel('Time difference [s]', 'Interpreter', 'Latex', 'FontSize', 16)
% %         %
% %         title({['ROXSI 2022: ' list_SmartMoorings{i1}(1:3) ' SN ' ...
% %                list_SmartMoorings{i1}(end-3:end) ': unique ' ...
% %                '(short) time'];'differences in processed data'}, 'Interpreter', 'Latex', 'FontSize', 20)

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

% %     % Save Level 1 data
% %     save(fullfile(repo_dirpath(), ['smart_mooring_' spotsmart.mooringID '_' spotsmart.SN '_L1.mat']), 'spotsmart');

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


%% Plot to make sure there are no clock timestamps going backwards

%
figure
    %
    plot(24*3600*diff(datenum(spotsmart.dtime)), '.-')
    hold on
    
    %
    ylim([-0.5, 2])
    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'XGrid', 'on', 'YGrid', 'on')
    %
    xlim([-0.02, 1.02].*length(spotsmart.dtime))
    %
    plot(xlim, [0, 0], '-k')

    %
    title([list_SmartMoorings{1}(1:3) ' - SN ' list_SmartMoorings{1}(end-3:end)], 'FontSize', 22)


%%

return

%%

% % %
% % inds_lim_procdata = [2492600, 2492700];
% % % inds_lim_procdata = inds_lim_procdata + 1e6;

% For E08
inds_lim_procdata = [2372700, 2372790];

%
inds_segment_procdata = inds_lim_procdata(1):inds_lim_procdata(2);

% Get corresponding time limits
time_lims = datenum(spotsmart.dtime(inds_lim_procdata));

% Now find the corresponding indices in the "unprocessed" data
time_uncorrected = 719529 + (spotterSmartdata.unixEpoch./86400) - (7/24);

% Do this in a way that I don't think will fail 
% even with potential weird stuff in the clock
inds_lim_uncorrectedtime = [find( (time_uncorrected > (time_lims(1) - 3/(24*3600))) & (time_uncorrected < (time_lims(1) + 3/(24*3600))), 1, 'first'), ...
                            find( (time_uncorrected > (time_lims(2) - 3/(24*3600))) & (time_uncorrected < (time_lims(2) + 3/(24*3600))), 1, 'last')];

inds_segment_uncorrectedtime = inds_lim_uncorrectedtime(1):inds_lim_uncorrectedtime(2);


%
figure
    %
    set(gcf, 'Units', 'normalized')
    set(gcf, 'Position', [0.47, 0.23, 0.19, 0.51])
    %
    haxs_1 = axes('Position', [0.1, 0.6, 0.8, 0.3]);
    haxs_2 = axes('Position', [0.1, 0.2, 0.8, 0.3]);


    %
    plot(haxs_1, 0.5 + inds_segment_uncorrectedtime(1:end-1), 24*3600*diff(time_uncorrected(inds_segment_uncorrectedtime)), '.-')
    plot(haxs_2, 0.5 + inds_segment_procdata(1:end-1), 24*3600*diff(datenum(spotsmart.dtime(inds_segment_procdata))), '.-')
    
  
% %     %
% %     plot(haxs_1, datetime(time_uncorrected(inds_segment_uncorrectedtime(1:end-1)), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles'), ...
% %                  24*3600*diff(time_uncorrected(inds_segment_uncorrectedtime)), '.-')
% %     
% %     %
% %     plot(haxs_2, spotsmart.dtime(inds_segment_procdata(1:end-1)), ...
% %                  24*3600*diff(datenum(spotsmart.dtime(inds_segment_procdata))), '.-')


%
set([haxs_1, haxs_2], 'FontSize', 16, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
%
axis(haxs_1, 'tight')
axis(haxs_2, 'tight')
% % axis(haxs_3, 'tight')

%
linkallaxes('xy')




%%

figure
    hold on
    %
    plot(spotsmart.dtime, spotsmart.pressure, '.-')
    plot(datetime(719529 + (spotterSmartdata.unixEpoch./86400) - (7/24), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles'), ...
         spotterSmartdata.pressure, '.-')




