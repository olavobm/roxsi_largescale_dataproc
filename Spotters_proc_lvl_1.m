%% Script to convert the RAW Spotter data (in *.mat files)
% to Level 1. This script does:
%   - Trim edges.
%   - Grid some variables and deal with minor gaps.
%   - May recompute a1, a2, b1, and b2.
%   - Output in better format
%
% The RAW *.mat files were created by the script read_spotter_ROXSI.m,
% that you can find in the RAW data folder. It's read_spotter_ROXSI.m
% that converts time in UTC recorded by the Spotter to local time
% (-7 hours).


clear
close all


%%
% --------------------------------------
% -------- SET DIRECTORY PATHS ---------
% --------------------------------------


%%

% %
% dir_data_parent = "/Volumes/ROXSI_Data/LargeScale_Data_2022/RAW/";
dir_data_parent = '/project/CSIDE/ROXSI/LargeScale_Data_2022/RAW/';


%% Define output directories

%
% dir_outlvl1 = "/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1_new/";
dir_outlvl1 = pwd;
%
% dir_QCfig = fullfile(dir_outlvl1, 'figs_QC');
dir_QCfig = dir_outlvl1;


%%
% --------------------------------------
% --- DEFINE VARIABLES FOR PROCESSING --
% --------------------------------------


%% Directory of the RAW data and list of Spotters
% that will be processed

% See right below the dynamic definition
% % dir_data_level_1 = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';
% dir_data_parsed = '/Volumes/ROXSI_Data/LargeScale_Data_2022/RAW/Spotters/SDcards/';


% All Spotters and Smart Moorings
list_spotters = {'B01_spot1150', 'B01_spot1158', ...
                 'B03_spot1152', 'B05_spot1153', ...
                 'X01_spot1151', 'X03_spot1157', 'X04_spot1155', ...
                 'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
                 'E07_spot1855', 'E07_spot1857', ...
                 'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
                 'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

% % % Test on a few Spotters
% % list_spotters = {'B01_spot1150', 'B01_spot1158', ...
% %                  'B05_spot1153', ...
% %                  'E02_spot1859', 'E05_spot1853'};
% % list_spotters = {'B01_spot1150'};

%
dir_data_parsed = strings(length(list_spotters), 1);
%
for i = 1:length(list_spotters)
    %
    if strcmp(list_spotters{i}(1), 'E')
        dir_data_parsed(i) = fullfile(dir_data_parent, 'Spotters_Smart', 'SDcards');
    else
        dir_data_parsed(i) = fullfile(dir_data_parent, 'Spotters', 'SDcards');
    end
end


%% Parameters of the data and for computing bulk statistics
% (CODE HAS NOT BEEN THOROUGHLY TESTED FOR OTHER PARAMETERS!)

% Sampling period (in seconds)
dt = 0.4;

% Bulk statistics period for which bulk statistics will be
% computed for (in minutes). Sofar does it hourly.
% dt_bulkstats = 30;
dt_bulkstats = 60;

% Parameters for bulk statistics,
% (which match what Sofar uses!)
nfft = 256;
% 
noverlap = nfft/2;
window = hanning(nfft);

% Frequency limits for computing bulk statistics
% freq_lims = [0.029, 1.245];
freq_lims = [0.045, 0.3];


%% Load Spotter deployment info.

%
file_spotter_deployment = 'deploymentInfo_Spotters_ROXSI2022.mat';
%
spotter_dplt = load(fullfile(repo_dirpath(), file_spotter_deployment));
spotter_dplt = spotter_dplt.dployInfo_Spotters;


%% List of tables in the data

list_tables_totrim = ["a1", "a2", "b1", "b2", ...
                      "bulkparameters", "displacement", "location", "sst"];


%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
% ------------------------------------------------------------
% -------------------- DO DATA PROCESSING --------------------
% ------------------------------------------------------------

%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_Spotter_procL1_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_outlvl1, log_file_name))

%
totalRunTime = tic;


%%

disp(' '), disp(' ')
disp('------------------------------ Processing data from Spotters: ------------------------------')
for i = 1:length(list_spotters)
    disp([num2str(i) ' - ' list_spotters{i}])
end


%%

% Loop over Spotters
for i = 1:length(list_spotters)

    %
    disp(' '), disp(' ')
    disp(['---- Loading data from Spotter ' list_spotters{i} ' ----'])

    %%
    data_aux = load(fullfile(dir_data_parsed(i), list_spotters{i}, 'parsed', [list_spotters{i} '.mat']));
    data_aux = data_aux.s;

    %
    ind_match = find(strncmp(spotter_dplt.SN, list_spotters{i}(9:12), 4));
    
    %% A basic plot of Sofar's significant height
    %
    figure
        %
        hold on
        %
        plot(data_aux.bulkparameters.time, ...
             data_aux.bulkparameters.("Significant Wave Height"), '.-')

        %
        trim_edge_1 = datetime(datenum(spotter_dplt.time_begin_trim(ind_match), 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum');
        trim_edge_2 = datetime(datenum(spotter_dplt.time_end_trim(ind_match), 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum');

        %
        xaxs_lims = [datetime(2022, 06, 15), datetime(2022, 07, 23)];
        xlim(xaxs_lims)

        %
        grid on
        set(gca, 'FontSize', 16, 'Box', 'on')

        %
        ylim_aux = ylim;
        plot([trim_edge_1, trim_edge_1], ylim_aux, '--k')
        plot([trim_edge_2, trim_edge_2], ylim_aux, '--k')
        ylim(ylim_aux)

        %
        title([list_spotters{i}(1:3) '-' list_spotters{i}(9:12)], 'FontSize', 18)

        %
        set(gcf, 'Units', 'normalized')
        set(gcf, 'Position', [0.25, 0.25, 0.6, 0.4])


    %% Trim data removing data points before/after (or during) deployment/recovery

    %
    disp('----- Trimming data -----')

    %
    data_trimmed = data_aux;

    %
    for i2 = 1:length(list_tables_totrim)

        % Make sure 
        if isfield(data_aux, list_tables_totrim(i2))

            %
            lintrim_aux = (data_aux.(list_tables_totrim(i2)).time >= trim_edge_1) & ...
                          (data_aux.(list_tables_totrim(i2)).time <= trim_edge_2);
    
            %
            data_trimmed.(list_tables_totrim(i2)) = data_trimmed.(list_tables_totrim(i2))(lintrim_aux, :);

            %
            list_vars_intable_aux = data_aux.(list_tables_totrim(i2)).Properties.VariableNames;

            %
            for i3 = 1:length(list_vars_intable_aux)
                
                %
                if strcmp(list_vars_intable_aux{i3}, 'year') || strcmp(list_vars_intable_aux{i3}, 'month') || ...
                   strcmp(list_vars_intable_aux{i3}, 'day') || strcmp(list_vars_intable_aux{i3}, 'hour') || ...
                   strcmp(list_vars_intable_aux{i3}, 'min') || strcmp(list_vars_intable_aux{i3}, 'sec') || ...
                   strcmp(list_vars_intable_aux{i3}, 'msec') || strcmp(list_vars_intable_aux{i3}, 'milisec')
    
                    %
                    data_trimmed.(list_tables_totrim(i2)) = removevars(data_trimmed.(list_tables_totrim(i2)), list_vars_intable_aux{i3});

                end
            end
        %
        else
            warning(['Table ' char(list_tables_totrim(i2)) ' does not exist in the data for '])
        %
        end            

    end

    %% Rename table variables for consistency

    %
    if any(strcmp(data_trimmed.displacement.Properties.VariableNames, "x(m)"))
        data_trimmed.displacement = renamevars(data_trimmed.displacement, ...
                                               "x(m)", "x (m)");
    end

    %
    if any(strcmp(data_trimmed.displacement.Properties.VariableNames, "y(m)"))
        data_trimmed.displacement = renamevars(data_trimmed.displacement, ...
                                               "y(m)", "y (m)");
    end

    %
    if any(strcmp(data_trimmed.displacement.Properties.VariableNames, "z(m)"))
        data_trimmed.displacement = renamevars(data_trimmed.displacement, ...
                                               "z(m)", "z (m)");
    end


    %% Interpolate short gaps (according to a threshold set below)


    % Now interpolate short gaps -- actually this interpolation
    % creates timestamps with about 0.4s and linearly interpolate
    % the displacement quantities. These short gaps are only
    % interpolated because they shouldn't give significant errors
    % to the bulk statistics and we don't want to lose a bulk
    % statistics data point.

    %
    gapTH = 5;    % in seconds
    % Because the data is required to have about 0.4s sampling period,
    % gaps that are not a multple of 0.4s are a problem (which I will
    % not deal with at this point). So set another threshold to make
    % sure that the interpolation will create timestamps that are about
    % 0.4s apart. If beyong the TH, then gap won't be interpolated
    % and the code blocks below will skip the computation of bulk
    % statistics
    departuresamplingperiodTH = 0.002;    % in seconds

    %
    difftime_aux = diff(data_trimmed.displacement.time);

    % These are the indices of the data point just before the gap
    ind_short_gaps = find((difftime_aux > seconds(0.5)) & ...
                          (difftime_aux < seconds(gapTH)));
    %
    Ngaps = length(ind_short_gaps);
    Nfixedgaps = Ngaps;
        
    % Loop over short gaps
    if ~isempty(ind_short_gaps)

        %
        disp(['----- ' num2str(length(ind_short_gaps)) ' short gap(s) found ' ...
              '(<' num2str(gapTH) 's). Adding timestamps and NaN to ' ...
              'the displacement data over these gaps. -----'])

        %
        dummytable.time = data_trimmed.displacement.time;
        dummytable.x = data_trimmed.displacement.("x (m)");
        dummytable.y = data_trimmed.displacement.("y (m)");
        dummytable.z = data_trimmed.displacement.("z (m)");

        % Loop over short gaps
        for i2 = 1:length(ind_short_gaps)
            
            %
            inds_before_after_gap_aux = ind_short_gaps(i2) + [0, 1];
            %
            gap_in_seconds = 24*3600 * datenum(diff(dummytable.time(inds_before_after_gap_aux)));
    
            %
            npts_interp = (gap_in_seconds/0.4) - 1;
    
            % Round the number of points. If the gap is multiple,
            % then npts_interp is already an integer. If not, the
            % gap is not exactly a multiple of 0.4s, and the sampling
            % period TH will check below if it is close enough
            if abs(npts_interp - round(npts_interp))~=0
                npts_interp = round(npts_interp);
            end
    
            % Create timestamps (+2 to include the edges with the data points)
            timestamps_in_gap = linspace(dummytable.time(inds_before_after_gap_aux(1)), ...
                                         dummytable.time(inds_before_after_gap_aux(2)), ...
                                         (npts_interp + 2));
            timestamps_in_gap = timestamps_in_gap(:);
    
            % Check the sampling period TH, and
            % do or don't do the interpolation
            if abs(0.4 - (24*3600*datenum(diff(timestamps_in_gap(1:2))))) > departuresamplingperiodTH
                %
                warning('Short gap length is not a multiple of sampling period. Skipping interpolation.')
                %
                Nfixedgaps = Nfixedgaps - 1;
            else
    
                %
                x_ongap_aux = interp1(dummytable.time(inds_before_after_gap_aux), ...
                                      dummytable.x(inds_before_after_gap_aux), ...
                                      timestamps_in_gap);
                %
                y_ongap_aux = interp1(dummytable.time(inds_before_after_gap_aux), ...
                                      dummytable.y(inds_before_after_gap_aux), ...
                                      timestamps_in_gap);
                %
                z_ongap_aux = interp1(dummytable.time(inds_before_after_gap_aux), ...
                                      dummytable.z(inds_before_after_gap_aux), ...
                                      timestamps_in_gap);

                % Insert the interpolated gap in the data
                dummytable.time = [dummytable.time(1:inds_before_after_gap_aux(1)); timestamps_in_gap(2:end-1); dummytable.time(inds_before_after_gap_aux(2):end)];
                dummytable.x = [dummytable.x(1:inds_before_after_gap_aux(1)); x_ongap_aux(2:end-1); dummytable.x(inds_before_after_gap_aux(2):end)];
                dummytable.y = [dummytable.y(1:inds_before_after_gap_aux(1)); y_ongap_aux(2:end-1); dummytable.y(inds_before_after_gap_aux(2):end)];
                dummytable.z = [dummytable.z(1:inds_before_after_gap_aux(1)); z_ongap_aux(2:end-1); dummytable.z(inds_before_after_gap_aux(2):end)];

                % If there are more short gaps, then change their
                % indices because the size of the data has just
                % been changed
                if i2 < length(ind_short_gaps)
                    ind_short_gaps((i2+1):end) = ind_short_gaps((i2+1):end) + npts_interp;
                end
    
            end
            
    
        end
    end


    % First plot the diff(time) in the data and new
    % diff(time) in case interpolation was done
    hfig_aux = figure;
        %
        if ~isempty(ind_short_gaps)
            haxs_1 = axes('Position', [0.1, 0.55, 0.8, 0.35]);
            haxs_2 = axes('Position', [0.1, 0.10, 0.8, 0.35]);
            %
            haxs_all = [haxs_1, haxs_2];

            %
            plot(haxs_2, dummytable.time(1:end-1), 24*3600*datenum(diff(dummytable.time)))

            %
            title(haxs_2, ['diff(time) after interpolation over gap(s). ' ...
                           num2str(Ngaps) ' short gaps found and ' num2str(Nfixedgaps) ' were fixed.'], ...
                          'Interpreter', 'Latex', 'FontSize', 18)

        else
            haxs_1 = axes('Position', [0.15, 0.15, 0.7, 0.7]);
            %
            haxs_all = haxs_1;
        end

        %
        plot(haxs_1, data_trimmed.displacement.time(1:end-1), ...
                     24*3600*datenum(diff(data_trimmed.displacement.time)))

        %
        set(haxs_all, 'FontSize', 16, 'Box', 'on', ...
                      'XGrid', 'on', 'YGrid', 'on', ...
                      'XLim', data_trimmed.displacement.time([1, end]) + [-hours(1); +hours(1)]);
        linkaxes(haxs_all, 'xy')

        %
        ylabel(haxs_all, '[seconds]', 'Interpreter', 'Latex', 'FontSize', 18)
        %
        title(haxs_1, ['diff(time) in the loaded data. ' ...
                       'Spotter ' list_spotters{i}(1:3) ' - SN ' ...
                       list_spotters{i}(9:12) ' (' num2str(Ngaps) ' ' ...
                       'short gap(s) found)'], ...
                       'Interpreter', 'Latex', 'FontSize', 18)
       
        %
        set(hfig_aux, 'Units', 'normalized')
        set(hfig_aux, 'Position', [0.4227, 0.2632, 0.3773, 0.4063])
       
    %
    exportgraphics(hfig_aux, fullfile(dir_QCfig, [list_spotters{i} '_difftime_displacement_data.png']), 'Resolution', 300)


    % Copy interpolated data to the data structure and keep
    % the original names
    if ~isempty(ind_short_gaps)
        data_trimmed.displacement = struct2table(dummytable);
        data_trimmed.displacement = renamevars(data_trimmed.displacement, ...
                                               data_trimmed.displacement.Properties.VariableNames, ...
                                               ["time", "x (m)", "y (m)", "z (m)"]);
    end



    
    %% Define time grid to calculate bulk statistics for the
    % specific Spotter (the grid passes on whole hours for
    % simplicity)

    % These are the first and last timestamps that could work as grid
    % points (where there is sufficient data at the edges). Since these
    % timestamps are very likely not "on whole hours", the appropriate
    % grid points need to be taken as the timese adjacent to these
    first_indata = data_trimmed.bulkparameters.time(1) + minutes(dt_bulkstats/2);
    last_indata = data_trimmed.bulkparameters.time(end) - minutes(dt_bulkstats/2);

    % Find first and last grid points that are contained in the interval
    % [first_time_indata, last_time_indata] -- a possibly too conservative
    % requirement, but that should be OK.
    %
    % To expand this, you would do something similar, but defining the
    % FFT periods, and trimming based on that instead of the bulk
    % statistics time grid points
    %
    % For the first grid point
    if (first_indata.Second==0) && (mod(first_indata.Minute, dt_bulkstats)==0)
        % This is highly unlikely, but still possible
        first_time_gridpoint = first_indata;
    else
        
        %
        mins_to_next_point = mod(60 - first_indata.Minute, dt_bulkstats);

        %
        if mins_to_next_point==0
            add_minute_factor = dt_bulkstats;
        else
            add_minute_factor = mins_to_next_point;
        end

        %
        first_time_gridpoint = first_indata;
        %
        first_time_gridpoint.Second = 0;
        first_time_gridpoint.Minute = first_time_gridpoint.Minute + add_minute_factor;
    end
    %
    % For the last grid point
    if (last_indata.Second==0) && (mod(last_indata.Minute, dt_bulkstats)==0)
        last_time_gridpoint = last_indata;
    else
        
        %
        mins_to_previous_point = mod(last_indata.Minute - 60, dt_bulkstats);
        %
        if mins_to_next_point==0
            minus_minute_factor = 0;
        else
            minus_minute_factor = mins_to_previous_point;
        end

        %
        last_time_gridpoint = last_indata;
        %
        last_time_gridpoint.Second = 0;
        last_time_gridpoint.Minute = last_time_gridpoint.Minute - minus_minute_factor;
    end

    %
    time_grid_aux = first_time_gridpoint : minutes(dt_bulkstats) : last_time_gridpoint;
    

    %% Dummy (short) time grid just to test the code
    %
% % %     time_grid_aux = datetime(2022, 07, 01, 0, 0, 0) : hours(1) : datetime(2022, 07, 03, 0, 0, 0);


% % %     % To match B03 time grid on the minutes for a more apples-to-apples
% % %     % comparison
% % %     time_grid_aux = time_grid_aux - minutes(15);
% % %     time_grid_aux = time_grid_aux(2:end);
% % % 
% % %     % To match with Sofar's bulk statistics, which
% % %     % are calculated from the preceding (1 h) data
% % %     % before the timestamp
% % %     time_grid_aux = time_grid_aux - minutes(30);
% % %     time_grid_aux = time_grid_aux(2:end);

    %%

% %     data_out = data_trimmed;

    %% Recalculate Fourier coefficients and bulk parameters
    % for specified time grid points -- here we neglect that
    % the data is NOT EXACTLY on a time grid (there are only
    % milisecond variations though, which are small and provided
    % by the Spotter. But if there are gaps, it would be a problem)

    %
    disp('----- Recalculating spectra and bulk statistics -----')
     
    %
    prealloc_aux = NaN(length(time_grid_aux), 128);

    %
    data_out.spectra.dt = diff(time_grid_aux(1:2));
    data_out.spectra.dtime = time_grid_aux(:);
    %
    data_out.spectra.df = 1/(dt*nfft);
    data_out.spectra.frequency = [];
    %
    data_out.spectra.Ezz = prealloc_aux;
    data_out.spectra.See = [];
    % DOF assuming 50% and subtracting 2 because the overlap
    % makes the chunks not entirely independent
    data_out.spectra.DOF = 2*(2*floor((60*dt_bulkstats/0.4)/nfft) - 1) - 2;

    %
    a1_matrix_dummy = prealloc_aux;
    a2_matrix_dummy = prealloc_aux;
    b1_matrix_dummy = prealloc_aux;
    b2_matrix_dummy = prealloc_aux;
    %
    bulkpars_matrix_dummy = NaN(length(time_grid_aux), 9);

    % Loop over time grid points, get the appropriate
    % displacement data, and compute bulk statistics
    %
    tic
    for i2 = 1:length(time_grid_aux)

        %
        time_edge_grid_1 = (time_grid_aux(i2) - minutes(dt_bulkstats/2));
        time_edge_grid_2 = (time_grid_aux(i2) + minutes(dt_bulkstats/2));

        %
        lintime_forbulkstats = (data_trimmed.displacement.time >= time_edge_grid_1) & ...
                               (data_trimmed.displacement.time < time_edge_grid_2);

        %
        x_data_aux = data_trimmed.displacement.("x (m)")(lintime_forbulkstats, :);
        y_data_aux = data_trimmed.displacement.("y (m)")(lintime_forbulkstats, :);
        z_data_aux = data_trimmed.displacement.("z (m)")(lintime_forbulkstats, :);

        %
        if length(x_data_aux)~=9000
            warning(['Wrong number of displacement data points! There ' ...
                     'are ' num2str(length(x_data_aux)) ' instead of 9000, ' ...
                     'which is the correct number for a sampling period of ' ...
                     '0.4 seconds. Skipping time grid point ' datestr(time_grid_aux(i2)) ' ' ...
                     'for spotter ' list_spotters{i}])

        end

        % Compute wave statistics from Spotter displacement
        [f, a1, a2, b1, b2, ...
         Ezz, Hsig, ...
         T_mean, dir_mean, spread_mean, ...
         T_peak, dir_peak, spread_peak, ...
         a1_bar, b1_bar, DOF] = ...
                    wave_spec_fourier_displacement(x_data_aux, y_data_aux, z_data_aux, ...
                                                   dt, nfft, noverlap, window, freq_lims);

        % Sofar's frequencies are 0 : 0.009765625 : 1.240234375,
        % where the first 3 are NaNs. For default parameters, there
        % is one more (higher) frequency that comes out of
        % wave_spec_fourier_displacement (cpsd.m creates f inside of
        % the higher level function), but this is apparently
        % neglected by Sofar.

        % Put Spectra in output structure
        if i2==1
            data_out.spectra.frequency = f(1:end-1);
        end
        %
        data_out.spectra.Ezz(i2, 4:end) = Ezz(4:end-1);

        % Put Fourier coefficients in output structure
        a1_matrix_dummy(i2, 4:end) = a1(4:end-1);    % starts at 4 (frequencies higher than 0.029 Hz) because that's what Sofar does.
        a2_matrix_dummy(i2, 4:end) = a2(4:end-1);
        b1_matrix_dummy(i2, 4:end) = b1(4:end-1);
        b2_matrix_dummy(i2, 4:end) = b2(4:end-1);
        
        % Put bulk statistics in dummy array
        bulkpars_matrix_dummy(i2, :) = [Hsig, ...
                                        T_mean, T_peak, ...
                                        dir_mean, dir_peak, ...
                                        spread_mean, spread_peak, ....
                                        a1_bar, b1_bar];
        
    end
    toc
    

    %% Convert arrays to tables

% %     %
% %     data_out.a1 = array2table(a1_matrix_dummy);
% %     data_out.a2 = array2table(a2_matrix_dummy);
% %     data_out.b1 = array2table(b1_matrix_dummy);
% %     data_out.b2 = array2table(b2_matrix_dummy);

    %
% % % %     data_out.spectra = array2table([time_grid_aux(:), bulkpars_matrix_dummy]);    % this would be Sofar's format 
% %     data_out.spectra = array2table(bulkpars_matrix_dummy);
% % 
% %     % Rename variables in the table of bulk parameters
% %     data_out.spectra = renamevars(data_out.spectra, ...
% %                                          data_out.spectra.Properties.VariableNames, ...
% %                           ["Significant Wave Height", "Mean Period", "Peak Period", ...
% %                            "Mean Direction", "Peak Direction", "Mean Spreading", "Peak Spreading"]);


    % Keep it as an array instead
    data_out.spectra.a1 = a1_matrix_dummy;
    data_out.spectra.a2 = a2_matrix_dummy;
    data_out.spectra.b1 = b1_matrix_dummy;
    data_out.spectra.b2 = b2_matrix_dummy;

    %
    data_out.spectra.freqband = freq_lims;

    %
    data_out.spectra.Hsig = bulkpars_matrix_dummy(:, 1);
    %
    data_out.spectra.Tmean = bulkpars_matrix_dummy(:, 2);
    data_out.spectra.meandir = bulkpars_matrix_dummy(:, 4);
    data_out.spectra.meandirspread = bulkpars_matrix_dummy(:, 6);
    %
    data_out.spectra.Tpeak = bulkpars_matrix_dummy(:, 3);
    data_out.spectra.peakdir = bulkpars_matrix_dummy(:, 5);
    data_out.spectra.peakdirspread = bulkpars_matrix_dummy(:, 7);
    %
    data_out.spectra.a1_bar = bulkpars_matrix_dummy(:, 8);
    data_out.spectra.b1_bar = bulkpars_matrix_dummy(:, 9);


    %% Trim out zero frequency and the first two non-zero frequencies
    % that Sofar does not give in their processing (though the code
    % above can compute them)

    %
    data_out.spectra.frequency = f(4:end);

    %
    data_out.spectra.a1 = data_out.spectra.a1(:, 4:end);
    data_out.spectra.a2 = data_out.spectra.a2(:, 4:end);
    data_out.spectra.b1 = data_out.spectra.b1(:, 4:end);
    data_out.spectra.b2 = data_out.spectra.b2(:, 4:end);
    %
    data_out.spectra.See = data_out.spectra.Ezz(:, 4:end);

    %
    data_out.spectra = rmfield(data_out.spectra, 'Ezz');


    %% Rename variables to be consistent with Sofar

%     % 9 decimal points is the maximum nonzero number of nonzero digits
%     names_frequency_char = num2str(f, '%.10f');    % First convert to char -- it's a matrix!
%     %
%     names_frequency_string = strings(1, length(data_out.frequency));    
% 
%     %
%     for i2 = 1:length(data_out.frequency)
%         % First deal with 0 frequency
%         if strcmp(names_frequency_char(i2, :), num2str(0, '%.10f'))
%             names_frequency_string(i2) = "0.0";
%         %
%         else
%             %
%             inds_loc_0 = strfind(names_frequency_char(i2, :), '0');
%             %
%             if length(inds_loc_0)==1
%                 names_frequency_string(i2) = convertCharsToStrings(names_frequency_char(i2, 1:(inds_loc_0(end)-1)));
%             else
% 
%                 %
%                 diff_inds_0 = diff(inds_loc_0);
%                 %
%                 ind_first_final_0 = find(diff_inds_0 > 1, 1, 'last') + 1;
% 
%                 %
%                 if isempty(ind_first_final_0)
%                     ind_first_final_0 = 1;
%                 end
% 
%                 %
%                 names_frequency_string(i2) = ...
%                                 convertCharsToStrings(names_frequency_char(i2, 1:(inds_loc_0(ind_first_final_0) - 1)));
%             end
%             
%             
%         end
%     end
% 
%     % Rename variables in the tables of coeffcients
%     data_out.a1 = renamevars(data_out.a1, data_out.a1.Properties.VariableNames, names_frequency_string);
%     data_out.a2 = renamevars(data_out.a2, data_out.a2.Properties.VariableNames, names_frequency_string);
%     data_out.b1 = renamevars(data_out.b1, data_out.b1.Properties.VariableNames, names_frequency_string);
%     data_out.b2 = renamevars(data_out.b2, data_out.b2.Properties.VariableNames, names_frequency_string);



    %% PROBABLY NOT SUPER USEFUL, BUT IT MIGHT BE GOOD TO HAVE
    % THE HOURLY AVERAGE LOCATION TO GO ALONG WITH THE BULK STATISTICS
    

    %% Pass other tables to output structure

    % Pass displacement
%     data_out.displacement = data_trimmed.displacement;
    data_out.displacement.dtime = data_trimmed.displacement.time;
    data_out.displacement.x = data_trimmed.displacement.("x (m)");
    data_out.displacement.y = data_trimmed.displacement.("y (m)");
    data_out.displacement.z = data_trimmed.displacement.("z (m)");

    % Pass location
%     data_out.location = data_trimmed.location;
    data_out.location.dtime = data_trimmed.location.time;
    data_out.location.latitude = data_trimmed.location.("latitude (decimal degrees)");
    data_out.location.longitude = data_trimmed.location.("longitude (decimal degrees)");

    % Pass sst (if it exists)
    if isfield(data_trimmed, "sst")
% %         data_out.sst = data_trimmed.sst;
% %         data_out.sst.time.TimeZone = 'America/Los_Angeles';
        data_out.SST.dtime = data_trimmed.sst.time;
        data_out.SST.dtime.TimeZone = 'America/Los_Angeles';
        %
        data_out.SST.temperature = data_trimmed.sst.("T (deg. Celcius)");
    end


    %% Add timezones to other datetime variables

% %     %
% %     data_out.timestats.TimeZone = 'America/Los_Angeles';
% %     data_out.displacement.time.TimeZone = 'America/Los_Angeles';
% %     data_out.location.time.TimeZone = 'America/Los_Angeles';

    %
    data_out.spectra.dtime.TimeZone = 'America/Los_Angeles';
    data_out.displacement.dtime.TimeZone = 'America/Los_Angeles';
    data_out.location.dtime.TimeZone = 'America/Los_Angeles';


    %% Create and organize output structure (including metadata)

    %
    spotterL1.mooringID = list_spotters{i}(1:3);
    spotterL1.SN = list_spotters{i}(9:12);
    
    %
    loc_mooring = ROXSI_mooringlocation(list_spotters{i}(1:3));
    %
    spotterL1.site = loc_mooring.roxsiarray;

    % Instead of the reference value that
    % ROXSI_mooringlocation.m gives, use the data
    spotterL1.latitude = mean(data_out.location.latitude);
    spotterL1.longitude = mean(data_out.location.longitude);
    %
    [spotterL1.X, ...
     spotterL1.Y] = ROXSI_lltoxy(spotterL1.latitude, ...
                                 spotterL1.longitude, spotterL1.site);

    %
    list_fields_aux = fieldnames(data_out);

    %
    for i2 = 1:length(list_fields_aux)
        %
        spotterL1.(list_fields_aux{i2}) = data_out.(list_fields_aux{i2});

% %         %
% %         if strcmp(list_fields_aux{i2}, 'Ezz')
% %             %
% %             npts_perbulkstats = dt_bulkstats * 60/0.4;
% % 
% %             % Degrees of freedom with ffts computed with 50% overlap
% %             % (remove 2 just because the overlap makes the chunks
% %             % not entirely independent)
% %             spotterL1.DOF = 2*(2*floor(npts_perbulkstats/nfft) - 1) - 2;
% %         end
    end
 
    %
    time_dataproc = datetime('now', 'TimeZone', 'Local');
    time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');
    %
    spotterL1.README = ['Level 1 Spotter data from ROXSI 2022. Data processed by script ' ...
                        mfilename() '.m on ' time_dataproc_char ' (TimeZone ' time_dataproc.TimeZone '). ' ...
                        'Data in the Level 1 structure has been trimmed for the deployment ' ...
                        'period, which is defined in the table at ' file_spotter_deployment '. ' ...
                        'See is the spectra of vertical displacement of the sea surface (in m2/Hz) ' ...
                        'and was computed as analogous as possible to Sofar''s calculations. These ' ...
                        'spectra, and all corresponding statistics, are calculated from ' ...
                        'data centered at the timestamps of the corresponding time vector (dtime).'];


    %% Make a plot of some of the basic variables
    %
    % Plot vertical displacement, significant wave height (both Sofar's
    % and computed by this code), and mean period.

    %
    hfig_aux = figure;
        %
        haxs_1 = axes('Position', [0.15, 0.7, 0.7, 0.15]);
        haxs_2 = axes('Position', [0.15, 0.5, 0.7, 0.15]);
        haxs_3 = axes('Position', [0.15, 0.3, 0.7, 0.15]);
        haxs_4 = axes('Position', [0.15, 0.1, 0.7, 0.15]);
        %
        haxs_all = [haxs_1, haxs_2, haxs_3, haxs_4];
        %
        hold(haxs_all, 'on')

            %
            time_raw_displacement = data_aux.displacement.time;
            time_raw_displacement.TimeZone = spotterL1.displacement.dtime.TimeZone;
            %
            time_raw_stats = data_aux.bulkparameters.time;
            time_raw_stats.TimeZone = spotterL1.spectra.dtime.TimeZone;

            % Plot vertical displacement
            if any(strcmp(data_aux.displacement.Properties.VariableNames, "z(m)"))
                field_z_aux = "z(m)";
            else
                field_z_aux = "z (m)";
            end
            plot(haxs_1, time_raw_displacement, data_aux.displacement.(field_z_aux))
            
            % Plot significant wave hieght
            plot(haxs_2, time_raw_stats, data_aux.bulkparameters.("Significant Wave Height"), '.-')
            plot(haxs_2, spotterL1.spectra.dtime, spotterL1.spectra.Hsig, '.-k')
            
            % Plot mean period
            plot(haxs_3, time_raw_stats, data_aux.bulkparameters.("Mean Period"), '.-')
            plot(haxs_3, spotterL1.spectra.dtime, spotterL1.spectra.Tmean, '.-k')
            % Plot mean direction
            plot(haxs_4, time_raw_stats, data_aux.bulkparameters.("Mean Direction"), '.-')
            plot(haxs_4, spotterL1.spectra.dtime, spotterL1.spectra.meandir, '.-k')

        %
        time_trim_1 = trim_edge_1;
        time_trim_2 = trim_edge_2;
        %
        time_trim_1.TimeZone = spotterL1.spectra.dtime.TimeZone;
        time_trim_2.TimeZone = spotterL1.spectra.dtime.TimeZone;
        %
        time_lims_plt = [time_trim_1; time_trim_2] + [-hours(24); +hours(24)];

        %
        set(haxs_all, 'FontSize', 16, 'Box', 'on', ...
                      'XGrid', 'on', 'YGrid', 'on', ...
                      'XLim', time_lims_plt);
        set(haxs_all(1:end-1), 'XTickLabel', [])
        % Set ylims of plots 2-4
        set(haxs_all(2), 'YLim', [0, 1.1*max(spotterL1.spectra.Hsig)])
        set(haxs_all(3), 'YLim', [min(spotterL1.spectra.Tmean), max(spotterL1.spectra.Tmean)])
        set(haxs_all(4), 'YLim', [min(spotterL1.spectra.meandir), max(spotterL1.spectra.meandir)])
        %
        linkaxes(haxs_all, 'x')

        % Plot trimming edges with dashed lines
        for i2 = 1:length(haxs_all)
            %
            ylims_aux = get(haxs_all(i2), 'YLim');
            %
            plot(haxs_all(i2), [time_trim_1, time_trim_1], ylims_aux, '--r')
            plot(haxs_all(i2), [time_trim_2, time_trim_2], ylims_aux, '--r')
            %
            set(haxs_all(i2), 'YLim', ylims_aux)
        end


        %
        ylabel(haxs_1, '[m]', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel(haxs_2, '[m]', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel(haxs_3, '[seconds]', 'Interpreter', 'Latex', 'FontSize', 18)
        ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 18)
        %
        title(haxs_1, {['ROXSI 2022: Spotter ' list_spotters{i}(1:3) ...
                       ' SN ' list_spotters{i}(9:12) ''];['Vertical ' ...
                       'displacement, $H$sig, mean period, and mean direction']; ...
                       '(processed data in black)'}, ...
                       'Interpreter', 'Latex', 'FontSize', 18)
       
        %
        set(hfig_aux, 'Units', 'normalized')
        set(hfig_aux, 'Position', [0.5023, 0.0778, 0.2844, 0.3917])



    %% Save the data and QC figure

    % -----------------------------------
    % Save QC figure
    %
    disp('----- Save level 1 data plot at: -----')
    disp(dir_QCfig)
    exportgraphics(hfig_aux, fullfile(dir_QCfig, ['spotter_L1_' list_spotters{i} '_data.png']), 'Resolution', 300)


    % -----------------------------------
    % Save data
    %
    str_filename = ['roxsi_spotter_L1_' list_spotters{i}(1:3) '_' list_spotters{i}(9:12)];
    str_fullpath_file = fullfile(dir_outlvl1, [str_filename '.mat']);
    %
    disp(['----- Saving level 1 data from Spotter ' list_spotters{i} ' at:-----'])
    disp(str_fullpath_file)

    %
    save(str_fullpath_file, 'spotterL1', '-v7.3')


    % -----------------------------------
    %
    disp(['----------------- Done with level 1 data processing for Spotter ' list_spotters{i} ' -----------------'])

    keyboard
    %
    clear spotterL1


end


%
disp('###################### Done with data processing for all Spotters ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');


