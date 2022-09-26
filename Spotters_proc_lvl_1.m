%% Script to convert the RAW Spotter data (in *.mat files)
% to Level 1. This script does:
%   - Trim edges.
%   - Grid some variables and deal with minor gaps.
%   - May recompute a1, a2, b1, and b2.
%   - Output in better format
%

clear
close all


%%


%
% dir_data_level_1 = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';
dir_data_parsed = '/Volumes/ROXSI_Data/LargeScale_Data_2022/RAW/Spotters/SDcards/';

 
% % All Spotters and Smart Moorings
% list_Spotters = {'B01_spot1150', 'B01_spot1158', 'B03_spot1152', ...
%               'B05_spot1153', 'X01_spot1151', 'X03_spot1157', ...
%               'X04_spot1155', ...
%               'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
%               'E07_spot1855', 'E07_spot1857', ...
%               'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
%               'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

% All Spotters and Smart Moorings
list_spotters = {'B01_spot1150', 'B01_spot1158', ...
                 'B03_spot1152', 'B05_spot1153', ...
                 'X01_spot1151', 'X03_spot1157', 'X04_spot1155', ...
                 'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
                 'E07_spot1855', 'E07_spot1857', ...
                 'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
                 'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};


%%

% Sampling period (in seconds)
dt = 0.4;    % this is 

% Bulk statistics period for which bulk
% statistics will be computed for (in minutes).
% Sofar does it hourly.
% dt_bulkstats = 30;
dt_bulkstats = 60;

% Parameters for bulk statistics,
% (which match what Sofar uses!)
nfft = 256;
% 
noverlap = nfft/2;
window = hanning(nfft);



%%

spotter_dplt = load(fullfile(repo_dirpath(), 'deploymentInfo_Spotters_ROXSI2022.mat'));
spotter_dplt = spotter_dplt.dployInfo_Spotters;


%%

list_tablefields = ["a1", "a2", "b1", "b2", "bulkparameters", "displacement", "location", "sst"];


%%
% ------------------------------------------------------------
% ------------------------------------------------------------
% ------------------------------------------------------------


%% First plot bulk parameters for all Spotters to
% check trimming edges

%%


%
for i = 3%1:length(list_spotters)

    %%
    data_aux = load(fullfile(dir_data_parsed, list_spotters{i}, 'parsed', [list_spotters{i} '.mat']));
    data_aux = data_aux.s;


    %
    ind_match = find(strncmp(spotter_dplt.SN, list_spotters{i}(9:12), 4));
    
    %%
    %
    figure
        plot(data_aux.bulkparameters.time, ...
             data_aux.bulkparameters.("Significant Wave Height"), '.-')

        %
        trim_edge_1 = datetime(datenum(spotter_dplt.time_begin_trim(ind_match), 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum');
        trim_edge_2 = datetime(datenum(spotter_dplt.time_end_trim(ind_match), 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum');

        %
        overlayline('v', trim_edge_1, '--k')
        overlayline('v', trim_edge_2, '--k')

        %
        xaxs_lims = [datetime(2022, 06, 15), datetime(2022, 07, 23)];
        xlim(xaxs_lims)

        %
        grid on
        set(gca, 'FontSize', 16, 'Box', 'on')

        %
        title([list_spotters{i}(1:3) '-' list_spotters{i}(9:12)], 'FontSize', 18)

        %
        set(gcf, 'Units', 'normalized')
        set(gcf, 'Position', [0.25, 0.25, 0.6, 0.4])


    %% Trim data removing data points before/after (or during) deployment/recovery

    %
    data_trimmed = data_aux;

    %
    for i2 = 1:length(list_tablefields)

        % Make sure 
        if isfield(data_aux, list_tablefields(i2))

            %
            lintrim_aux = (data_aux.(list_tablefields(i2)).time >= trim_edge_1) & ...
                          (data_aux.(list_tablefields(i2)).time <= trim_edge_2);
    
            %
            data_trimmed.(list_tablefields(i2)) = data_trimmed.(list_tablefields(i2))(lintrim_aux, :);

            %
            list_vars_intable_aux = data_aux.(list_tablefields(i2)).Properties.VariableNames;

            %
            for i3 = 1:length(list_vars_intable_aux)
                
                %
                if strcmp(list_vars_intable_aux{i3}, 'year') || strcmp(list_vars_intable_aux{i3}, 'month') || ...
                   strcmp(list_vars_intable_aux{i3}, 'day') || strcmp(list_vars_intable_aux{i3}, 'hour') || ...
                   strcmp(list_vars_intable_aux{i3}, 'min') || strcmp(list_vars_intable_aux{i3}, 'sec') || ...
                   strcmp(list_vars_intable_aux{i3}, 'msec') || strcmp(list_vars_intable_aux{i3}, 'milisec')
    
                    %
                    data_trimmed.(list_tablefields(i2)) = removevars(data_trimmed.(list_tablefields(i2)), list_vars_intable_aux{i3});

                end
            end
        %
        else
            warning(['Table ' list_tablefields(i2) ' does not exist in the data for '])
        %
        end            

    end

    %%

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


    %% Interpolate short (3s) gaps? deal with longer gaps?


    
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
        
% %         %
% %         mins_to_next_point = mod(2*30 - first_indata.Minute, 30);

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
    

    %%

% %     data_out = data_trimmed;

    %% Recalculate Fourier coefficients and bulkparameters
    % for specified time grid points -- here we neglect that
    % the data is NOT EXACTLY on a time grid (there are only
    % milisecond variations though, which are small and provided
    % by the Spotter)

    %
    data_out.timestats = time_grid_aux;
    %
    data_out.frequency = [];
    data_out.Ezz = NaN(128, length(time_grid_aux));

    % Loop over grid points, get the appropriate displacement
    % data, and compute bulk statistics
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
            warning(['Wrong number of data points!!! There are ' num2str(length(x_data_aux)) ' instead'])
        end

        %

        %
        [a1, a2, b1, b2, ...
              Ezz, Hsig, T_mean, ...
              dir_mean, spread_mean, f_peak, T_peak, dir_peak, spread_peak, ...
              f, DoF] = ...
                    wave_spec_fourier_displacement(x_data_aux, y_data_aux, z_data_aux, ...
                                                   dt, nfft, noverlap, window);

        % Sofar's coefficients are 0 : 0.009765625 : 1.240234375,
        % where the first 3 are NaNs. For default parameters, there
        % is one more (higher) frequency that comes out of
        % wave_spec_fourier_displacement (cpsd.m creates f inside of
        % the higher level function), but this is apparently
        % neglected by Sofar.

        % Put Spectra in output structure
        if i2==1
            data_out.frequency = f(1:end-1);
        end
        %
        data_out.Ezz(:, i2) = Ezz(1:end-1);

        % Put data in spectrum in output structure

        % Put Fourier coefficients in output structure

        % Put bulk statistics


    end
    toc

    %%

    % Pass displacement
    data_out.displacement = data_trimmed.displacement;

    % Pass location
    data_out.location = data_trimmed.location;

    % Pass sst (if it exists)
    if isfield(data_trimmed, "sst")
        data_out.sst = data_trimmed.sst;
    end

    
    %
    keyboard


end


%%
    
% %     %
% %     figure
% %         plot(data_aux.bulkparameters.time, ...
% %              data_aux.bulkparameters.("Significant Wave Height"), '.-')
% % 


%%
% ------------------------------------------------------------
% ------------------------------------------------------------
% ------------------------------------------------------------

% % % % From Sofar Spotter parsing script
% % 
% % The first columns indicate the time (year, month etc.) and dof is the 
% % degrees of freedom (dof) used to calculate the spectra. After 
% % the degrees of freedom, each subsequent entry corresponds to the variance 
% % density at the frequency indicated by the header line (E0 is the energy in
% % the mean, E1 at the first frequency f1 etc). The Spotter records
% % at an equidistant spectral resolution of df=0.009765625 and there are
% % nf=128 spectral entries, given by f(j) = df * j (with 0<=j<128). Frequencies are
% % in Hertz, and spectral entries are given in squared meters per Hz (m^2/Hz) or 
% % are dimensionless (for the directional moments a1,a2,b1,b2).



