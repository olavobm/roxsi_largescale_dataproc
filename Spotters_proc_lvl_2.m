%% Compute directional Spectra for Spotters in ROXSI 2022
% 
%


clear
close all


%%
% --------------------------------------
% -------- SET DIRECTORY PATHS ---------
% --------------------------------------


%%

%
% % dir_data_level_1 = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';
dir_data_level_1 = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';


%%

% Output directory
% % dir_output_level_2 = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level2_Data/Spotter_Level2/';
% % dir_output_level_2 = pwd;
dir_output_level_2 = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level2_Data/Spotter_Level2/';


%% Add WAFO toolbox to Matlab path

% % %
% % % addpath(genpath(fullfile(repo_dirpath(), 'wafo')))
% % addpath(genpath(fullfile(repo_dirpath(), 'Spotter_DirectionalSpectra', 'wafo')))
%
roxsi_add_libraries()


%%
% -------------------------------------------
% --- DEFINE VARIABLES FOR DATA PROCESSING --
% -------------------------------------------

%%



 
% % All Spotters and Smart Moorings
% list_Spotters = {'B01_spot1150', 'B01_spot1158', 'B03_spot1152', ...
%               'B05_spot1153', 'X01_spot1151', 'X03_spot1157', ...
%               'X04_spot1155', ...
%               'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
%               'E07_spot1855', 'E07_spot1857', ...
%               'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
%               'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

% All Spotters and Smart Moorings
list_Spotters = {'B03_spot1152', 'B05_spot1153', ...
                 'X01_spot1151', 'X03_spot1157', 'X04_spot1155', ...
                 'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
                 'E07_spot1855', 'E07_spot1857', ...
                 'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
                 'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

% Just a few to test
% % list_Spotters = {'B03_spot1152'};
% % list_Spotters = {'B03_spot1152', 'B05_spot1153', 'E02_spot1859'};
% % list_Spotters = {'E02_spot1859'};
% list_Spotters = {'X03_spot1157'};

% list_Spotters = {'B01_spot1150', 'B01_spot1158'};
list_Spotters = {'B01_spot1158'};


%% Load the location file with bathymetry
% (actually, elevation to mean sea level)

spotter_location = load(fullfile(dir_data_level_1, 'Spotter_all_location_depth.mat'));
spotter_location = spotter_location.spotter_location;

% % %% Load trimming table for Spotters
% % %
% % % The trimming for Smart Moorings was defined in terms of
% % % pressure data, which had some issues for a few
% % % Smart Moorings. Looking at the definition, they work
% % % well for the Spotter buoy data, except for E09-SN1850.
% % % The end trim for this one should be 06/24-06:00:00.
% % 
% % 
% % %
% % dplySpotters = load(['/Users/olavobm/Documents/ROXSI_Postdoc' ...
% %                      '/MyResearch/ROXSI/Common_Code' ...
% %                      '/LargeScale_Data_2022/code_proc' ...
% %                      '/deploymentInfo_Spotters_ROXSI2022.mat']);
% % dplySpotters = dplySpotters.dployInfo_Spotters;


%% Define time limits for computing directional spectra
% (or comment this out if you just want to use all of
% the data available).

timelims_L2proc = [datetime(2022, 06, 20, 00, 00, 00), ...
                   datetime(2022, 07, 10, 00, 00, 00)];
timelims_L2proc.TimeZone = 'America/Los_Angeles';


%%

%
dt = 0.4;    % the sampling interval in seconds
nfft = 256;    % the number of points in each window (nfft)
dtheta = 1;    % the angle bin width for the spectral analysis
% analysis_period_hours = 0.5;    % the analysis period in hours.
analysis_period_hours = 1;


%%

% % dspec_method = ["EMEM","IMLM","MLM","MEM"]; % the directional algorithms to use
% dspec_method = ["EMEM","IMLM","MLM"]; % the directional algorithms to use

%
dspec_method = "EMEM";

% CHECK that the index for the displacement and location files are
% consistent (line 64)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The directional algorithms available are listed in dat2spec.m 
% 'BDM'  Bayesian Directional Spectrum Estimation Method 
% 'MLM'  Maximum Likelihood Method (default)
% 'IMLM' Iterative  Maximum Likelihood Method 
% 'MEM'  Maximum Entropy Method   (slow)
% 'EMEM' Extended Maximum Entropy Method


% For the reduced data structure (i.e. that does not
% have the full directional spectra). Use only 1
% method to compute the time-averaged directional
% spectrum
dspec_method_reduced = "EMEM";


%%

% % % index the displacement and depth files
% % % The location data must have been processed to include depth 
% % % % cd(displacement_data_path)
% % spotter_disp = dir(fullfile(displacement_data_path, "*displacement_deployed.mat"));
% % spotter_loc = dir(fullfile(displacement_data_path, "*location_deployed.mat"));


%%

%
analysis_period_seconds = analysis_period_hours * 3600; % hourly estimates are being generated
N = analysis_period_seconds/dt;
% t = 0:dt:(N-1)*dt;    % a time vector that will be passed to wafo 
t = [0:1:(N-1)].*dt;    % a time vector that will be passed to wafo 
t = t(:);

%
Nt = 360/dtheta + 1; % the number of anglular bins, dtheta makes more sense

%
T = nfft.*dt; % the window length in seconds
df = 1/T;    % frequency resolution, a consequence of T
fN = 0.5/dt; % Nyquist frequency
f = 0:df:fN;
f = f(:);

%
pos = [0,  0,  0; ...
       0,  0,  0; ...
       0,  0,  0; ...
       1, 16, 17; ...
       1,  0,  0].';    % data to the wafo, look at the
                        % dat2dspec.m and tran.m for more details 


%%
% -------------------------------------------------------------------------
% --------------------------- DO DATA PROCESSING --------------------------
% --------------------- (COMPUTE DIRECTIONAL SPECTRA) ---------------------
% -------------------------------------------------------------------------

%%

%
log_file_name = ['log_Spotter_procL2_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
diary(fullfile(dir_outlvl1, log_file_name))
%
totalRunTime = tic;


%%

disp(' '), disp(' ')
disp('------------------------------ L2 data processing for Spotters: ------------------------------')
for i = 1:length(list_Spotters)
    disp([num2str(i) ' - ' list_Spotters{i}])
end


%%

%
lprogress_switch = true;

%%

%
for i = 1:length(list_Spotters)
    
    %% --------- LOAD DISPLACEMENT FILE ---------

    % Print message to the screen
    disp(' ')
    disp(' ')
    %
    disp(['------------------ Load Spotter data from ' list_Spotters{i} ' ------------------'])
% %     %
% %     data_buoy = load(fullfile(dir_data_level_1, [list_Spotters{i} '.mat']));
% %     data_buoy = data_buoy.s;
% %     %
% %     data_buoy.displacement.time.TimeZone = 'America/Los_Angeles';
    %
    spotterL1 = load(fullfile(dir_data_level_1, [list_Spotters{i} '.mat']));
    spotterL1 = spotterL1.s;


    %% --------- GET LOCATION FILE ---------
    %
    for i2 = 1:length(spotter_location)
        %
        if strcmp(spotter_location(i2).dataID, list_Spotters{i})
            ind_match = i2;
        end
    end


    %% --------- GET ONLY THE NECESSARY VARIABLES ---------

    %
    vars2dirspectra.timeedges = spotter_location(ind_match).timetrimedges;

    %
    vars2dirspectra.spotterloc.dtime = spotter_location(ind_match).location.time(spotter_location(ind_match).ltrimedges);
    vars2dirspectra.spotterloc.latitude = spotter_location(ind_match).location.("latitude (decimal degrees)")(spotter_location(ind_match).ltrimedges);
    vars2dirspectra.spotterloc.longitude = spotter_location(ind_match).location.("longitude (decimal degrees)")(spotter_location(ind_match).ltrimedges);
    vars2dirspectra.spotterloc.depth = -spotter_location(ind_match).location.z_msl(spotter_location(ind_match).ltrimedges);

    %
    lintrim_edges = (spotterL1.displacement.time >= vars2dirspectra.timeedges(1)) & ...
                    (spotterL1.displacement.time <= vars2dirspectra.timeedges(2));
    % 
    vars2dirspectra.spotterdisp.dtime = spotterL1.displacement.time(lintrim_edges);
    vars2dirspectra.spotterdisp.x = spotterL1.displacement.("x (m)")(lintrim_edges);
    vars2dirspectra.spotterdisp.y = spotterL1.displacement.("y (m)")(lintrim_edges);
    vars2dirspectra.spotterdisp.z = spotterL1.displacement.("z (m)")(lintrim_edges);

    
    %% Check clock  of the data -- THESE SHOULD BE GRIDDED!

    %
% %     %
% %     figure
% %         plot(vars2dirspectra.spotterdisp.dtime(1:end-1), ...
% %              diff(vars2dirspectra.spotterdisp.dtime), '.-')
% %     %
% %     grid on
% %     set(gca, 'FontSize', 16)
% %     lims_dates = [datetime(2022, 06, 15), datetime(2022, 07, 23)];
% %     lims_dates.TimeZone = 'America/Los_Angeles';
% %     xlim(lims_dates)
% %     %
% %     title([list_Spotters{i}(1:3) ' SN:' list_Spotters{i}(end-3:end)], 'FontSize', 16)

    figure
        plot(vars2dirspectra.spotterloc.dtime(1:end-1), ...
             diff(vars2dirspectra.spotterloc.dtime), '.-')
        %
        grid on
        set(gca, 'FontSize', 16)
        lims_dates = [datetime(2022, 06, 15), datetime(2022, 07, 23)];
        lims_dates.TimeZone = 'America/Los_Angeles';
        xlim(lims_dates)
        %
        title([list_Spotters{i}(1:3) ' SN:' list_Spotters{i}(end-3:end)], 'FontSize', 16)


% %     keyboard
% %     continue
    

    %% Select the analysis time vector

    %
    if ~exist('timelims_L2proc', 'var')
        %
        dtime_proc_aux = spotterL1.spectra.dtime;
    else
        %
        dtime_proc_aux = timelims_L2proc(1) : hours(analysis_period_hours) : timelims_L2proc(2);
    end

    %
    analysis_periods = length(dtime_proc_aux);


% % % %     !!!TO BE DELETED!!!
% % 
% %     %% Select the analysis periods
% % 
% %     % pull out an analysis period of data
% %     % ind is the index of the start of the first full hour
% %     ind_start = find(minute(vars2dirspectra.spotterdisp.dtime) == 0, 1);
% % 
% %     dtime = vars2dirspectra.spotterdisp.dtime(ind_start) : ...
% %                     hours(analysis_period_hours) : ...
% %             (vars2dirspectra.spotterdisp.dtime(end) - hours(analysis_period_hours));
% % 
% %     analysis_periods = length(dtime); % this is the number of analysis period in the data


    %% Do the calculation

    % preallocate for speed
    S_f_theta_temp = NaN(length(f), Nt, analysis_periods, length(dspec_method));
    D_f_theta_temp = NaN(length(f), Nt, analysis_periods, length(dspec_method));
    S_f_temp = NaN(length(f), analysis_periods);
    depth = NaN(1, analysis_periods);
    lat = NaN(1, analysis_periods);
    lon = NaN(1, analysis_periods);

    % Print message to the screen
    disp(' '), disp(' ')
    %
    disp(['--------- Computing directional spectrum ' ...
          'for Spotter ' list_Spotters{i} ' ---------'])
    %
    disp(['The total number of analysis ' ...
          'periods is: ' num2str(analysis_periods)])
    
    % Loop over analysis periods
    tic
    for sample = 1:(analysis_periods)    % - 1 until I check what could be a minor bug/feature in the demo script
        
% %         % Old stuff
% %         data_index = ind_start + (sample -1) *N : ...
% %                      ind_start + sample *N -1;
        % 
        linanalysis_disp = (vars2dirspectra.spotterdisp.dtime >= (dtime_proc_aux(sample) - (hours(analysis_period_hours)/2))) & ...
                           (vars2dirspectra.spotterdisp.dtime  < (dtime_proc_aux(sample) + (hours(analysis_period_hours)/2)));
        % This is slower than the index approach, but shouldn't
        % make much of a difference because the directional 
        % spectra takes a lot more time
        %
        xt = vars2dirspectra.spotterdisp.x(linanalysis_disp);
        yt = vars2dirspectra.spotterdisp.y(linanalysis_disp);
        zt = vars2dirspectra.spotterdisp.z(linanalysis_disp);
        dtime_sample  = vars2dirspectra.spotterdisp.dtime(linanalysis_disp); 
     

        % the index of the depth record with the same dtime as the displacement 
        linanalysis_location = (vars2dirspectra.spotterloc.dtime >= (dtime_proc_aux(sample) - (hours(analysis_period_hours)/2))) & ...
                               (vars2dirspectra.spotterloc.dtime  < (dtime_proc_aux(sample) + (hours(analysis_period_hours)/2)));
        %
        lat(sample) = mean(vars2dirspectra.spotterloc.latitude(linanalysis_location), 'omitnan');
        lon(sample) = mean(vars2dirspectra.spotterloc.longitude(linanalysis_location), 'omitnan');
        depth(sample) = mean(vars2dirspectra.spotterloc.depth(linanalysis_location), 'omitnan');
        h = abs(depth(sample));
    
        %
        Data = [t, zt, xt, yt];
        
        % Loop over methods to compute directional spectra
        for jj = 1 : length(dspec_method)
            if jj == 1
                [Sd,D,Sw] = dat2dspec(Data, pos, h, nfft, Nt, dspec_method(jj));
            else
                [Sd,D] = dat2dspec(Data, pos, h, nfft, Nt, dspec_method(jj));
            end
            S_f_theta_temp(:, :, sample, jj) = Sd.S'.*(2*pi); %scaled because I want direction in degrees
            D_f_theta_temp(:, :, sample, jj) = D.S';
        end
        %       
        S_f_temp(:, sample) = Sw.S.*(2*pi);    
        
        % ------------------------------
        % Print progress message to the screen (the
        % second argument is the percentage step when
        % the progress will be printed to the screen)
        if mod(round(100*sample/analysis_periods), 25)==0
            %
            if lprogress_switch

                %
                disp(' '), disp(' ')
                disp(['----- Done with analysis period ' num2str(sample) ' ' ...
                      'out of ' num2str(analysis_periods) ' -----'])
                toc
                %
                lprogress_switch = false;
            end

        else
            lprogress_switch = true;
        end
        % ------------------------------
        
    end

    % from the wafo documentation
%             theta  = angle vector -pi..pi of length Nt 
%                      (theta = 0 -> + x-axis, theta = pi/2 -> + y-axis) 
% so in this angle definition, theta = 0 mean waves heading towards the
% east, (from the west) 

    % change the angular definition from cartesian to nautical
    dir = rad2deg(Sd.theta(1:end));
    dir_naut_unsorted= mod(270 - rad2deg(Sd.theta(1:end-1)),360);
    [dir_naut, ind] = sort(dir_naut_unsorted);


    %% Organize directional spectra output


    % ---------------------------------
    % First copy fields from L1 to L2 data strcture
    spotterL2.mooringID = spotterL1.mooringID;
    spotterL2.SN = spotterL1.SN;
    spotterL2.site = spotterL1.site;
    spotterL2.latitude = spotterL1.latitude;
    spotterL2.latitude = spotterL1.longitude;
    spotterL2.X = spotterL1.X;
    spotterL2.Y = spotterL1.Y;


    % ---------------------------------
% %     %
% %     spotterL2.latitude = lat(:);
% %     spotterL2.longitude = lon(:);
% %     %
% %     spotterL2.depth = depth(:);
% %     
% %     %
% %     spotterL2.dthour = analysis_period_hours;
% %     spotterL2.dtime = dtime(:);
    
    % ---------------------------------
    spotterL2.dtime =   dtime_proc_aux(:);

    % 
    spotterL2.frequency = f;
    
    % REMOVE 0'TH FREQUENCY???

    % The sea surface elevation spectra from WAFO
    % (L1 already has these spectra, but NOT computed
    % by WAFO)
    spotterL2.S_f = S_f_temp;

    %
    spotterL2.nfft = nfft;
    spotterL2.df = f(2) - f(1);


    % ---------------------------------
    % Directional spectra

    %
    spotterL2.direction_nautical = dir_naut;

    % The list of methods used in the calculation
    % of the directional spectra
    spotterL2.list_methods_dirspec = dspec_method;

    % Save each dspec method to the spotterL2 structure
    for jj = 1 : length(dspec_method)

        %
        spotterL2.(dspec_method(jj)).S_f_theta = S_f_theta_temp(:, ind, :, jj);
        spotterL2.(dspec_method(jj)).D_f_theta = D_f_theta_temp(:, ind, :, jj);

    end

    % Compute time-averaged direction spectrum



    
    %% Compute bulk parameters from wave spectrum

    %
    [freq_peak, freq_mean, Hsig] = bulkstats_from_wave_spectrum(spotterL2.frequency, spotterL2.S_f);
    %
    spotterL2.peak_f = freq_peak(:);
    spotterL2.mean_f = freq_mean(:);
    spotterL2.Hsig = Hsig(:);

    %
    for i2 = 1 : length(dspec_method)

        %
        bulkstats_aux = bulkstats_from_wave_dirspectrum(spotterL2.frequency, spotterL2.direction_nautical, ...
                                                        spotterL2.(dspec_method(i2)).D_f_theta, spotterL2.(dspec_method(i2)).S_f_theta, ...
                                                        spotterL2.S_f);
    
        %
        spotterL2.(dspec_method(i2)).mean_dir_f = bulkstats_aux.dir_mean_f;
        spotterL2.(dspec_method(i2)).mean_spread_f = bulkstats_aux.spread_mean_f;
        %
        spotterL2.(dspec_method(i2)).mean_dir = bulkstats_aux.dir_mean(:);
        spotterL2.(dspec_method(i2)).mean_spread = bulkstats_aux.spread_mean(:);
        %
        spotterL2.(dspec_method(i2)).peak_dir = bulkstats_aux.dir_peak(:);
        spotterL2.(dspec_method(i2)).peak_spread = bulkstats_aux.spread_peak(:);
    end


    %% Save directional spectra
    %
    disp(['--- Saving full directional spectra for ' list_Spotters{i} ' ---'])
    %
    fname = fullfile(dir_output_level_2, [list_Spotters{i} '_dspec.mat']);
    save(fname, "spotterL2" , '-v7.3')

    %
    disp(['--- Saving reduced (averaged) spectra for ' list_Spotters{i} ' ---'])
    %
    for i2 = 1:length(dspec_method)
        spotterL2.(dspec_method(i2)) = rmfield(spotterL2.(dspec_method(i2)), 'S_f_theta');
        spotterL2.(dspec_method(i2)) = rmfield(spotterL2.(dspec_method(i2)), 'D_f_theta');
    end
    %
    fname = fullfile(dir_output_level_2, [list_Spotters{i} '_dspec_reduced.mat']);
    save(fname, "spotterL2" , '-v7.3')


    %% Plot displacement spectrum (averaged in direcion, 
    % as a function of time and frequency)

    %
    hfig_spec = figure;
        pcolor(spotterL2.dtime, spotterL2.f, log10(spotterL2.S_f))
        shading flat

        %
        caxis([-4.5, 0.5])
        %
        hcb = colorbar;
            hcb.Label.Interpreter = 'Latex';
            hcb.Label.String = '$\log_{10}$ of displacement variance [m$^2$ Hz$^{-1}$]';
            hcb.Label.FontSize = 18;

    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'XGrid', 'on', 'YGrid', 'on', ...
             'YScale', 'log')
    set(gca, 'XLim', [datetime(2022, 06, 15, 06, 0, 0, 'TimeZone', 'America/Los_Angeles'), ...
                      datetime(2022, 07, 23, 12, 0, 0, 'TimeZone', 'America/Los_Angeles')], ...
             'YLim', [0, 1.25])
    %
    xlabel('Time [PDT]', 'Interpreter', 'Latex', 'FOntSize', 22)
    ylabel('Frequency [Hz]', 'Interpreter', 'Latex', 'FOntSize', 22)
    %
    title(['ROXSI 2022: displacement spectra. Spotter ' list_Spotters{i}(1:3) ...
           ' SN ' list_Spotters{i}(9:12)], 'Interpreter', 'Latex', 'FOntSize', 22)

    %
    set(gcf, 'units', 'normalized')
    set(gcf, 'Position', [0.35, 0.28, 0.41, 0.27])

    %
    disp(['--- Saving spectra figure for ' list_Spotters{i} ' ---'])
    %
    exportgraphics(hfig_spec, fullfile(dir_output_level_2, ['spectra_' list_Spotters{i} '.png']), 'Resolution', 300)
    %
    pause(5)
    %
    close(hfig_spec)

    %%
    %
    disp(['------------------ Done with directional spectrum ' ...
          'for Spotter ' list_Spotters{i} ' ------------------'])
end



%% End log file

%
disp('###################### Done with L2 data processing for all Spotters ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');



