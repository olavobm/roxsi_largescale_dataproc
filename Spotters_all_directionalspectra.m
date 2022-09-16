%% Compute directional Spectra for Spotters in ROXSI 2022
% 
%


clear
close all


%%
% -----------------------------------------------------
% ----------------- PRELIMINARY STUFF -----------------
% -----------------------------------------------------

%%

%
dir_data_level_1 = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';

 
% % All Spotters and Smart Moorings
% list_Spotters = {'B01_spot1150', 'B01_spot1158', 'B03_spot1152', ...
%               'B05_spot1153', 'X01_spot1151', 'X03_spot1157', ...
%               'X04_spot1155', ...
%               'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
%               'E07_spot1855', 'E07_spot1857', ...
%               'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
%               'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

% Just a few to test
% % list_Spotters = {'B03_spot1152', 'B05_spot1153', 'E02_spot1859'};
list_Spotters = {'B05_spot1153', 'E02_spot1859'};

% Output directory
dir_output_level_2 = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level2_Data/Spotter_Level2/';


% Print message to the screen:
disp(' '), disp(' ')
%
disp('Will compute directional spectrum for Spotters:')
list_Spotters
%
disp('Output directory is')
dir_output_level_2


%% Load the location file with bathymetry
% (actually, elevation to mean sea level)

spotter_location = load(fullfile(dir_data_level_1, 'Spotter_all_location_depth.mat'));
spotter_location = spotter_location.spotter_location;

% % %% Load trimming table for Spotters
% % %
% % % The trimming for Smart Moorings was defined in terms of
% % % pressure data, which had some issues for a few
% % % Smart Moorings. Looking at ther definition, they work
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


%% Add WAFO toolbox to Matlab path

%
% addpath(genpath(fullfile(repo_dirpath(), 'wafo')))
addpath(genpath(fullfile(repo_dirpath(), 'Spotter_DirectionalSpectra', 'wafo')))


%%
% ------------------------------------------------------
% --------- DEFINE PARAMETERS FOR DIRECTIONAL ----------
% ----------------- SPECTRA CALCULATION ----------------
% ------------------------------------------------------


%%
%
dt = 0.4; %the sampling interval in seconds
nfft = 256; %the number of points in each window (nfft)
dtheta = 1; %the angle bin width for the spectral analysis
analysis_period_hours = 0.5; %the analysis period in hours. 


%%
% % dspec_method = ["EMEM","IMLM","MLM","MEM"]; % the directional algorithms to use

% dspec_method = ["EMEM","IMLM","MLM"]; % the directional algorithms to use

dspec_method = "IMLM";

% CHECK that the index for the displacement and location files are
% consistent (line 64)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The directional algorithms available are listed in dat2spec.m 
% 'BDM'  Bayesian Directional Spectrum Estimation Method 
% 'MLM'  Maximum Likelihood Method (default)
% 'IMLM' Iterative  Maximum Likelihood Method 
% 'MEM'  Maximum Entropy Method   (slow)
% 'EMEM' Extended Maximum Entropy Method

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
t = 0:dt:(N-1)*dt;    % a time vector that will be passed to wafo 
t = t(:);

%
Nt = 360/dtheta + 1; % the number of anglular bins, dtheta makes more sense

%
T = nfft.*dt; % the window length in seconds
df = 1/T; % frequency resolution, a consequence of T
fN = 0.5/dt; % nyquist frequency
f = 0:df:fN; f = f(:);

%
pos = [0,  0,  0; ...
       0,  0,  0; ...
       0,  0,  0; ...
       1, 16, 17; ...
       1,  0,  0].';    % data to the wafo, look at the
                        % dat2dspec.m and tran.m for more details 


%%
% ------------------------------------------------------
% ------------ COMPUTE DIRECTIONAL SPECTRA -------------
% ------------------------------------------------------

%%

%
lprogress_switch = true;

%%

%
for i = 1:length(list_Spotters)
    
    %% --------- LOAD DISPLACEMENT FILE ---------

    %
    data_buoy = load(fullfile(dir_data_level_1, [list_Spotters{i} '.mat']));
    data_buoy = data_buoy.s;
    %
    data_buoy.displacement.time.TimeZone = 'America/Los_Angeles';


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
    lintrim_edges = (data_buoy.displacement.time >= vars2dirspectra.timeedges(1)) & ...
                    (data_buoy.displacement.time <= vars2dirspectra.timeedges(2));
    % 
    vars2dirspectra.spotterdisp.dtime = data_buoy.displacement.time(lintrim_edges);
    vars2dirspectra.spotterdisp.x = data_buoy.displacement.("x (m)")(lintrim_edges);
    vars2dirspectra.spotterdisp.y = data_buoy.displacement.("y (m)")(lintrim_edges);
    vars2dirspectra.spotterdisp.z = data_buoy.displacement.("z (m)")(lintrim_edges);

    
    %%
    
% %     load(spotter_disp(i).name);
% %     %%%%%%%%% This needs to correct for the displacement file name 
% %     site_name = spotter_disp(i).name(1:end-26);
% %     load(spotter_loc(i).name); % check that the index is the same for the displacement and location file... it should ne

    %%

    % pull out an analysis period of data
    % ind is the index of the start of the first full hour
    ind_start = find(minute(vars2dirspectra.spotterdisp.dtime) == 0, 1);

    dtime = vars2dirspectra.spotterdisp.dtime(ind_start) : ...
                    hours(analysis_period_hours) : ...
            (vars2dirspectra.spotterdisp.dtime(end) - hours(analysis_period_hours));

    analysis_periods = length(dtime); % this is the number of analysis period in the data


    %%

    % preallocate for speed
    S_f_theta_temp = NaN(length(f), Nt, analysis_periods, length(dspec_method));
    D_f_theta_temp = NaN(length(f), Nt, analysis_periods, length(dspec_method));
    S_f_temp = NaN(length(f), analysis_periods);
    depth = NaN(1, analysis_periods);
    lat = NaN(1, analysis_periods);
    lon = NaN(1, analysis_periods);

    % Print message to the screen
    disp(' ')
    disp(' ')
    %
    disp(['------------------ Computing directional spectrum ' ...
          'for Spotter ' list_Spotters{i} ' ------------------'])
    %
    disp(['----- The total number of analysis ' ...
          'periods is: ' num2str(analysis_periods) ' -----'])
    
    % Loop over analysis periods
    tic
    for sample = 1:(analysis_periods - 1)    % - 1 until I check what could be a minor bug/feature in the demo script
        
        %
        data_index = ind_start + (sample -1) *N : ...
                     ind_start + sample *N -1;
        xt = vars2dirspectra.spotterdisp.x(data_index);
        yt = vars2dirspectra.spotterdisp.y(data_index);
        zt = vars2dirspectra.spotterdisp.z(data_index);
        dtime_sample  = vars2dirspectra.spotterdisp.dtime(data_index); 
     

        % the index of the depth record with the same dtime as the displacement 
        good = vars2dirspectra.spotterloc.dtime >= dtime_sample(1) & ...
               vars2dirspectra.spotterloc.dtime <= dtime_sample(end);
        %
        lat(sample) = mean(vars2dirspectra.spotterloc.latitude(good),'omitnan');
        lon(sample) = mean(vars2dirspectra.spotterloc.longitude(good),'omitnan');
        depth(sample) = mean(vars2dirspectra.spotterloc.depth(good), 'omitnan');
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
        S_f_temp(:,sample) = Sw.S.*(2*pi);    
        
        % ------------------------------
        % Print progress message to the screen (the
        % second argument is the percentage for which
        % the progress will be printed to the screen)
        if mod(round(100*sample/analysis_periods), 5)==0
            %
            if lprogress_switch

                %
                disp(' ')
                disp(' ')
                %
                toc
                %
                disp(['----- Done with analysis period ' num2str(sample) ' ' ...
                      'out of ' num2str(analysis_periods) ' -----'])
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


    %% Now save directional spectra output

    %
    dspec.site = list_Spotters{i};

    % Save each dspec method to the dspec structure
    for jj = 1 : length(dspec_method)
        eval("dspec.S_f_theta_" + dspec_method(jj) + "= S_f_theta_temp(:,ind,:,jj);")
        eval("dspec.D_f_theta_" + dspec_method(jj) + "= D_f_theta_temp(:,ind,:,jj);")
    end
    %
    dspec.S_f = S_f_temp;
    dspec.f = f;
    dspec.direction_nautical = dir_naut;
    dspec.dtime = dtime(:);
    dspec.depth = depth(:);
    dspec.lon   = lon(:);
    dspec.lat   = lat(:);
    dspec.nfft = nfft;
    dspec.analysis_period_hours = analysis_period_hours;

    %
    disp(['--- Saving directional spectra ' list_Spotters{i} ' ---'])
    %
    fname = fullfile(dir_output_level_2, [list_Spotters{i} '_dspec.mat']);
    save(fname, "dspec" , '-v7.3')

    %
    disp(['------------------ Done with directional spectrum ' ...
          'for Spotter ' list_Spotters{i} ' ------------------'])
end




