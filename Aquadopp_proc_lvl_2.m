%% Do L2 processing for Aquadopp data
%
% L2 consists of:
%   - hourly depth-averaged currents
%   - depth-averaged currents at 1 Hz (if it is 1 Hz sampling)
%   - Spectra of u, v, w averaged over 1 hours (if it is 1 Hz sampling)
%       
%
%
%
%


clear
close all


%%
% --------------------------------------
% -------------- PREAMBLE --------------
% --------------------------------------

%% Directory of L1 data

%
% % dirparent_data = data_dirpath();
%
dirparent_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/';
dir_dataL1 = fullfile(dirparent_data, 'Level1_Data', 'Aquadopp_Level1');


%% Directory to output L2 data and figures

%
dir_output_parent = dirparent_data;
%
dir_output_data_L2 = fullfile(dir_output_parent, 'Level2_Data', 'Aquadopp_Level2');
dir_output_figs_L2 = fullfile(dir_output_data_L2, 'figs_QC');


% % %% Load ADCP deployment information
% % 
% % % Just to be clear: file and variable have
% % % the same name (though not a requirement)
% % % dir_coderepo = repo_dirpath();
% % dir_coderepo = '/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc';
% % %
% % load(fullfile(dir_coderepo, 'deploymentInfo_ROXSI2022.mat'), 'deploymentInfo_ROXSI2022')


%% List of Aquadopps that will be processed

% All Aquadopps
list_Aquadopp = {'A03_5380', ...
                 'B02_12507', 'B04_2147', 'B07_2141', 'B08_13288', 'B11_12280', ...
                 'C03_0709', ...
                 'D01_12346', 'D02_0653', ...
                 'E03_13300', 'E04_13172', 'E06_9736', 'E12_11150', ...
                 'F01_9995', 'F02_5838', 'F03_5384', 'F04_5401', 'F05_14032', ...
                 'X06_13290', 'X13_9945'};

% % A subset with a few
% list_Aquadopp = {'A03_5380', ...
%                  'B02_12507', ...
%                  'C03_0709', ...
%                  'E03_13300', ...
%                  'F01_9995', 'F02_5838', 'F03_5384'};

% A subset with a few
list_Aquadopp = {'B02_12507', 'E03_13300'};

%
Naquadopps = length(list_Aquadopp);


%%
% The Aquadopps sampling at 1 Hz were:
% 
%   - B08_13288, B11_12280
%   - E03_13300, E04_13172, E06_9736, E12_11150
%   - X06_13290, X13_9945
%   - F01_9995, F02_5838, F03_5384, F04_5401, F05_14032


%%
% --------------------------------------
% ------- DEFINE PARAMETERS AND --------
% ---- VARIABLES FOR DATA PROCESSING ---
% --------------------------------------

%%

% Time window for statistical averaging (in seconds)
windowavg = 60*60;

% Time window for individual fft (in seconds)
windowfft = 6*60;
    
% Sea-swell frequency band
freq_lims = [0.045, 0.3];


%%

list_copyvars = {'SN', 'mooringID', 'latitude', 'longitude', ...
                 'site', 'X', 'Y', ...
                 'samplingtime', 'averagingtime', 'binsize', 'zhab'};


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------------- DO DATA PROCESSING ----------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------



%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_Aquadopp_procL2_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_output_data_L2, log_file_name))
%
totalRunTime = tic;


%% Display on the screen:

%
disp(' '), disp(' ')
disp('------------------------------ Processing L2 data from Aquadopps ------------------------------')
disp('List of Aquadopps being processed:')
%
for i = 1:Naquadopps
    disp([num2str(i) ' - ' list_Aquadopp{i}])
end


%% Do L2 data processing

% Loop over Aquadopps in the list
for i1 = 1:Naquadopps

    % ----------------------------------------------------
    %
    disp(' '), disp(' ')
    disp(['--- Start L2 data processing for Aquadopp : ' list_Aquadopp{i1} ' ---'])

    % ----------------------------------------------------
    %
    disp('--- Loading L1 data ---')

    %
    aquadoppL1 = load(fullfile(dir_dataL1, ['roxsi_aquadopp_L1_' list_Aquadopp{i1}(1:3) '_' list_Aquadopp{i1}(5:end) '.mat']));
    aquadoppL1 = aquadoppL1.aquadoppL1;

    % Get rid of the substructure with smoothed/low-passed velocity
    aquadoppL1 = rmfield(aquadoppL1, 'averaged');

    % ----------------------------------------------------
    % Start populating L2 data structure

    %
    for i2 = 1:length(list_copyvars)
        aquadoppL2.(list_copyvars{i2}) = aquadoppL1.(list_copyvars{i2});
    end

    % ----------------------------------------------------
    % Compute hourly mean velocities
    %
    disp(['--- Computing averaged velocities over ' num2str(windowavg/60, '%.1f') ' min ---'])

    % Get first and last time stamps of hourly statistics where there is
    % sufficient data at the edges on either half of the averaging windo
    %
    time_edgedata_1 = aquadoppL1.dtime(1) + (hours(windowavg/3600)/2);
    time_edgedata_2 = aquadoppL1.dtime(end) - (hours(windowavg/3600)/2);

    %
    time_edge_1 = time_edgedata_1;
    time_edge_1.Minute = 0;
    time_edge_1.Second = 0;
    time_edge_1.Hour = time_edge_1.Hour + 1;
    %
    time_edge_2 = time_edgedata_2;
    time_edge_2.Minute = 0;
    time_edge_2.Second = 0;
    
    %
    timestatslims = [time_edge_1, time_edge_2];
    
    %
    aquadoppL2.dtime = timestatslims(1) : hours(windowavg/3600) : timestatslims(2);


    % Compute hourly mean velocities
    aquadoppL2.u = NaN(length(aquadoppL2.zhab), length(aquadoppL2.dtime));
    aquadoppL2.v = aquadoppL2.u;
    aquadoppL2.w = aquadoppL2.u;

    tic
    % Loop over ADCP bins
    for i2 = 1:length(aquadoppL2.zhab)
        %
        [aquadoppL2.u(i2, :), time_aux] = time_smooth_reg(aquadoppL1.dtime, aquadoppL1.u(i2, :), windowavg, timestatslims);
        aquadoppL2.v(i2, :) = time_smooth_reg(aquadoppL1.dtime, aquadoppL1.v(i2, :), windowavg, timestatslims);
        aquadoppL2.w(i2, :) = time_smooth_reg(aquadoppL1.dtime, aquadoppL1.w(i2, :), windowavg, timestatslims);
    end
    toc
   
    % Compute depth-averaged hourly mean velocity
    aquadoppL2.udepthavg = mean(aquadoppL2.u, 1, 'omitnan');
    aquadoppL2.vdepthavg = mean(aquadoppL2.v, 1, 'omitnan');
    aquadoppL2.wdepthavg = mean(aquadoppL2.w, 1, 'omitnan');
    
    keyboard

    %%

% %     % ----------------------------------------------------
% %     % Compute spectra for u, v, and w if sampling is 1 Hz
% %     if aquadoppL1.samplingtime == 1
% % 
% %         %
% %         disp(['--- Sampling rate of Aquadopp ' list_Aquadopp{i1} ' is ' num2str(aquadoppL1.samplingtime) ' second(s). Computing velocity spectra ---'])
% % 
% %         tic
% %         % Loop over ADCP bins
% %         for i2 = 1:length(aquadoppL1.zhab)
% % % %             %
% % % %             [Spp, timespec, freqvec, dof, avgpres] = ...
% % % %                         spectra_scalar_reg(spotsmartL1.dtime, spotsmartL1.pressure, ...
% % % %                                            windowfft, windowavg, timespeclims);
% %         end
% %         toc
% %     end
% % 
% %     % ----------------------------------------------------
% %     % Compute variance in the sea-swell bands if sampling is 1 Hz
% %     if aquadoppL1.samplingtime == 1
% %         %
% %         disp('--- Integrating spectra to compute variance for each velocity component in the sea-swell band ---')
% % 
% %         %
% %         aquadoppL2.spectra.freqband = freq_lims;
% % 
% %         %
% %         linfreqband_seaswell = (aquadoppL2.frequency >= aquadoppL2.freqband(1)) & ...
% %                                (aquadoppL2.frequency  < aquadoppL2.freqband(2));
% % 
% %         %
% %         aquadoppL2.spectra.uvar = trapz(aquadoppL2.frequency(linfreqband_seaswell), ...
% %                                         aquadoppL2.Suu(linfreqband_seaswell, :), 1);
% %         %
% %         aquadoppL2.spectra.vvar = trapz(aquadoppL2.frequency(linfreqband_seaswell), ...
% %                                         aquadoppL2.Svv(linfreqband_seaswell, :), 1);
% %         %
% %         aquadoppL2.spectra.wvar = trapz(aquadoppL2.frequency(linfreqband_seaswell), ...
% %                                         aquadoppL2.Sww(linfreqband_seaswell, :), 1);
% % 
% %     end
% % 
% %     %%
% %     % ----------------------------------------------------
% %     % Organize L2 data structure
% % 
% %     aquadoppL2.README = ['Level 2 Aquadopp data from ROXSI 2022'];
% % 
% %     % ----------------------------------------------------
% %     % Make QC plot
% % 
% %     % ----------------------------------------------------
% %     % Save L2 data structure
% % 
% %     %
% %     disp('----- Saving level 2 data -----')
% %     %
% %     str_filename = ['roxsi_aquadopp_L2_' char(spotsmartL2.mooringID) '_' char(spotsmartL2.SN)];
% %     save(fullfile(dir_output_data_L2, [str_filename '.mat']), 'aquadoppL2', '-v7.3')
% %     %
% %     disp('----- Done saving data -----')


% % %     %
% % %     disp('----- Saving level 2 QC plots -----')
% % %     %
% % %     str_filename = ['**typeofplot**_' char(spotsmartL2.mooringID) '_' char(spotsmartL2.SN)];
% % %     % Save figure as *.png
% % %     exportgraphics(fig_QC, fullfile(dir_output_data_L2, [str_filename '.png']), 'Resolution', 300)


    %% Print progress message

    %
    disp(' ')
    disp(['----- Done with L2 data processing for Aquadopp ' list_Aquadopp{i1} ' -----'])

    %
    toc(totalRunTime)


end


%% Print message the processing is done and close log file

%
disp('###################### Done with L2 data processing for all Aquadopps in the list ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');