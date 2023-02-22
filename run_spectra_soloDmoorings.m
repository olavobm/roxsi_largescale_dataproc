%% Computes surface elevation spectra from pressure for
% all SoloD's -- in a consistent manner with the directional
% spectra product.

clear
close all


%%
% ----------------------------------------
% --------------- PREAMBLE ---------------
% ----------------------------------------

%%

% % addpath(genpath('/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc/'))


%% Directory with data files of data that
% has been prepared for this script

% % % 
% % dir_data = ['/project/CSIDE/ROXSI/'];
%
dir_output = ['/home/omarques/Documents/obm_ROXSI/Analyses' ...
              '/dirspectra_product/new_results_3/data_with_dirspectra/'];


%% Moorings to use for calculation

%
list_moorings = {'A02', 'A04', 'A05', 'A06', 'A07', 'A09', ...
                 'B04', 'B06', 'B07', 'B09', 'B12', 'B14', 'B16', 'B18', ...
                 'C02', 'C03', 'C05', 'C06', 'C07', 'C08', 'C09', ...
                 'D01', 'D02', ...
                 'X07', 'X08', 'X09', 'X10'};


%% Define time grid for computing -- get other L2 data
% or make sure that it is the same time grid as the Spotter

%
dt_spec = hours(1);

% % % Full timeseries (excluding remaining data after
% % % recovery started, and a bit on the first day)
% % dtimegridlims = [datetime(2022, 06, 15, 20, 00, 00), ...
% %                  datetime(2022, 07, 20, 05, 00, 00)];
% %              
             
% A sample
dtimegridlims = [datetime(2022, 06, 28, 12, 00, 00), ...
                 datetime(2022, 06, 30, 12, 00, 00)];
             

%
dtimegrid = dtimegridlims(1) : dt_spec : dtimegridlims(2);
dtimegrid.TimeZone = 'America/Los_Angeles';


%% Define FFT properties -- time length of
% individual FFTs, averaging time window,
% and high-frequency cut-off

%
dt_fftsegment = 120;    % in seconds
% % dt_fftsegment = 240;    % in seconds
dt_avgfft = seconds(dt_spec);    % in seconds

% The cut-off frequency (in Hz) to use because
% the transfer function blows up 
% freq_cutoff = 0.15;
% freq_cutoff = 0.2;
freq_cutoff = 0.3;


%% Define frequency bands for bulk calculations

%
bulkbands(1).ID = 'swell';
bulkbands(1).freqlims = [0.06, 0.1];

%
bulkbands(2).ID = 'sea';
bulkbands(2).freqlims = [0.10001, 0.2];

%
bulkbands(3).ID = 'seaswell';
bulkbands(3).freqlims = [0.06, 0.2];


%%
% ---------------------------------------------
% ------ DO CALCULATIONS FOR EACH SOLOD -------
% ---------------------------------------------


%% Print progress message

%
disp(' '), disp(' ')
%
disp('----- Recomputing surface elevation spectra from pressure data. Doing calculations for moorings:  -----')
for i = 1:length(list_moorings)
    disp([num2str(i) ') ' list_moorings{i}])
end

%
totalRunTime = tic;


%%

% Loop over SoloDs
for i1 = 1:length(list_moorings)


    %% --------------------------------------
    % Print progress message
    %
    disp(' ')
    %
    disp(['--- Starting computations for mooring ' list_moorings{i1} ' ---'])

    %% --------------------------------------
    % Load data

    %
    dtimelims_load = [(dtimegridlims(1) - dt_spec), (dtimegridlims(2) + dt_spec)];

    %
    dataL1 = ROXSI_load_pressureL1(list_moorings{i1}, dtimelims_load);

    %% --------------------------------------
    % Compute bottom depth from hydrostatic pressure

    %
    dataL1.bottomdepthfrompres = dataL1.pressure + dataL1.zhab;


    %% --------------------------------------
    % Get sampling rate and set number of points for FFT
    
    %
    dtsampling = seconds(diff(dataL1.dtime(1:2)));
    nfft = (1/dtsampling) * dt_fftsegment;


    %% --------------------------------------
    % Add fields to output data structure

    %
    solodL2.mooringID = list_moorings{i1};
    solodL2.instrument = "SoloD";
    solodL2.latitude = dataL1.latitude;
    solodL2.longitude = dataL1.longitude;
    solodL2.site = dataL1.site;
    solodL2.X = dataL1.X;
    solodL2.Y = dataL1.Y;

    %
    solodL2.fftwindow = dt_fftsegment;
    solodL2.fftavgwindow = dt_avgfft;
    %
    solodL2.freqcutoff = freq_cutoff;

    %
    solodL2.dt = dt_spec;
    solodL2.dtime = dtimegrid(:);


    %% --------------------------------------
    % Use pressure to compute frequency spectra
    % of elevation.

    %
    dtime_lims = dtimegrid([1, end]);

    %
    disp('--- Computing elevation spectra ---')
    %
    [See_aux, ...
     frequency_aux, dtimespec_aux, ...
     k_aux, bottomdepth_aux] = ...
                    presdata_to_See(dataL1.dtime, dataL1.pressure, ...
                                    dataL1.zhab, ...
                                    dataL1.bottomdepthfrompres, ...
                                    solodL2.fftwindow, solodL2.fftavgwindow, ...
                                    dtimegridlims);
    disp('--- Done computing elevation spectra from pressure ---')

    %
    lbelowcutoff = (frequency_aux <= freq_cutoff);

    %
    solodL2.frequency = frequency_aux(lbelowcutoff);
    solodL2.frequency = solodL2.frequency(:);

    %
    solodL2.See = See_aux(:, lbelowcutoff);

    %
    solodL2.bottomdepth = bottomdepth_aux(:);
    solodL2.k = k_aux(:, lbelowcutoff);
    
    
%     %% Get the frequency limits for each band for bulk statistics
%     
%     for i2 = 1:length(bulkbands)
%         %
%         bulkbands(i2).linlims = (solodL2.frequency >= bulkbands(i2).freqlims(1)) & ...
%                                 (solodL2.frequency <= bulkbands(i2).freqlims(2));
%     end


    %% --------------------------------------
    % Add kH, cp, and cg do data structure
    
	%
    H_array_aux = repmat(solodL2.bottomdepth, 1, length(solodL2.frequency));
    
	%
    solodL2.kH = solodL2.k .* H_array_aux;
    
	%
    solodL2.cp = wave_cp(solodL2.k, H_array_aux);
    solodL2.cg = wave_cg(solodL2.k, H_array_aux);
    
    
    %% --------------------------------------
    % Compute energy Flux without taking into account direction
    solodL2.Flux = (1025 * 9.8) * solodL2.See .* solodL2.cg;


    %% --------------------------------------
    % Compute mean frequency, peak frequency, and integrate
    % flux at different frequency bands

    %
    for i2 = 1:length(bulkbands)
        %
        [meanfreq_aux, peakfreq_aux] = wave_bulk_frequency(solodL2.frequency, ...
                                                           solodL2.See, ...
                                                           bulkbands(i2).freqlims, 2);
        %
        solodL2.(bulkbands(i2).ID).meanfreq = meanfreq_aux;
        solodL2.(bulkbands(i2).ID).peakfreq = peakfreq_aux(:);

        %
        linlims_aux = (solodL2.frequency >= bulkbands(i2).freqlims) & ...
                      (solodL2.frequency <= bulkbands(i2).freqlims);

        %
        solodL2.(bulkbands(i2).ID).Flux = trapz(solodL2.frequency(linlims_aux), ...
                                                solodL2.Flux(:, linlims_aux), 2);
    end


    %% --------------------------------------
    % Save data

    %
    disp(['--- Saving L2 SoloD data from mooring ' list_moorings{i1} ' ---'])
    
    %
    save(fullfile(dir_output, ['roxsi_solodL2_' list_moorings{i1} '.mat']), 'solodL2', '-v7.3')

    %%

    clear dataL1 solodL2


end

