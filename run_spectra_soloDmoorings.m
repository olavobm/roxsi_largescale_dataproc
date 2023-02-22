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
% % dir_dataaux = ['/Users/olavobm/Documents/ROXSI_Postdoc' ...
% %                '/MyResearch/figures_bydate/2023_02_01/'];
% % %
% % dir_outputdata = ['/Users/olavobm/Documents/ROXSI_Postdoc' ...
% %                   '/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022' ...
% %                   '/code_proc/directional_properties/'];

% % % 
% % dir_dataaux = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results/data_for_dirspectra/';
% % %
% % % % dir_outputdata = dir_dataaux;
% % dir_outputdata = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results/data_with_dirspectra/';


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
    dataSoloD = load(fullfile(dir_dataaux, ['adcpdata_' list_moorings{i1} '.mat']));
    dataSoloD = dataSoloD.dataADCP;

    %
    adcpL2.mooringID = list_moorings{i1};

end


%%
return







    
%%

% Loop over ADCPs
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
    dataSoloD = load(fullfile(dir_dataaux, ['adcpdata_' list_moorings{i1} '.mat']));
    dataSoloD = dataSoloD.dataADCP;

    %
    adcpL2.mooringID = list_moorings{i1};

    %% --------------------------------------
    % Set FFT parameters based on the sampling
    % rate of the i1'th instrument
    
    %
    dtsampling = seconds(diff(dataSoloD.dtime(1:2)));
    nfft = (1/dtsampling) * dt_fftsegment;

    %% --------------------------------------
    % Add fields to output data structure

    %
    adcpL2.instrument = dataSoloD.instrument;
    adcpL2.latitude = dataSoloD.latitude;
    adcpL2.longitude = dataSoloD.longitude;
    adcpL2.site = dataSoloD.site;
    adcpL2.X = dataSoloD.X;
    adcpL2.Y = dataSoloD.Y;

    %
    adcpL2.fftwindow = dt_fftsegment;
    adcpL2.fftavgwindow = dt_avgfft;
    %
    adcpL2.freqcutoff = freq_cutoff;

    %
    adcpL2.dt = dt_spec;
    adcpL2.dtime = dtimegrid(:);


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
                    presdata_to_See(dataSoloD.dtime, dataSoloD.pressure, ...
                                    dataSoloD.transducerHAB, ...
                                    dataSoloD.bottomdepthfrompres, ...
                                    adcpL2.fftwindow, adcpL2.fftavgwindow, ...
                                    dtime_lims);
    disp('--- Done computing elevation spectra from pressure ---')

    %
    lbelowcutoff = (frequency_aux <= freq_cutoff);

    %
    adcpL2.frequency = frequency_aux(lbelowcutoff);
    adcpL2.frequency = adcpL2.frequency(:);

    %
    adcpL2.See = See_aux(:, lbelowcutoff);

    %
    adcpL2.bottomdepth = bottomdepth_aux;
    adcpL2.k = k_aux(:, lbelowcutoff);
    
    
    %% Get the frequency limits for each band for bulk statistics
    
    for i2 = 1:length(bulkbands)
        %
        bulkbands(i2).linlims = (adcpL2.frequency >= bulkbands(i2).freqlims(1)) & ...
                                (adcpL2.frequency <= bulkbands(i2).freqlims(2));
    end

end








%%

return


%%

%
windows_spectra = [(6*60), (3600)];

% % %
% % freq_lims_spectra = [0.01, 0.325];
% % freq_lims_energetics = [0.05, 0.25];

%
% % freq_lims_spectra = [0.01, 1];
% % freq_lims_energetics = [0.05, 1.25];

%
freq_lims_spectra = [0.01, 0.5];
freq_lims_energetics = [0.05, 0.3];


%% Mooring on the B line

% % % Test
% % list_moorings = {'B06', 'B08'};
% % %list_moorings = {'B06', 'B08', 'B03', 'E08'};
% % 
% % 
% % %
% % list_moorings = {'B03', 'B04', 'B05', 'B06', ...
% %                  'E08', ...
% %                  'B08', 'B09', 'B10', 'B11', 'B12', ...
% %                  'B13', 'B14', 'B15', 'B16'};
          
% % %
% % list_moorings = {'B03', 'B04', 'B05', 'B06', ...
% %                  'E08', ...
% %                  'B08', 'B09', 'B10', 'B11', 'B12', ...
% %                  'B13', 'B14', 'B15', 'B16', ...
% %                  'A01', 'A02', 'A04', 'A05', 'A06', 'A07', ...
% %                  'C01', 'C02', 'C03', 'C05', 'C06', 'C07', 'C08', ...
% %                  'X01', 'X03', 'X04', 'X05', 'X06', 'X07', 'X08', 'X09', 'X10', 'X11', 'X12', 'X13'};

%
% list_moorings = {'B17', 'B18'};

% % list_moorings = {'B03'};
             
%             
list_moorings = {'B03', 'B04', 'B05', 'B06', ...
                 'E08', ...
                 'B08', 'B09', 'B10', 'B11', 'B12', ...
                 'B13', 'B14', 'B15'};

%              %             
% list_moorings = {'B03', 'B06', 'B10'};


%%

list_smartmoorings = {'E08'};


%%

% % time_lims_spectra = [datetime(2022, 06, 15, 12, 0, 0), ...
% %                      datetime(2022, 07, 20, 04, 0, 0)];

time_lims_spectra = [datetime(2022, 06, 18, 0, 0, 0), ...
                     datetime(2022, 07, 20, 04, 0, 0)];
                 
%%
% ----------------------------------------
% ----------- RECOMPUTE SPECTRA ----------
% ----------------------------------------

%%

%
totalRunTime = tic;


%%

%
for i = 1:length(list_moorings)
    
    %
    disp(['--- Recomputing spectra for ' list_moorings{i} ' ---'])
    
    %
    if ~any(ismember(list_smartmoorings, list_moorings{i}))
        %
        spectra_output = wavestats_recomputingspectra_roxsi2022(list_moorings{i}, windows_spectra, freq_lims_spectra, time_lims_spectra);
% %         spectra_output = rmfield(spectra_output, 'data');
        %
        spectraMooring.(list_moorings{i}) = spectra_output;
    else
        
        % Using Smart Mooring buoy
        spectra_output_buoy = wavestats_recomputingspectra_roxsi2022(list_moorings{i}, windows_spectra, freq_lims_spectra, time_lims_spectra);

        % Using Smart Mooring pressure
        spectra_output_pressure = wavestats_recomputingspectra_roxsi2022(list_moorings{i}, windows_spectra, freq_lims_spectra, time_lims_spectra, true);
        
        %
%         spectra_output_buoy = rmfield(spectra_output_buoy, 'data');
%         spectra_output_pressure = rmfield(spectra_output_pressure, 'data');
        
        %
        spectraMooring.(list_moorings{i}) = spectra_output_buoy;
        spectraMooring.([list_moorings{i} 'pres']) = spectra_output_pressure;
        
    end
  
    
    %
    disp(['--- Done with ' list_moorings{i} ' ---'])
    toc(totalRunTime)
end


%%
disp('--- DONE WITH RECOMPUTING SPECTRA FROM ALL MOORINGS ---')

%%
% ----------------------------------------
% ------- ORGANIZE SPECTRA IN NEW --------
% ------------ DATA STRUCTURE ------------
% ----------------------------------------

%%

list_fields = fieldnames(spectraMooring);

%%


%
mooringFlux.Nmoorings = length(list_fields);

%
mooringFlux.mooringID = list_fields;
mooringFlux.mooringID = mooringFlux.mooringID(:).';


%%

%
mooringFlux.instrument = strings(1, mooringFlux.Nmoorings);

%
for i = 1:mooringFlux.Nmoorings
    mooringFlux.instrument(i) = convertCharsToStrings(spectraMooring.(list_fields{i}).instrument);
end


%% Get mooring locations

% Load mooring table
mooringtable = load(['/home/omarques/Documents/MATLAB' ...
                     '/roxsi_largescale_dataproc/ROXSI2022_mooringtable.mat']);
mooringtable = mooringtable.mooringtable;

%
mooringFlux.latitude = NaN(mooringFlux.Nmoorings, 1);
mooringFlux.longitude = mooringFlux.latitude;

%
for i = 1:length(mooringFlux.mooringID)
    
    %
    lloop_aux = true;
    indmatch_aux = 0;
    %
    while (indmatch_aux < size(mooringtable, 1)) && lloop_aux
        %
        indmatch_aux = indmatch_aux + 1;
        %
        mooringID_ontable_aux = char(mooringtable.mooringID(indmatch_aux));
        
        %
        if strcmp(mooringFlux.mooringID{i}(1:3), mooringID_ontable_aux(1:3))
            %
            if ~strcmp(mooringtable.instrument(indmatch_aux), "tchain")
                lloop_aux = false;
                
            end
        end
    end
    
    %
    mooringFlux.latitude(i) = mooringtable.latitude(indmatch_aux);
    mooringFlux.longitude(i) = mooringtable.longitude(indmatch_aux);
    
end

%
[mooringFlux.X, mooringFlux.Y] = ROXSI_lltoxy(mooringFlux.latitude, mooringFlux.longitude);


%%

%
linlims_aux = spectraMooring.(list_fields{1}).spectra.frequency <= freq_lims_spectra(2);

%
mooringFlux.dtime = spectraMooring.(list_fields{1}).spectra.dtime;
mooringFlux.frequency = spectraMooring.(list_fields{1}).spectra.frequency(linlims_aux);
mooringFlux.df = spectraMooring.(list_fields{1}).spectra.df;



%% Add bottom depth and elevation spectra

%
mooringFlux.bottomdepth = NaN(length(list_fields), length(mooringFlux.dtime));

%
prealloc_matrix_aux = NaN(mooringFlux.Nmoorings, length(mooringFlux.dtime), length(mooringFlux.frequency));
mooringFlux.See = prealloc_matrix_aux;

%
for i = 1:mooringFlux.Nmoorings
    %
    mooringFlux.bottomdepth(i, :) = spectraMooring.(list_fields{i}).spectra.bottomdepth(:)';
    
    linlims_aux = spectraMooring.(list_fields{i}).spectra.frequency <= freq_lims_spectra(2);
    
    %
    mooringFlux.See(i, :, :) = spectraMooring.(list_fields{i}).spectra.See(linlims_aux, :).';
    
end


%%

% %
% mooringFlux.freqlims_Flux = freq_lims_energetics;

%
prealloc_matrix_aux = NaN(mooringFlux.Nmoorings, length(mooringFlux.dtime), length(mooringFlux.frequency));

%
mooringFlux.k = prealloc_matrix_aux;
mooringFlux.kH = prealloc_matrix_aux;
mooringFlux.cp = prealloc_matrix_aux;
mooringFlux.cg = prealloc_matrix_aux;
mooringFlux.E = prealloc_matrix_aux;
mooringFlux.Flux = prealloc_matrix_aux;


%%
% -------------------------------------
% --- COMPUTE ENERGETICS SPECTRALLY ---
% -------------------------------------

%%

%
rho0 = 1030;
g = 9.8;

% (PS: frequency in Hz and k in radians per meter)
disp_rel = @(k, freq, H) g*k*tanh(k*H) - (2*pi*freq)^2;


%% Compute energetics

%
disp('--- Starting to compute wavenumber ---')
toc(totalRunTime)

%
for i1 = 1:mooringFlux.Nmoorings
    
% %     %
% %     waveStats.(list_mooring_IDs{i1}).spectrally.k = NaN(length(waveStats.(list_mooring_IDs{i1}).frequency), ...
% %                                                         length(waveStats.(list_mooring_IDs{i1}).dtime));
% %     waveStats.(list_mooring_IDs{i1}).spectrally.cg = waveStats.(list_mooring_IDs{i1}).spectrally.k;
% %     waveStats.(list_mooring_IDs{i1}).spectrally.E = waveStats.(list_mooring_IDs{i1}).spectrally.cg;
% %     waveStats.(list_mooring_IDs{i1}).spectrally.Flux = waveStats.(list_mooring_IDs{i1}).spectrally.cg;
% % % 
% %     % Identify type of instrument
% %     ind_underscore_aux = strfind(list_alldata{i1}, '_');
% %     %
% %     strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));
% % 
% %     % Make bottom depth a column vector and positive
% %     dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth) = ...
% %             abs(dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth)(:));
    

    %
    bottomdepth_aux = mooringFlux.bottomdepth(i1, :);
    
    %
    for i2 = 1:length(mooringFlux.dtime)
        %
        if (bottomdepth_aux(i2) > 0.2)
            %
            for i3 = 1:length(mooringFlux.frequency)
                %
                if (mooringFlux.frequency(i3) >= freq_lims_energetics(1)) && ...
                   (mooringFlux.frequency(i3) < freq_lims_energetics(2))    
                    
                    %
                    disp_rel_eval = @(k) disp_rel(k, mooringFlux.frequency(i3), bottomdepth_aux(i2));
                    %
                    k_aux = fzero(disp_rel_eval, [(2*pi/(5000)), (2*pi/(1))]);

                    %
                    mooringFlux.k(i1, i2, i3) = k_aux;
        
                end
            end
        end
    end
    
    
    % Compute group velocity
    mooringFlux.kH = mooringFlux.k .* ...
                     repmat(mooringFlux.bottomdepth, 1, 1, length(mooringFlux.frequency));

    %
    mooringFlux.cp = sqrt(g * tanh(mooringFlux.kH) ./ mooringFlux.k);
    
    %
    mooringFlux.cg = 0.5 * mooringFlux.cp .* ...
                     (1 + (2*mooringFlux.kH./sinh(2*mooringFlux.kH)));
                               
	% Compute Energy
    mooringFlux.E = rho0 * g * mooringFlux.df * mooringFlux.See;

    % Compute Energy Flux per frequency band
    mooringFlux.Flux = mooringFlux.cg .* mooringFlux.E;
    
% %     %
% %     linlims_aux = waveStats.(list_mooring_IDs{i1}).linlims;
% %     %
% %     waveStats.(list_mooring_IDs{i1}).spectrally.inband.freq_lims = waveStats.(list_mooring_IDs{i1}).freq_lims;
% %     waveStats.(list_mooring_IDs{i1}).spectrally.inband.linlims = linlims_aux;
% %     
% %     %
% %     waveStats.(list_mooring_IDs{i1}).spectrally.inband.Hsig = trapz(waveStats.(list_mooring_IDs{i1}).frequency(linlims_aux), ...
% %                                                                     waveStats.(list_mooring_IDs{i1}).See(linlims_aux, :), 1);
% % 	waveStats.(list_mooring_IDs{i1}).spectrally.inband.Hsig = 4*sqrt(waveStats.(list_mooring_IDs{i1}).spectrally.inband.Hsig);
% %     
% %     % Compute Energy Flux in one frequency band
% %     waveStats.(list_mooring_IDs{i1}).spectrally.inband.Flux = trapz(waveStats.(list_mooring_IDs{i1}).frequency(linlims_aux), ...
% %                                                                     waveStats.(list_mooring_IDs{i1}).spectrally.E(linlims_aux, :) .* ...
% %                                                                     waveStats.(list_mooring_IDs{i1}).spectrally.cg(linlims_aux, :), 1);
% %     
% % 	%
% %     waveStats.(list_mooring_IDs{i1}).spectrally.cumulativeFlux = NaN(size(waveStats.(list_mooring_IDs{i1}).spectrally.E));
% % 	%
% %     waveStats.(list_mooring_IDs{i1}).spectrally.cumulativeFlux(linlims_aux, :) = cumsum(waveStats.(list_mooring_IDs{i1}).spectrally.E(linlims_aux, :) .* ...
% %                                                                                         waveStats.(list_mooring_IDs{i1}).spectrally.cg(linlims_aux, :), 1) .* ...
% %                                                                                  waveStats.(list_mooring_IDs{i1}).df;
                                                         
	%
    disp(['--- Done with spectral energy calculation for mooring ' mooringFlux.mooringID{i1} ' ---'])
    toc(totalRunTime)
end

%
disp('--- Done with energetics calculation for all moorings ---')
toc(totalRunTime)


%%

toc(totalRunTime)
