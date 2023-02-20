%% Signature1000_proc_lvl_2
%
% Since:
%   - This L2 processing for Signatures is basically
%     the same as for Aquadopps and;
%   - there is a lot of Signature data that may require
%     doing analyses on different files;
% I will write this script and a higher level function
% to do the L2 processing on Signature (or maybe any ADCP) data.
%
% If I take depth-averaged velocity, or in the first bin,
% I can probably patch all of the timeseries together
% and then do the analysis.


clear
close all


%%
% --------------------------------------------
% ----------------- PREAMBLE -----------------
% --------------------------------------------

%% Directory of L1 data

%
dirparent_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/';
dir_dataL1 = fullfile(dirparent_data, 'Level1_Data', 'Signature_Level1');


%% Directory to output L2 data and figures

%
dir_output_parent = dirparent_data;
dir_dataL2 = fullfile(dir_output_parent, 'Level2_Data', 'Signature_Level2');


%%
% --------------------------------------------
% ---------- PROCESSING PARAMETERS -----------
% --------------------------------------------

%%

% All Signatures
list_Signature = {'B10_103045', ...    % Do B first, because they have less
                  'B13_103046', ...    % data than A01, so if something is
                  'B15_103056', ...    % wrong with the code, the error
                  'B17_101923', ...    % pops up sooner
                  'A01_103043', ...
                  'C01_102128', ...
                  'X05_100231', ...
                  'X11_101941'};

Nsignatures = length(list_Signature);


%%

% Time window for statistical averaging (in seconds)
windowavg = 60*60;

% Time window for individual fft (in seconds)
windowfft = 6*60;
    
% Sea-swell frequency band
freq_lims = [0.045, 0.3];
freqbands.SS = [0.045, 0.3];
freqbands.IG = [0.005, 0.03];

%
frequencytrimTH = 0.35;

%%

list_copyvars = {'SN', 'mooringID', 'latitude', 'longitude', ...
                 'site', 'X', 'Y', ...
                 'transducerHAB', 'binsize', 'zhab'};


%%
% --------------------------------------------
% ------------- DEFINE TIME GRID -------------
% --------------------------------------------

%% Define this instead of using the data
% (it's simpler this way)

% All ADCPs programmed to start at 2022/06/21 18:00:00
% (B10 not in the water at start time)
timestatslims = [datetime(2022, 06, 21, 19, 00, 00), ...
                 datetime(2022, 07, 25, 10, 00, 00)];
% most ADCPs recovered by the 21st at 9:00AM.

%
timegrid_L2 = timestatslims(1) : hours(1) : timestatslims(2);
timegrid_L2.TimeZone = 'America/Los_Angeles';


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------------- DO DATA PROCESSING ----------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------

%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_signature_procL2_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_dataL2, log_file_name))
%
totalRunTime = tic;


%% Display on the screen:

%
disp(' '), disp(' ')
disp('------------------------------ Processing L2 data from Signatures ------------------------------')
disp('List of Signatures being processed:')
%
for i = 1:Nsignatures
    disp([num2str(i) ' - ' list_Signature{i}])
end


%% Do L2 processing

% Loop over instruments
for i1 = 1:Nsignatures


    %% Load data and get reduced data
    
    %
    list_dirsegments = dir(fullfile(dir_dataL1, 'segment_*'));
    
    % Loop over segments
    for i2 = 1:length(list_dirsegments)

        % Check if there is a file in the i2'th segment
        % (it may not have if the data from different
        % ADCPs are trimmed differently)
        filedir_aux = dir(fullfile(list_dirsegments(i2).folder, ...
                                   list_dirsegments(i2).name, ...
                                   ['roxsi_signature_L1_' list_Signature{i1} '.mat']));

        % If there is no file, skip/continue to the next loop iteration
        if isempty(filedir_aux)
            continue
        end
        
        %
        sigL1 = load(fullfile(filedir_aux.folder, filedir_aux.name));
        sigL1 = sigL1.sigL1;

        %
        if ~exist('sigL2', 'var')
            %
            for i3 = 1:length(list_copyvars)
                sigL2.(list_copyvars{i3}) = sigL1.(list_copyvars{i3});
            end
    % % %             sigL2.SN = sigL1.SN;
    % % %             sigL2.mooringID = sigL1.mooringID;
    % % %             % ... others ...
    % % %             %
    % % %             sigL2.transducerHAB = sigL1.transducerHAB;
    % % %             %
    % % %             sigL2.zhab = sigL1.zhab;

            %
            sigL2.dtimedata = sigL1.dtime;
            %
            sigL2.pressure = sigL1.pressure;
            %
            sigL2.udepthavg = mean(sigL1.u, 1, 'omitnan');
            sigL2.vdepthavg = mean(sigL1.v, 1, 'omitnan');
            sigL2.wdepthavg = mean(sigL1.w, 1, 'omitnan');
            %
            sigL2.ubin1 = sigL1.u(1, :);
            sigL2.vbin1 = sigL1.v(1, :);
            sigL2.wbin1 = sigL1.w(1, :);

            %
            sigL2.udepthavg = sigL2.udepthavg(:);
            sigL2.vdepthavg = sigL2.vdepthavg(:);
            sigL2.wdepthavg = sigL2.wdepthavg(:);
            %
            sigL2.ubin1 = sigL2.ubin1(:);
            sigL2.vbin1 = sigL2.vbin1(:);
            sigL2.wbin1 = sigL2.wbin1(:);
        end
        
        %
        if i2==1
            

        %
        else

            % Find matching grid points to concatenate data
            ind_match_aux = find(sigL1.dtime == sigL2.dtimedata(end));
            ind_start_aux = ind_match_aux + 1;

            %
            sigL2.dtimedata = [sigL2.dtimedata; sigL1.dtime(ind_start_aux:end)];
            %
            sigL2.pressure = [sigL2.pressure; sigL1.pressure(ind_start_aux:end)];
            %
            sigL2.udepthavg = [sigL2.udepthavg; mean(sigL1.u(:, ind_start_aux:end), 1, 'omitnan').'];
            sigL2.vdepthavg = [sigL2.vdepthavg; mean(sigL1.v(:, ind_start_aux:end), 1, 'omitnan').'];
            sigL2.wdepthavg = [sigL2.wdepthavg; mean(sigL1.w(:, ind_start_aux:end), 1, 'omitnan').'];
            %
            sigL2.ubin1 = [sigL2.ubin1; sigL1.u(1, ind_start_aux:end).'];
            sigL2.vbin1 = [sigL2.vbin1; sigL1.v(1, ind_start_aux:end).'];
            sigL2.wbin1 = [sigL2.wbin1; sigL1.w(1, ind_start_aux:end).'];
        end

        % Remove sigL1 variable
        clear sigL1
    end
  
    
    %% Compute hourly depth-averaged quantities

    %
    sigL2.dtime = timegrid_L2;

    %
    sigL2.pressuremean = time_smooth_reg(sigL2.dtimedata, sigL2.pressure, windowavg, sigL2.dtime([1, end]));
    %
    sigL2.udepthavgmean = time_smooth_reg(sigL2.dtimedata, sigL2.udepthavg, windowavg, timestatslims);
    sigL2.vdepthavgmean = time_smooth_reg(sigL2.dtimedata, sigL2.vdepthavg, windowavg, timestatslims);
    sigL2.wdepthavgmean = time_smooth_reg(sigL2.dtimedata, sigL2.wdepthavg, windowavg, timestatslims);
    %
    sigL2.ubin1mean = time_smooth_reg(sigL2.dtimedata, sigL2.ubin1, windowavg, timestatslims);
    sigL2.vbin1mean = time_smooth_reg(sigL2.dtimedata, sigL2.vbin1, windowavg, timestatslims);
    sigL2.wbin1mean = time_smooth_reg(sigL2.dtimedata, sigL2.wbin1, windowavg, timestatslims);
    

    %% Compute spectra
    

    tic
    %
    [sigL2.Suu, ~, sigL2.frequency, sigL2.DOF_uvw] = ...
                spectra_scalar_reg(sigL2.dtimedata, sigL2.udepthavg, ...
                                   windowfft, windowavg, timestatslims);
    %
    sigL2.Svv = spectra_scalar_reg(sigL2.dtimedata, sigL2.vdepthavg, ...
                                   windowfft, windowavg, timestatslims);
    %
    sigL2.Sww = spectra_scalar_reg(sigL2.dtimedata, sigL2.wdepthavg, ...
                                   windowfft, windowavg, timestatslims);
    %
    sigL2.Spp = spectra_scalar_reg(sigL2.dtimedata, sigL2.pressure, ...
                                   windowfft, windowavg, timestatslims);
    toc
    
    %
    sigL2.frequency = sigL2.frequency(:);


    %% Trim out unnecessary very high frequencies

    %
    lkeepfrequency = (sigL2.frequency <= frequencytrimTH);
    %
    sigL2.frequency = sigL2.frequency(lkeepfrequency);
    sigL2.Suu = sigL2.Suu(lkeepfrequency, :);
    sigL2.Svv = sigL2.Svv(lkeepfrequency, :);
    sigL2.Sww = sigL2.Sww(lkeepfrequency, :);
    %
    sigL2.Spp = sigL2.Spp(lkeepfrequency, :);


    %% Compute elevation spectra from pressure
    
    %
    g = 9.8;
    rho0 = 1030;
    
    % Computes bottom depth from average pressure (which is
    % hydrostatic over 1 hour)
    sigL2.bottomdepthmean = 1e4*sigL2.pressuremean ./ (rho0*g);
    
    % Convert frequencies to wavenumbers using
    % linear wave theory
    %
    % Function handle to compute wavenumbers
    % (PS: frequency in Hz and k in radians per meter)
    disp_rel = @(k, H, freq) g*k*tanh(k*H) - (2*pi*freq)^2;
    
    
    % Compute wavenumbers (in radians per meter) from linear wave theory
    k_matrix = NaN(length(sigL2.frequency), length(sigL2.dtime));
    
    % Loop over time
    tic
    for i2 = 1:length(sigL2.dtime)
    
        % Only do calculation if there is a valid
        % average bottom depth at time timespec(i2)
        if ~isnan(sigL2.bottomdepthmean(i2))

            % Loop over frequencies
            for i3 = 2:length(sigL2.frequency)
        
                %
                disp_rel_solver = @(k) disp_rel(k, sigL2.bottomdepthmean(i2), sigL2.frequency(i3));
        
                % In radians per meter
                k_matrix(i3, i2) = fzero(disp_rel_solver, [(2*pi/(5000)), (2*pi/(1))]);
        
            end
        end
    end
    toc
    
    % Compute transfer function
    h_matrix = repmat(sigL2.bottomdepthmean, length(sigL2.frequency), 1);
    correction = cosh(k_matrix .* h_matrix) ./ ...
                 cosh(k_matrix * sigL2.transducerHAB);
    sigL2.See = sigL2.Spp .* correction.^2;


    %% Compute statistics in the sea-swell and infragravity bands

    %
    disp('--- Integrating spectra to compute variance for each velocity component in the sea-swell band ---')

    %
    list_freqbands = fieldnames(freqbands);
    %
    for i2 = 1:length(list_freqbands)
        %
        freqlims_aux = freqbands.(list_freqbands{i2});
        %
        sigL2.(['freqband' list_freqbands{i2}]) = freqlims_aux;

        %
        linband_aux = (sigL2.frequency >= freqlims_aux(1)) & ...
                      (sigL2.frequency  < freqlims_aux(2));
% % %         %
% % %         sigL2.(['uvar' list_freqbands{i2}]) = ...
% % %                           trapz(aquadoppL2.frequency(linband_aux), ...
% % %                                 aquadoppL2.Suu(linband_aux, :), 1);
% % %         sigL2.(['vvar' list_freqbands{i2}]) = ...
% % %                           trapz(aquadoppL2.frequency(linband_aux), ...
% % %                                 aquadoppL2.Svv(linband_aux, :), 1);
% % %         sigL2.(['wvar' list_freqbands{i2}]) = ...
% % %                           trapz(aquadoppL2.frequency(linband_aux), ...
% % %                                 aquadoppL2.Sww(linband_aux, :), 1);

        %
        m0_aux = trapz(sigL2.frequency(linband_aux), sigL2.See(linband_aux, :), 1);
        m1_aux = trapz(sigL2.frequency(linband_aux), ...
                                    (sigL2.See(linband_aux, :) .* ...
                                     repmat(sigL2.frequency(linband_aux), 1, length(sigL2.dtime))), 1);
        %
        sigL2.(['Tmean' list_freqbands{i2}]) = m0_aux./m1_aux;
        sigL2.(['Hsig' list_freqbands{i2}]) = 4.*sqrt(m0_aux);

    end
    
% %         %
% %         linfreqband_seaswell = (aquadoppL2.frequency >= aquadoppL2.freqband(1)) & ...
% %                                (aquadoppL2.frequency  < aquadoppL2.freqband(2));
% % 
% %         %
% %         aquadoppL2.uvar = trapz(aquadoppL2.frequency(linfreqband_seaswell), ...
% %                                 aquadoppL2.Suu(linfreqband_seaswell, :), 1);
% %         %
% %         aquadoppL2.vvar = trapz(aquadoppL2.frequency(linfreqband_seaswell), ...
% %                                 aquadoppL2.Svv(linfreqband_seaswell, :), 1);
% %         %
% %         aquadoppL2.wvar = trapz(aquadoppL2.frequency(linfreqband_seaswell), ...
% %                                 aquadoppL2.Sww(linfreqband_seaswell, :), 1);

    
    %% Organize L2 data structure

    %
    sigL2.README = ['Level 2 Signature1000 data from the large-scale ' ...
                    'array in ROXSI 2022. The ' ...
                    'L2 data structure includes hourly statistics' ...
                    'of pressure, depth-averaged velocity, and velocity in the first bin. ' ...
                    'Hourly statistics are given at times dtime, ' ...
                    'which is the center time at which the hourly ' ...
                    'statistics are computed. ' ...
                    'Bulk statistics for sea-swell (SS) and infragravity ' ...
                    '(IG) bands are also included'];

    %%
    % ----------------------------------------------------
    % Save L2 data structure
    disp('----- Saving level 2 data -----')
    %
    str_filename = ['roxsi_signature_L2_' char(sigL2.mooringID) '_' char(sigL2.SN)];
    save(fullfile(dir_dataL2, [str_filename '.mat']), 'sigL2', '-v7.3')
    %
    disp('----- Done saving data -----')

% %     % ----------------------------------------------------
% %     % Save QC plots
% %     if aquadoppL1.samplingtime == 1
% %         disp('----- Saving level 2 QC plots -----')
% %         %
% %         str_filename = ['uvw_spectra_' char(aquadoppL2.mooringID) '_' char(aquadoppL2.SN)];
% %         % Save figure as *.png
% %         exportgraphics(hfig_spec_pcolor, fullfile(dir_output_data_L2, [str_filename '.png']), 'Resolution', 300)
% %     end


    %% Print progress message

    %
    disp(['----- Done with L2 data processing for Signature ' list_Signature{i1} ' -----'])

    %
    toc(totalRunTime)

    %
    clear sigL2


end
 

%% Print message the processing is done and close log file

%
disp('###################### Done with L2 data processing for all Signature1000 in the list ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');