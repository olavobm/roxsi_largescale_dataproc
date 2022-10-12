%% Script to process the level 2 Smart Mooring data.

clear
close all

%%
% --------------------------------------
% -------------- PREAMBLE --------------
% --------------------------------------

%%

% % roxsi_add_libraries()

%% Directory where data will be loaded from

%
% % dirparent_data = data_dirpath();
dirparent_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
%
dir_L1data_parent = fullfile(dirparent_data, 'Level1_Data', 'Spotters_Smart', 'gridded');


%% Directory where output will be saved

%
% % dir_output_parent = data_dirpath();
dir_output_parent = pwd;
%
dir_output_data_L2 = dir_output_parent;


%%
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------

%%

% list_smartmoorings = {'E01_1851', ...
%                       'E02_1859', ...
%                       'E05_1853', ...
%                       'E07_1855', ...
%                       'E07_1857', ...
%                       'E08_1852', ...
%                       'E09_1850', ...
%                       'E09_1856', ...
%                       'E10_1848', ...
%                       'E11_1860', ...
%                       'E13_1849'};

% Just a test
list_smartmoorings = {'E02_1859'};

%
Nspotsmart = length(list_smartmoorings);


%%

%
disp(' '), disp(' ')
disp('------------------------------ Processing data from Smart Moorings to L2 ------------------------------')
disp('List of Smart Moorings being processed:')
%
for i = 1:Nspotsmart
    disp([num2str(i) ' - ' list_smartmoorings{i}])
end


%% Do L2 processing

% Loop over Smart Moorings in the list
for i1 = 1:Nspotsmart

    %%
    %
    disp(' '), disp(' ')
    disp(['----- Start processing Smart Mooring: ' list_smartmoorings{i1} ' -----'])


    %% Load the data from the i1'th Smart Mooring

    %
    filename_aux = ['smart_mooring_' list_smartmoorings{i1}(1:3) 'sp_' ...
                    list_smartmoorings{i1}(end-3:end) '_L1_gridded.mat'];
    %
    spotsmartL1 = load(fullfile(dir_L1data_parent, filename_aux));
    spotsmartL1 = spotsmartL1.spotsmart;


    %% Copy some metadata from L1 to L2
    
    %
    L2out.mooringID = spotsmartL1.mooringID;
    L2out.SN = spotsmartL1.SN;
    
    %
    mooringLoc = ROXSI_mooringlocation('E01sp');
    
    L2out.site = mooringLoc.roxsiarray;
    
    %
    L2out.latitude = spotsmartL1.latitude;
    L2out.longitude = spotsmartL1.longitude;
    
    
    %% Copute x/y coordinates
    
    [L2out.X, L2out.Y] = ROXSI_lltoxy(L2out.latitude, L2out.longitude, L2out.site);
    
    %% Height of the sensor above the bottom
    
    %
    L2out.zhab = spotsmartL1.zhab;
    
    
    %% Parameters for computing spectrum
    
    % Both in seconds
    windowfft = 6*60;
    windowavg = 60*60;
    
    %
    timespeclims = [datetime(2022, 06, 17, 19, 0, 0), ...
                    datetime(2022, 07, 20, 05, 0, 0)];
    
    %% Frequency limits of wave bands of interest
    
    freq_bands.IG = [0.005, 0.03];
    freq_bands.SS = [0.045, 0.3];
    
    
    %%
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
    
    %% Initialize a log file with what is printed to the
    % command window and timer for running the whole script
    
    %
    log_file_name = ['log_SmartMooring_procL2_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
    %
    diary(fullfile(dir_output_data_L2, log_file_name))
    
    %
    totalRunTime = tic;
    
    
    %% Computes pressure spectra
    
    %
    tic
    [Spp, timespec, freqvec, dof, avgpres] = ...
                spectra_scalar_reg(spotsmart.dtime, spotsmart.pressure, ...
                                   windowfft, windowavg, timespeclims);
    toc
    
    
    %% Computes SSH spectra (or do it inside calculation above so that
    % the average pressure/depth is used at windowfft instead of windowavg?????).
    % Falk's code for SoloDs uses the hourly averaged pressure. The other option
    % would be to write a function similar to the code above, but include
    % the SSH conversion inside it
    
    %
    g = 9.8;
    rho0 = 1030;
    
    % Computes bottom depth from average pressure (which is
    % hydrostatic over 1 hour)
    avgbottomdepth = 1e4*avgpres ./ (rho0*g);
    
    % Convert frequencies to wavenumbers using
    % linear wave theory
    %
    % Function handle to compute wavenumbers
    % (PS: frequency in Hz and k in radians per meter)
    disp_rel = @(k, H, freq) g*k*tanh(k*H) - (2*pi*freq)^2;
    
    
    %Compute wavenumbers (in radians per meter) from linear wave theory
    k_matrix = NaN(length(freqvec), length(timespec));
    
    % Loop over time
    tic
    for i2 = 1:length(timespec)
    
        % Loop over frequencies
        for i3 = 2:length(freqvec)
    
            %
            disp_rel_solver = @(k) disp_rel(k, avgbottomdepth(i2), freqvec(i3));
    
            % In radians per meter
            k_matrix(i3, i2) = fzero(disp_rel_solver, [(2*pi/(5000)), (2*pi/(1))]);
    
        end
    end
    toc
    
    % Compute transfer function
    h_matrix = repmat(avgbottomdepth, length(freqvec), 1);
    correction = cosh(k_matrix .* h_matrix) ./ ...
                 cosh(k_matrix * spotsmartL1.zhab);
    See = Spp .* correction.^2;
    
    
    %% Now compute significant wave height
    % at sea-swell and infragravity wave bands
    
    % Frequency resolution (in Hz)
    df = freqvec(2) - freqvec(1);
    
    %
    lin_SSband = ((freqvec > freq_bands.SS(1)) & (freqvec < freq_bands.SS(2)));
    lin_IGband = ((freqvec > freq_bands.IG(1)) & (freqvec < freq_bands.IG(2)));
    
    %
    Hsig_SS = 4*sqrt(sum(See(lin_SSband, :), 1).*df);
    Hsig_IG = 4*sqrt(sum(See(lin_IGband, :), 1).*df);
    
    % Plot pcolor of SSH spectra up until upper band
    % of sea-swell wave, to check the spectra are not
    % blowing up
    ltrimspec = (freqvec <= freq_bands.SS(2));
    
    % % % % QC figures
    % % % %
    % % % figure
    % % %     %
    % % %     pcolor(timespec, freqvec(ltrimspec), log10(See(ltrimspec, :)))
    % % %     shading flat
    % % %     %
    % % %     hold on
    % % % 
    % % %     %
    % % %     plot(timespec([1, end]), freq_bands.SS(1).*[1, 1], '--r')
    % % %     plot(timespec([1, end]), freq_bands.SS(2).*[1, 1], '--r')
    % % %     plot(timespec([1, end]), freq_bands.IG(1).*[1, 1], '--k')
    % % %     plot(timespec([1, end]), freq_bands.IG(2).*[1, 1], '--k')
    % % % 
    % % %     %
    % % %     ylim_aux = ylim;
    % % %     ylim([ylim_aux(1), (freq_bands.SS(2) + 3*df)])
    % % % 
    % % % 
    % % %     %
    % % %     colorbar
    % % % 
    % % % %
    % % % figure
    % % %     %
    % % %     plot(freqvec, mean(See, 2), '-k')
    % % %     hold on
    % % %     plot(freqvec(ltrimspec), mean(See(ltrimspec, :), 2), '.-b')
    % % % 
    % % %     %
    % % %     set(gca, 'FontSize', 16, 'Box', 'on', ...
    % % %              'XGrid', 'on', 'YGrid', 'on', ...
    % % %              'XScale', 'log', 'YScale', 'log')
    % % %     %
    % % %     set(gca, 'XLim', [(df/1.5), (freq_bands.SS(2) + 10*df)])
    % % %     set(gca, 'YLim', [1e-3, 1.5])
    % % % 
    % % %     %
    % % %     ylim_aux = ylim;
    % % %     %
    % % %     plot(freq_bands.SS(1).*[1, 1], ylim_aux, '--r')
    % % %     plot(freq_bands.SS(2).*[1, 1], ylim_aux, '--r')
    % % %     plot(freq_bands.IG(1).*[1, 1], ylim_aux, '--k')
    % % %     plot(freq_bands.IG(2).*[1, 1], ylim_aux, '--k')
    
    
    % Combine the 2 plots above in the same figure
    %
    fig_QC = figure;
        %
        set(gcf, 'Units', 'normalized')
        set(gcf, 'Position', [0.2, 0.2, 0.4, 0.6])
        %
        haxs_1 = axes('Position', [0.15, 0.50, 0.7, 0.4]);
        haxs_2 = axes('Position', [0.15, 0.15, 0.7, 0.3]);
        hold(haxs_1, 'on')
        hold(haxs_2, 'on')
    
            %
            pcolor(haxs_1, freqvec(ltrimspec), timespec, log10(See(ltrimspec, :)).')
            shading(haxs_1, 'flat')
            %
            plot(haxs_1, freq_bands.SS(1).*[1, 1], timespec([1, end]), '--r')
            plot(haxs_1, freq_bands.SS(2).*[1, 1], timespec([1, end]), '--r')
            plot(haxs_1, freq_bands.IG(1).*[1, 1], timespec([1, end]), '--k')
            plot(haxs_1, freq_bands.IG(2).*[1, 1], timespec([1, end]), '--k')
        
            %
            plot(haxs_2, freqvec, mean(See, 2, 'omitnan'), '-k')        
            plot(haxs_2, freqvec(ltrimspec), mean(See(ltrimspec, :), 2), '.-b')
            %
            plot(haxs_2, freq_bands.SS(1).*[1, 1], ylim_aux, '--r')
            plot(haxs_2, freq_bands.SS(2).*[1, 1], ylim_aux, '--r')
            plot(haxs_2, freq_bands.IG(1).*[1, 1], ylim_aux, '--k')
            plot(haxs_2, freq_bands.IG(2).*[1, 1], ylim_aux, '--k')
    
            %
            hcb = colorbar(haxs_1);
                hcb.Position = [(sum(haxs_1.Position([1, 3])) + 0.01), ...
                                haxs_1.Position(2), 0.02, haxs_1.Position(4)];
                hcb.Ticks = -4:1:1;
                hcb.Label.Interpreter = 'Latex';
                hcb.Label.String = '$\log_{10}$ of m$^2$ Hz$^{-1}$';
    
            %
            set([haxs_1, haxs_2], 'FontSize', 16, 'Box', 'on', ...
                                  'XGrid', 'on', 'YGrid', 'on')
            set([haxs_1, haxs_2], 'XScale', 'log')
            set(haxs_2, 'XScale', 'log', 'YScale', 'log')
            %
            freq_plt_lims = [(df/1.5), (freq_bands.SS(2) + 10*df)];
            %
            set([haxs_1, haxs_2], 'XLim', freq_plt_lims)
            %
            set(haxs_1, 'YLim', timespec([1, end]))
            set(haxs_2, 'YLim', [1e-3, 1.5])
    
    
        %
        xlabel(haxs_2, 'Frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 12)
        ylabel(haxs_1, 'Time ', 'Interpreter', 'Latex', 'FontSize', 12)
        ylabel(haxs_2, '[m$^2$ Hz$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', 12)
        %
        title(haxs_1, {['ROXSI 2022 pressure-based surface elevation spectra ' ...
                        '(and time-mean spectrum)'];  ...
                       ['from smart mooring ' ...
                       char(L2out.mooringID) '. Frequencies higher than ' ...
                       'sea-swell upper bound were trimmed']}, ...
                       'Interpreter', 'Latex', 'FontSize', 12)
        
    
    
    
    %% Compute mean and peak frequency in sea-swell band
    
    %
    ind_SSband = find(lin_SSband);
    
    % 
    [~, ind_peak] = max(See(lin_SSband, :), [], 1);
    % Go from indices of See(lin_SSband, :) to indices of See
    ind_peak = ind_SSband(ind_peak);
    
    %
    peakfreqSS = freqvec(ind_peak);
    
    %
    freq_matrix = repmat(freqvec(:), 1, length(timespec));
    %
    meanfreqSS = sum(See(lin_SSband, :) .* freq_matrix(lin_SSband, :), 1) ./ ...
                 sum(See(lin_SSband, :), 1);
    
    
    %% Add variables to L2 data structure
    
    %
    L2out.dtime = timespec;
    L2out.frequency = freqvec;
    
    %
    L2out.mean_depth = avgbottomdepth;
    
    %
    % % L2out.nptsfft = windowfft / samplingfreq
    
    %
    L2out.Spp = Spp;
    L2out.See = See;
    
    %
    L2out.DOF = dof(:);
    
    %
    L2out.freqbandSS = freq_bands.SS;
    L2out.freqbandIG = freq_bands.IG;
    
    %
    L2out.HsigSS = Hsig_SS(:);
    L2out.HsigIG = Hsig_IG(:);
    
    %
    L2out.meanfreqSS = meanfreqSS(:);
    L2out.peakfreqSS = peakfreqSS(:);
    
    
    %% Add README
    
    %
    time_dataproc = datetime('now', 'TimeZone', 'Local');
    time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');
    
    %
    L2out.README = ['ROXSI Level 2 Smart Mooring data. Serial number ' ...
                    '(SN) is from the Spotter buoy of the smart mooring ' ...
                    'deployed at mooringID (it''s not  the SN from the ' ...
                    'pressure sensor). X and Y are cartesian coordinates ' ...
                    'in a local coordinate system of the deployment site ' ...
                    'Data contains hourly estimates of mean_depth (m), ' ...
                    'pressure spectra (dbar^2 / Hz), and sea surface ' ...
                    'elevation spectra (m^2 / Hz). A spectrum is computed ' ...
                    'every ' num2str(windowfft) ' seconds (with 50% ' ...
                    'overlap) and averaged over an hour. Significant ' ...
                    'wave height (m) is given for sea-swell (SS) and ' ...
                    'infragravity (IG) bands ' ...
                    '(and the frequency bands for each are in Hz). Mean ' ...
                    'and peak frequencies are also computed for the SS ' ...
                    'band. DOF is the number of degrees of freedom for ' ...
                    'the spectral estimates (smaller DOF indicates ' ...
                    'gaps in the original pressure data). Data processed ' ...
                    'by ' mfilename() '.m on ' time_dataproc_char ' (TimeZone ' ...
                    time_dataproc.TimeZone ').'];
    
    
    %% Save L2 data and QC figure

    %
    spotsmartL2 = L2out;

    % ---------------------------------------------
    %
    disp('----- Saving level 2 data -----')

    %
    str_filename = ['roxsi_smartmooring_L2_' char(spotsmartL2.mooringID) '_' char(spotsmartL2.SN)];
    %
    save(fullfile(dir_output_data_L2, [str_filename '.mat']), 'sig1000', '-v7.3')

    %
    disp('----- Done saving data -----')


    % ---------------------------------------------
    %
    disp('----- Saving level 2 QC plot -----')

    %
    str_filename = ['roxsi_smartmooring_L2_' char(spotsmartL2.mooringID) '_' char(spotsmartL2.SN) '_QC_fig'];
    % Save figure as *.png
    exportgraphics(fig_QC, fullfile(dir_output_data_L2, [str_filename '.png']), 'Resolution', 300)


    %% Print progress message

    %
    disp(' '), disp(' ')
    disp(['----- Done with processing Smart Mooring: ' list_smartmoorings{i1} ' -----'])


end


%%
%
disp('###################### Done with L2 data processing for all Smart Moorings ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');






