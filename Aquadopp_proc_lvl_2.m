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
                 'E03_13300', 'E04_13172', 'E06_9736', 'E12_11150', ...
                 'F01_9995', 'F02_5838', 'F03_5384', 'F04_5401', 'F05_14032', ...
                 'X06_13290', 'X13_9945'};
% %                  'D01_12346', 'D02_0653', ...


% % A subset with a few
% list_Aquadopp = {'A03_5380', ...
%                  'B02_12507', ...
%                  'C03_0709', ...
%                  'E03_13300', ...
%                  'F01_9995', 'F02_5838', 'F03_5384'};

% % % A subset with a few
% % list_Aquadopp = {'B02_12507', 'E03_13300'};

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
freqbands.SS = [0.045, 0.3];
freqbands.IG = [0.005, 0.03];


%%

list_copyvars = {'SN', 'mooringID', 'latitude', 'longitude', ...
                 'site', 'X', 'Y', ...
                 'samplingtime', 'averagingtime', ...
                 'transducerHAB', 'binsize', 'zhab'};


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

    % First compute depth-averaged velocity from the data
    %
    aquadoppL2.dtimedata = aquadoppL1.dtime;
    aquadoppL2.pressure = aquadoppL1.pressure;
    aquadoppL2.udepthavg = mean(aquadoppL1.u, 1, 'omitnan');
    aquadoppL2.vdepthavg = mean(aquadoppL1.v, 1, 'omitnan');
    aquadoppL2.wdepthavg = mean(aquadoppL1.w, 1, 'omitnan');
    % Turn into column vectors
    aquadoppL2.pressure = aquadoppL2.pressure(:);
    aquadoppL2.udepthavg = aquadoppL2.udepthavg(:);
    aquadoppL2.vdepthavg = aquadoppL2.vdepthavg(:);
    aquadoppL2.wdepthavg = aquadoppL2.wdepthavg(:);

    % Get first and last time stamps of hourly statistics where there is
    % sufficient data at the edges on either half of the averaging windo
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

    aquadoppL2.pressuremean = time_smooth_reg(aquadoppL1.dtime, aquadoppL2.pressure, windowavg, timestatslims);
    %
    aquadoppL2.umean = time_smooth_reg(aquadoppL1.dtime, aquadoppL2.udepthavg, windowavg, timestatslims);
    aquadoppL2.vmean = time_smooth_reg(aquadoppL1.dtime, aquadoppL2.vdepthavg, windowavg, timestatslims);
    aquadoppL2.wmean = time_smooth_reg(aquadoppL1.dtime, aquadoppL2.wdepthavg, windowavg, timestatslims);




% % %     % Compute hourly mean velocities
% % %     aquadoppL2.u = NaN(length(aquadoppL2.zhab), length(aquadoppL2.dtime));
% % %     aquadoppL2.v = aquadoppL2.u;
% % %     aquadoppL2.w = aquadoppL2.u;
% % %     %
% % % % %     aquadoppL2.nptsavg = aquadoppL2.u;
% % % 
% % %     tic
% % %     % Loop over ADCP bins
% % %     for i2 = 1:length(aquadoppL2.zhab)
% % %         %
% % % % %         [aquadoppL2.u(i2, :), ~, aquadoppL2.nptsavg(i2, :)] = time_smooth_reg(aquadoppL1.dtime, aquadoppL1.u(i2, :), windowavg, timestatslims);
% % %         aquadoppL2.u(i2, :) = time_smooth_reg(aquadoppL1.dtime, aquadoppL1.u(i2, :), windowavg, timestatslims);
% % %         aquadoppL2.v(i2, :) = time_smooth_reg(aquadoppL1.dtime, aquadoppL1.v(i2, :), windowavg, timestatslims);
% % %         aquadoppL2.w(i2, :) = time_smooth_reg(aquadoppL1.dtime, aquadoppL1.w(i2, :), windowavg, timestatslims);
% % %     end
% % %     toc
% % %    
% % %     % Compute depth-averaged hourly mean velocity
% % %     aquadoppL2.udepthavg = mean(aquadoppL2.u, 1, 'omitnan');
% % %     aquadoppL2.vdepthavg = mean(aquadoppL2.v, 1, 'omitnan');
% % %     aquadoppL2.wdepthavg = mean(aquadoppL2.w, 1, 'omitnan');
    

    %%

    % ----------------------------------------------------
    % Compute spectra for DEPTH-AVERAGED u, v, and w if sampling is 1 Hz
    if aquadoppL1.samplingtime == 1

        %%
        %
        disp(['--- Sampling rate of Aquadopp ' list_Aquadopp{i1} ' is ' num2str(aquadoppL1.samplingtime) ' second(s). Computing spectra from depth-averaged velocity ---'])

        %%

        %
        aquadoppL2.ubin1 = aquadoppL1.u(1, :);
        aquadoppL2.vbin1 = aquadoppL1.v(1, :);
        aquadoppL2.wbin1 = aquadoppL1.w(1, :);
        %
        aquadoppL2.ubin1 = aquadoppL2.ubin1(:);
        aquadoppL2.vbin1 = aquadoppL2.vbin1(:);
        aquadoppL2.wbin1 = aquadoppL2.wbin1(:);
        %
        aquadoppL2.ubin1mean = time_smooth_reg(aquadoppL1.dtime, aquadoppL2.ubin1, windowavg, timestatslims);
        aquadoppL2.vbin1mean = time_smooth_reg(aquadoppL1.dtime, aquadoppL2.vbin1, windowavg, timestatslims);
        aquadoppL2.wbin1mean = time_smooth_reg(aquadoppL1.dtime, aquadoppL2.wbin1, windowavg, timestatslims);


        %%

        tic
        %
        [Suu, ~, freqvec, dof_uvw] = ...
                    spectra_scalar_reg(aquadoppL2.dtimedata, aquadoppL2.ubin1, ...
                                       windowfft, windowavg, timestatslims);
        %
        Svv = spectra_scalar_reg(aquadoppL2.dtimedata, aquadoppL2.vbin1, ...
                                 windowfft, windowavg, timestatslims);
        %
        Sww = spectra_scalar_reg(aquadoppL2.dtimedata, aquadoppL2.wbin1, ...
                                 windowfft, windowavg, timestatslims);

        toc

        %
        disp('--- Computing spectra from pressure ---')
        tic
        %
        [Spp, ~, ~, ~, avgpres] = ...
                    spectra_scalar_reg(aquadoppL2.dtimedata, aquadoppL2.pressure, ...
                                       windowfft, windowavg, timestatslims);
        toc

        % Add variables to L2 output data structure
        aquadoppL2.frequency = freqvec(:);
        %
        aquadoppL2.Suu = Suu;
        aquadoppL2.Svv = Svv;
        aquadoppL2.Sww = Sww;
        %
        aquadoppL2.Spp = Spp;

    end


    
    %% Compute elevation spectra from pressure
    if aquadoppL1.samplingtime == 1

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
        k_matrix = NaN(length(freqvec), length(aquadoppL2.dtime));
        
        % Loop over time
        tic
        for i2 = 1:length(aquadoppL2.dtime)
        
            % Only do calculation if there is a valid
            % average bottom depth at time timespec(i2)
            if ~isnan(avgbottomdepth(i2)) && (avgbottomdepth(i2)>0.2)
    
                % Loop over frequencies
                for i3 = 2:length(freqvec)
            
                    %
                    disp_rel_solver = @(k) disp_rel(k, avgbottomdepth(i2), freqvec(i3));
            
                    % In radians per meter
                    k_matrix(i3, i2) = fzero(disp_rel_solver, [(2*pi/(5000)), (2*pi/(1))]);
            
                end
            end
        end
        toc
        
        % Compute transfer function
        h_matrix = repmat(avgbottomdepth, length(freqvec), 1);
        correction = cosh(k_matrix .* h_matrix) ./ ...
                     cosh(k_matrix * aquadoppL2.transducerHAB);
        aquadoppL2.See = aquadoppL2.Spp .* correction.^2;


    end


    %% Compute variance in the sea-swell bands if sampling is 1 Hz
    if aquadoppL1.samplingtime == 1
        %
        disp('--- Integrating spectra to compute variance for each velocity component in the sea-swell band ---')

        %
        list_freqbands = fieldnames(freqbands);
        %
        for i2 = 1:length(list_freqbands)
            %
            freqlims_aux = freqbands.(list_freqbands{i2});
            %
            aquadoppL2.(['freqband' list_freqbands{i2}]) = freqlims_aux;

            %
            linband_aux = (aquadoppL2.frequency >= freqlims_aux(1)) & ...
                          (aquadoppL2.frequency  < freqlims_aux(2));


            %
            aquadoppL2.(['uvar' list_freqbands{i2}]) = ...
                              trapz(aquadoppL2.frequency(linband_aux), ...
                                    aquadoppL2.Suu(linband_aux, :), 1);
            aquadoppL2.(['vvar' list_freqbands{i2}]) = ...
                              trapz(aquadoppL2.frequency(linband_aux), ...
                                    aquadoppL2.Svv(linband_aux, :), 1);
            aquadoppL2.(['wvar' list_freqbands{i2}]) = ...
                              trapz(aquadoppL2.frequency(linband_aux), ...
                                    aquadoppL2.Sww(linband_aux, :), 1);

            %
            m0_aux = trapz(aquadoppL2.frequency(linband_aux), aquadoppL2.See(linband_aux, :), 1);
            m1_aux = trapz(aquadoppL2.frequency(linband_aux), ...
                                        (aquadoppL2.See(linband_aux, :) .* ...
                                         repmat(aquadoppL2.frequency(linband_aux), 1, length(aquadoppL2.dtime))), 1);
            %
            aquadoppL2.(['Tmean' list_freqbands{i2}]) = m0_aux./m1_aux;
            aquadoppL2.(['Hsig' list_freqbands{i2}]) = 4.*sqrt(m0_aux);
            
% %             %
% %             aquadoppL2.(['Hsig' list_freqbands{i2}]) = ...
% %                               trapz(aquadoppL2.frequency(linband_aux), ...
% %                                     aquadoppL2.See(linband_aux, :), 1);
% %             aquadoppL2.(['Hsig' list_freqbands{i2}]) = ...
% %                           4*sqrt(aquadoppL2.(['Hsig' list_freqbands{i2}]));
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

    end

    % ----------------------------------------------------
    % Organize L2 data structure
    aquadoppL2.README = ['Level 2 Aquadopp data from ROXSI 2022. The ' ...
                         'L2 data structure includes hourly statistics' ...
                         'of pressure and depth-averaged velocity. ' ...
                         'Hourly statistics are given at times dtime, ' ...
                         'which is the center time at which the hourly ' ...
                         'statistics are computed. ' ...
                         'If the Aquadopp had a sampling rate of 1 Hz, ' ...
                         'then hourly mean, hourly spectra, and variance ' ...
                         'in the sea-swell (SS) and infragravity (IG) ' ...
                         'bands are included. If instead the Aquadopp ' ...
                         'had low memory storage and recorded averaged ' ...
                         'data, then only the hourly mean is computed.'];


    %%
    % ----------------------------------------------------
    % Make QC plot of spectra
    if aquadoppL1.samplingtime == 1
        %
        hfig_spec_pcolor = figure;
    
            %
            set(hfig_spec_pcolor, 'units', 'normalized')
            set(hfig_spec_pcolor, 'Position', [0.4887, 0.1049, 0.2828, 0.4833])
            %
            haxs_1 = axes('Position', [0.15, 0.65, 0.6, 0.2]);
            haxs_2 = axes('Position', [0.15, 0.40, 0.6, 0.2]);
            haxs_3 = axes('Position', [0.15, 0.15, 0.6, 0.2]);
            %
            haxs_all = [haxs_1, haxs_2, haxs_3];
            %
            for i = 1:length(haxs_all)
                hold(haxs_all(i), 'on')
            end
    
                %
                pcolor(haxs_1, aquadoppL2.dtime, aquadoppL2.frequency, log10(aquadoppL2.Suu))
                pcolor(haxs_2, aquadoppL2.dtime, aquadoppL2.frequency, log10(aquadoppL2.Sww))
                %
                hp_u = plot(haxs_3, aquadoppL2.dtime, aquadoppL2.uvarSS, '-');
                hp_v = plot(haxs_3, aquadoppL2.dtime, aquadoppL2.vvarSS, '-');
                hp_w = plot(haxs_3, aquadoppL2.dtime, 10.*aquadoppL2.wvarSS, '-');
    
            %
            shading(haxs_1, 'flat')
            shading(haxs_2, 'flat')

            %
            hcb_1 = colorbar(haxs_1);
                hcb_1.Position = [(sum(haxs_1.Position([1, 3]))+0.02), haxs_1.Position(2), 0.025, haxs_1.Position(4)];
                hcb_1.Label.Interpreter = 'Latex';
                hcb_1.Label.String = '[$\log_{10}$ m$^2$ s$^{-2}$ Hz$^{-1}$]';
            hcb_2 = colorbar(haxs_2);
                hcb_2.Position = [(sum(haxs_2.Position([1, 3]))+0.02), haxs_2.Position(2), 0.025, haxs_2.Position(4)];
                hcb_2.Label.Interpreter = 'Latex';
                hcb_2.Label.String = '[$\log_{10}$ m$^2$ s$^{-2}$ Hz$^{-1}$]';
                
            %
            caxis(haxs_1, [-3, 1])
            caxis(haxs_2, [-4, 0])

            %
            set(haxs_all, 'FontSize', 12, 'Box', 'on', ...
                          'XGrid', 'on', 'YGrid', 'on', ...
                          'XLim', aquadoppL2.dtime([1, end]))
            set(haxs_all(1:2), 'YScale', 'log', 'YLim', aquadoppL2.frequency([2, end]))
    
            %
            ylabel(haxs_1, 'Frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 10)
            ylabel(haxs_2, 'Frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 10)
            ylabel(haxs_3, '[m$^2$ s$^{-2}$]', 'Interpreter', 'Latex', 'FontSize', 10)
    
            %
            title(haxs_1, {['ROXSI 2022: Aquadopp ' list_Aquadopp{i1}(1:3) ' SN ' list_Aquadopp{i1}(5:end) '']; ...
                           'Spectra of depth-averaged u and w.'; 'Sea-swell variance of $u$, $v$, $10 \times w$'}, ...
                           'Interpreter', 'Latex', 'FontSize', 10)
    
            %
            linkaxes(haxs_all, 'x')
    
            %
            for indhaxs = [1, 2]
                xlim_aux = xlim(haxs_all(indhaxs));
                plot(haxs_all(indhaxs), xlim_aux, aquadoppL2.freqbandSS(1).*[1, 1], '--k')
                plot(haxs_all(indhaxs), xlim_aux, aquadoppL2.freqbandSS(2).*[1, 1], '--k')
                xlim(haxs_all(indhaxs), xlim_aux)
            end

            %
            hleg = legend([hp_u, hp_v, hp_w], 'var($u$)', 'var($v$)', '$10 \times$ var($w$)', 'Location', 'EastOutside');
                hleg.Interpreter = 'Latex';
                hleg.Position = [0.85, 0.215, 0.05, 0.05];
    end


    %%
    % ----------------------------------------------------
    % Save L2 data structure
    %
    %
    disp('----- Saving level 2 data -----')
    %
    str_filename = ['roxsi_aquadopp_L2_' char(aquadoppL2.mooringID) '_' char(aquadoppL2.SN)];
    save(fullfile(dir_output_data_L2, [str_filename '.mat']), 'aquadoppL2', '-v7.3')
    %
    disp('----- Done saving data -----')

    % ----------------------------------------------------
    % Save QC plots
    if aquadoppL1.samplingtime == 1
        disp('----- Saving level 2 QC plots -----')
        %
        str_filename = ['uvw_spectra_' char(aquadoppL2.mooringID) '_' char(aquadoppL2.SN)];
        % Save figure as *.png
        exportgraphics(hfig_spec_pcolor, fullfile(dir_output_data_L2, [str_filename '.png']), 'Resolution', 300)
    end


    %% Print progress message

    %
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
