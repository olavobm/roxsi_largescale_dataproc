%% Compute directional spectra, mean direction, and directional
% spread from ADCP data.

clear
close all

%
% % addpath(genpath(['/Users/olavobm/Documents/ROXSI_Postdoc' ...
% %                  '/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/']))


%%
% ---------------------------------------------
% ---------------------------------------------
% ---------------------------------------------


%% Directory with data files of data that
% has been prepared for this script

% % % 
% % dir_dataaux = ['/Users/olavobm/Documents/ROXSI_Postdoc' ...
% %                '/MyResearch/figures_bydate/2023_02_01/'];
% % %
% % dir_outputdata = ['/Users/olavobm/Documents/ROXSI_Postdoc' ...
% %                   '/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022' ...
% %                   '/code_proc/directional_properties/'];

% 
dir_dataaux = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results_3/data_for_dirspectra/';
%
% dir_outputdata = dir_dataaux;
dir_outputdata = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results_3/data_with_dirspectra/';


%% Moorings to use for calculation

% All moorings
list_moorings = {'B10', 'B13', 'B15', 'B17', 'A01', 'C01', 'X05', 'X11', ...   % Signatures
                 'B11', 'X13', ...    % Aquadopps, 2 MHz
                 'B08', 'E03', 'E04', 'E06', 'E12', 'X06'};     % Aquadopps, 1 MHz
             
             
%% Define time grid for computing spectra

%
dt_spec = hours(1);

% % %
% % dtimegridlims = [datetime(2022, 06, 27, 12, 00, 00), ...
% %                  datetime(2022, 06, 30, 18, 00, 00)];

% Full timeseries (excluding remaining data after recovery started)
dtimegridlims = [datetime(2022, 06, 21, 20, 00, 00), ...
                 datetime(2022, 07, 20, 05, 00, 00)];

%
dtimegrid = dtimegridlims(1) : dt_spec : dtimegridlims(2);
dtimegrid.TimeZone = 'America/Los_Angeles';

% % % Small segment 
% % dtimegrid = dtimegrid(1:40);


%% Define FFT properties -- time length of
% individual FFTs, averaging time window,
% and high-frequency cut-off

%
dt_fftsegment = 120;    % in seconds
% % dt_fftsegment = 240;    % in seconds
dt_avgfft = seconds(dt_spec);    % in seconds

% The cut-off frequency (in Hz) to use because
% the transfer function blows up 
% freq_cutoff = 0.2;
freq_cutoff = 0.3;    % this is probably good for most moorings
                      % (and it also may give results at harmonics
                      % in the range 0.2-0.3 Hz)


%% Parameters for directional spectra.
% Angle resolution and list of methods
% for directional spectra.

% Angle resolution for directional spectra
dtheta = 1;    % in degrees

%
% dspec_method = ["EMEM","IMLM","MLM","MEM"];    % the directional algorithms to use
dspec_method = "EMEM";    % the directional algorithms to use


%% Define frequency bands for bulk calculations

%
bulkbands(1).ID = 'swell';
bulkbands(1).freqlims = [0.06, 0.1];
% % bulkbands(1).linlims = (adcpdirspectra.frequency >= bulkbands(1).freqlims(1)) & ...
% %                        (adcpdirspectra.frequency <= bulkbands(1).freqlims(2));
%
bulkbands(2).ID = 'sea';
bulkbands(2).freqlims = [0.10001, 0.2];
% % bulkbands(2).linlims = (adcpdirspectra.frequency >= bulkbands(2).freqlims(1)) & ...
% %                        (adcpdirspectra.frequency <= bulkbands(2).freqlims(2));
%
bulkbands(3).ID = 'seaswell';
bulkbands(3).freqlims = [0.06, 0.2];
% % bulkbands(3).linlims = (adcpdirspectra.frequency >= bulkbands(3).freqlims(1)) & ...
% %                        (adcpdirspectra.frequency <= bulkbands(3).freqlims(2));


%%

% % list_sources

%
list_moments = {'a1', 'b1', 'a2', 'b2'};


%%
% ---------------------------------------------
% ------- DO CALCULATIONS FOR EACH ADCP -------
% ---------------------------------------------

%% Print progress message

%
disp(' '), disp(' ')
%
disp('----- Running code for directional calculations from ADCP data. Doing calculations for moorings:  -----')
for i = 1:length(list_moorings)
    disp([num2str(i) ') ' list_moorings{i}])
end

%
totalRunTime = tic;
    
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
    dataADCP = load(fullfile(dir_dataaux, ['adcpdata_' list_moorings{i1} '.mat']));
    dataADCP = dataADCP.dataADCP;

    %
    adcpL2.mooringID = list_moorings{i1};

    %% --------------------------------------
    % Set FFT parameters based on the sampling
    % rate of the i1'th instrument
    
    %
    dtsampling = seconds(diff(dataADCP.dtime(1:2)));
    nfft = (1/dtsampling) * dt_fftsegment;

    %% --------------------------------------
    % Add fields to output data structure

    %
    adcpL2.instrument = dataADCP.instrument;
    adcpL2.latitude = dataADCP.latitude;
    adcpL2.longitude = dataADCP.longitude;
    adcpL2.site = dataADCP.site;
    adcpL2.X = dataADCP.X;
    adcpL2.Y = dataADCP.Y;

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
                    presdata_to_See(dataADCP.dtime, dataADCP.pressure, ...
                                    dataADCP.transducerHAB, ...
                                    dataADCP.bottomdepthfrompres, ...
                                    adcpL2.fftwindow, adcpL2.fftavgwindow, ...
                                    dtime_lims);
    disp('--- Done computing elevation spectra from pressure ---')
    toc(totalRunTime)
    
    %
    lbelowcutoff = (frequency_aux <= freq_cutoff);

    %
    adcpL2.frequency = frequency_aux(lbelowcutoff);
    adcpL2.frequency = adcpL2.frequency(:);

    %
    adcpL2.See = See_aux(:, lbelowcutoff);

    %
    adcpL2.bottomdepth = bottomdepth_aux(:);
    adcpL2.k = k_aux(:, lbelowcutoff);
    
    
    %% Get the frequency limits for each band for bulk statistics
    
    for i2 = 1:length(bulkbands)
        %
        bulkbands(i2).linlims = (adcpL2.frequency >= bulkbands(i2).freqlims(1)) & ...
                                (adcpL2.frequency <= bulkbands(i2).freqlims(2));
    end

    
    %% --------------------------------------
% %     % Compute a1, b2, a2, b2
% %    
% %     %
% %     adcpL2.a1 = NaN(length(adcpL2.dtime), length(adcpL2.frequency));   % ### HARD-CODED 26 ###
% %     adcpL2.b1 = adcpL2.a1;
% %     adcpL2.a2 = adcpL2.a1;
% %     adcpL2.b2 = adcpL2.a1;
% %     
% %     %
% %     disp('--- Computing first and second directional moments ---')
% %     %
% %     for i2 = 1:length(adcpL2.dtime)
% %         
% %         %
% %         linsubchunk_data = (dataADCP.dtime >= (dtimegrid(i2) - (dt_spec/2))) & ...
% %                            (dataADCP.dtime < (dtimegrid(i2) + (dt_spec/2)));
% %         
% %         % -----------------------------
% %         % Skip dtime(i2) if there is no data around the i2 grid point
% %         if ~any(linsubchunk_data)
% %             continue
% %         end
% %                        
% %         % Skip dtime(i2) if NaN is found
% %         if any(isnan(dataADCP.u(linsubchunk_data)))
% %             continue
% %         end
% %         % -----------------------------
% %     
% %                        
% %         %
% %         [freq_moments, a1, a2, b1, b2] = ...
% %                     wave_dirmoments_uvp(double(dataADCP.u(linsubchunk_data)), ...   % if data is not double, the frequency is single and approximation makes it a slightly different vector than what it should be
% %                                         double(dataADCP.v(linsubchunk_data)), ...
% %                                         double(dataADCP.pressure(linsubchunk_data)), ...
% %                                         dtsampling, ...
% %                                         nfft, nfft/2, ...
% %                                         hanning(nfft));
% %         
% %         %
% %         if i2==1
% %             lbelowcutoff = (freq_moments<=freq_cutoff);
% %         end
% %         %
% %         adcpL2.a1(i2, :) = a1(lbelowcutoff);
% %         adcpL2.b1(i2, :) = b1(lbelowcutoff);
% %         adcpL2.a2(i2, :) = a2(lbelowcutoff);
% %         adcpL2.b2(i2, :) = b2(lbelowcutoff);
% %     end
% %     disp('--- Done computing a1, b1, a2, and b2 ---')
    
    
    %% --------------------------------------
    % Compute a1, b2, a2, b2 -- similar as above
    % taking advantage that the data is gridded
    % to do the calculation efficiently
   
    %
    adcpL2.a1 = NaN(length(adcpL2.dtime), length(adcpL2.frequency));   % ### HARD-CODED 26 ###
    adcpL2.b1 = adcpL2.a1;
    adcpL2.a2 = adcpL2.a1;
    adcpL2.b2 = adcpL2.a1;
    
    %
    disp('--- Computing first and second directional moments ---')
    
    %
    [indsub, reshapeNdims, indgridlims] = reshapeintoWindows(dataADCP.dtime, dtimegrid);
    
    %
    indfill = indgridlims(1) : 1 : indgridlims(2);
    
    %
    u_array = reshape(double(dataADCP.u(indsub)), reshapeNdims);
    v_array = reshape(double(dataADCP.v(indsub)), reshapeNdims);
    p_array = reshape(double(dataADCP.pressure(indsub)), reshapeNdims);
    
    
% %     % A 2xN matrix
% %     time_bounds = [(dtimegrid(:) - (dt_spec/2)).'; ...
% %                    (dtimegrid(:) + (dt_spec/2)).'];
% %     
% %     %
% %     ind_first_wholeinterval = find(time_bounds(1, :) >= dataADCP.dtime(1), 1, 'first');
% %     ind_last_wholeinterval = find(time_bounds(2, :) < dataADCP.dtime(end), 1, 'last');
% %     
% %     %
% %     ind_data_start_firstinterval = find(dataADCP.dtime == time_bounds(1, ind_first_wholeinterval));
% %     ind_data_end_lastinterval = find(dataADCP.dtime == time_bounds(2, ind_last_wholeinterval));
% %     ind_data_end_lastinterval = ind_data_end_lastinterval - 1;
% %     %
% %     ind_allsub = ind_data_start_firstinterval : ind_data_end_lastinterval;
% %     
% %     %
% %     npts_spec = (1/dtsampling)*seconds(dt_spec);
% %     npts_intervals = length(ind_allsub)/npts_spec;

    
	% The function wave_dirmoments_uvp.m currently uses cpsd, which
    % does not accept NaNs. For now, I'll ignore the whole hour that
    % has even one NaN
    lgoodcolumns = ~any(isnan(u_array), 1) & ~any(isnan(p_array), 1);
    
    %
    tic
    [freq_moments, a1, a2, b1, b2] = ...
                wave_dirmoments_uvp(u_array(:, lgoodcolumns), ...
                                    v_array(:, lgoodcolumns), ...
                                    p_array(:, lgoodcolumns), ...   % if data is not double, the frequency is single and approximation makes it a slightly different vector than what it should be
                                    dtsampling, ...
                                    nfft, nfft/2, ...
                                    hanning(nfft));
    toc
    %
    lbelowcutoff = (freq_moments<=freq_cutoff);
    %
    adcpL2.a1(indfill(lgoodcolumns), :) = a1(lbelowcutoff, :).';
    adcpL2.b1(indfill(lgoodcolumns), :) = b1(lbelowcutoff, :).';
    adcpL2.a2(indfill(lgoodcolumns), :) = a2(lbelowcutoff, :).';
    adcpL2.b2(indfill(lgoodcolumns), :) = b2(lbelowcutoff, :).';

    disp('--- Done computing a1, b1, a2, and b2 ---')
    toc(totalRunTime)
    
   
    
    
% %     ind_match_last = find(dataADCP.dtime == dtimegrid(end));
% %     
% %     %
% %     if isempty(ind_match_first) || isempty(ind_match_last)
% %         error()
% %     end
% %     
% %     %
% %     for i2 = 1:length(adcpL2.dtime)
% %         
% %         %
% %         linsubchunk_data = (dataADCP.dtime >= (dtimegrid(i2) - (dt_spec/2))) & ...
% %                            (dataADCP.dtime < (dtimegrid(i2) + (dt_spec/2)));
% %         
% %         % -----------------------------
% %         % Skip dtime(i2) if there is no data around the i2 grid point
% %         if ~any(linsubchunk_data)
% %             continue
% %         end
% %                        
% %         % Skip dtime(i2) if NaN is found
% %         if any(isnan(dataADCP.u(linsubchunk_data)))
% %             continue
% %         end
% %         % -----------------------------
% %     
% %                        
% %         %
% %         [freq_moments, a1, a2, b1, b2] = ...
% %                     wave_dirmoments_uvp(double(dataADCP.u(linsubchunk_data)), ...   % if data is not double, the frequency is single and approximation makes it a slightly different vector than what it should be
% %                                         double(dataADCP.v(linsubchunk_data)), ...
% %                                         double(dataADCP.pressure(linsubchunk_data)), ...
% %                                         dtsampling, ...
% %                                         nfft, nfft/2, ...
% %                                         hanning(nfft));
% %         
% %         %
% %         if i2==1
% %             lbelowcutoff = (freq_moments<=freq_cutoff);
% %         end
% %         %
% %         adcpL2.a1(i2, :) = a1(lbelowcutoff);
% %         adcpL2.b1(i2, :) = b1(lbelowcutoff);
% %         adcpL2.a2(i2, :) = a2(lbelowcutoff);
% %         adcpL2.b2(i2, :) = b2(lbelowcutoff);
% %     end
% %     disp('--- Done computing a1, b1, a2, and b2 ---')
    
    %% --------------------------------------
    % Add kH, cp, and cg do data structure

	%
    H_array_aux = repmat(adcpL2.bottomdepth, 1, length(adcpL2.frequency));
    
	%
    adcpL2.kH = adcpL2.k .* H_array_aux;
    
	%
    adcpL2.cp = wave_cp(adcpL2.k, H_array_aux);
    adcpL2.cg = wave_cg(adcpL2.k, H_array_aux);
    
% %     %
% %     adcpdirspectra.E = 1025 * 9.8 * adcpdirspectra.See;
    

    %% --------------------------------------
    % Integrate flux with no direction in frequency bands

    % Flux without taking into account direction
    adcpL2.nodirection.Flux = (1025 * 9.8) * adcpL2.See .* adcpL2.cg;

    %
    for i2 = 1:length(bulkbands)
        %
        lin_aux = bulkbands(i2).linlims;
        %
        adcpL2.nodirection.(bulkbands(i2).ID).Flux = ...
                            trapz(adcpL2.frequency(lin_aux), ...
                                  adcpL2.nodirection.Flux(:, lin_aux), 2);
    end
    
    
    %% --------------------------------------
    % Compute directional spectra:
    %
    % The direction in the output is in degrees,
    % it follows the trigonometric convention,
    % and says the direction where waves are going to.

    
    dataADCP.pressure = single(dataADCP.pressure);
% %     dataADCP.bottomdepthfrompres = single(dataADCP.bottomdepthfrompres);

    %
    disp('--- Computing directional spectra ---')
    %
    directional_spec = dirspectra_ADCP_uvp(dt_fftsegment, freq_cutoff, ...
                                           dtimegrid, seconds(dt_spec), ...
                                           [dataADCP.zhab(end), dataADCP.transducerHAB], ...
                                           dataADCP.dtime, ...
                                           dataADCP.u, dataADCP.v, ...
                                           dataADCP.pressure, dataADCP.bottomdepthfrompres, ...
                                           dspec_method);
    disp('--- Done Computing directional spectra ---')
    toc(totalRunTime)
    
    % good to check that adcpL2.frequency = directional_spec.frequency,
    % which I did and it is right.

    % Add to the data structure
    for i2 = 1:length(dspec_method)
        %
        adcpL2.(dspec_method{i2}).direction = directional_spec.direction;
        %
        adcpL2.(dspec_method{i2}).See = directional_spec.(dspec_method{i2}).See;
        adcpL2.(dspec_method{i2}).D = directional_spec.(dspec_method{i2}).D;

        % Compute directional spectra of flux
        adcpL2.(dspec_method{i2}).Flux = 1025 * 9.8 * adcpL2.(dspec_method{i2}).See;

    end

    
    %% --------------------------------------
    % Spectrally compute directional properties
    % from moments, and net flux components

    %
    adcpL2.moments.meandir = atan2(adcpL2.b1, adcpL2.a1);

    %
    adcpL2.moments.dirspread = sqrt(0.5.*(1 - sqrt(adcpL2.a2.^2 + adcpL2.b2.^2)));


    %
    adcpL2.moments.Fx = 1025 * 9.8 * adcpL2.See .* adcpL2.cg .* adcpL2.a1;
    adcpL2.moments.Fy = 1025 * 9.8 * adcpL2.See .* adcpL2.cg .* adcpL2.b1;


    %% --------------------------------------
    % Compute mean direction and spread from moments
    % at different frequency bands
    
    % Then loop over bands and average a1, b1, a2, and b2.
    % Then use the standard formulas to compute mean direction
    % and directional spread
    for i2 = 1:length(bulkbands)
        %
        adcpL2.moments.(bulkbands(i2).ID).freqlims = bulkbands(i2).freqlims;
        
        %
        See_int_aux = trapz(adcpL2.frequency(bulkbands(i2).linlims), ...
                            adcpL2.See(:, bulkbands(i2).linlims), 2);

        %
        for i3 = 1:length(list_moments)
            %
            adcpL2.moments.(bulkbands(i2).ID).(list_moments{i3}) = ...
                       trapz(adcpL2.frequency(bulkbands(i2).linlims), ...
                             (adcpL2.(list_moments{i3})(:, bulkbands(i2).linlims) .* ...
                              adcpL2.See(:, bulkbands(i2).linlims)), 2);
            %
            adcpL2.moments.(bulkbands(i2).ID).(list_moments{i3}) = ...
                        (1./See_int_aux) .* adcpL2.moments.(bulkbands(i2).ID).(list_moments{i3});
        end

        % Compute mean direction and directional
        % spread for the i2'th frequency band
        %
        adcpL2.moments.(bulkbands(i2).ID).meandir = atan2(adcpL2.moments.(bulkbands(i2).ID).b1, ...
                                                          adcpL2.moments.(bulkbands(i2).ID).a1);
        %
        adcpL2.moments.(bulkbands(i2).ID).dirspread = ...
                        sqrt(0.5.*(1 - sqrt(adcpL2.moments.(bulkbands(i2).ID).a2.^2 + ...
                                            adcpL2.moments.(bulkbands(i2).ID).b2.^2)));
    end


    %% --------------------------------------
    % Compute net flux components at different
    % frequency bands

    %
    for i2 = 1:length(bulkbands)
        %
        lin_aux = bulkbands(i2).linlims;
        %
        adcpL2.moments.(bulkbands(i2).ID).Fx = trapz(adcpL2.frequency(lin_aux), adcpL2.moments.Fx(:, lin_aux), 2);
        adcpL2.moments.(bulkbands(i2).ID).Fy = trapz(adcpL2.frequency(lin_aux), adcpL2.moments.Fy(:, lin_aux), 2);
    end

    
    %%
    % ----------------------------------------------------
    % ----------------------------------------------------
    % --- NOW DO CALCULATIONS FROM DIRECIONAL SPECTRUM ---
    % ----------------------------------------------------
    % ----------------------------------------------------

    %% Select the sectors for integrating the directional
    % spectra for waves going to different directions
    %
    % Don't integrate with trapz!!! Just use sum and don't
    % worry about wrapping the angle

    %
    alloc_selectsector = false(length(adcpL2.EMEM.direction), 1);
    
    %
    sector.fullcircle.lindir = alloc_selectsector;
% %     sector.onshore.lindir = alloc_selectsector;
% %     sector.offshore.lindir = alloc_selectsector;

    %
    sector.fullcircle.lindir = ~sector.fullcircle.lindir;
    %
    sector.onshore.lindir((adcpL2.EMEM.direction >= (-90)) & (adcpL2.EMEM.direction < (90))) = true;
    %
    sector.offshore.lindir((adcpL2.EMEM.direction < (-90)) | (adcpL2.EMEM.direction >= (90))) = true;

    %
    list_sectors = fieldnames(sector);
    Nsectors = length(list_sectors);


    %% --------------------------------------
    % Compute mean direction and directional spread
    % by integrating the directional spectrum

    %
    direction_radians_array = (pi/180) * adcpL2.EMEM.direction;
    direction_radians_array = direction_radians_array.';
    %
    direction_radians_array = repmat(direction_radians_array, [length(adcpL2.frequency), 1, length(adcpL2.dtime)]);

    % Loop over directional spectrum methods
    for i2 = 1:length(dspec_method)

        % Pre-allocate space
        prealloc_aux = NaN(length(adcpL2.dtime), Nsectors);
        %
        for i3 = 1:length(bulkbands)
            %
            adcpL2.EMEM.(bulkbands(i3).ID).freqlims = bulkbands(i3).freqlims;
            %
            adcpL2.EMEM.(bulkbands(i3).ID).Hs = prealloc_aux;
            adcpL2.EMEM.(bulkbands(i3).ID).meandir = prealloc_aux;
            adcpL2.EMEM.(bulkbands(i3).ID).dirspread = prealloc_aux;
            %
            adcpL2.EMEM.(bulkbands(i3).ID).Fx = prealloc_aux;
            adcpL2.EMEM.(bulkbands(i3).ID).Fy = prealloc_aux;
        end

        %
        for i3 = 1:Nsectors
            %
            lin_aux = sector.(list_sectors{i3}).lindir;
            %
            See_cos_1_aux = (dtheta*pi/180) * squeeze(sum(cos(direction_radians_array(:, lin_aux, :)) .* adcpL2.(dspec_method(i2)).See(:, lin_aux, :), 2));    % 2 == sum over direction
            See_sin_1_aux = (dtheta*pi/180) * squeeze(sum(sin(direction_radians_array(:, lin_aux, :)) .* adcpL2.(dspec_method(i2)).See(:, lin_aux, :), 2));

            % Integrate within frequency band
            for i4 = 1:length(bulkbands)
                %
                linf_aux = bulkbands(i4).linlims;

                %
                See_int_1_aux = (dtheta*pi/180) * squeeze(sum(adcpL2.(dspec_method(i2)).See(:, lin_aux, :), 2));
                See_int_2_aux = trapz(adcpL2.frequency(linf_aux), See_int_1_aux(linf_aux, :), 1);

                % Integrate over frequency frequency
                b1bar_fullcircle_aux = trapz(adcpL2.frequency(linf_aux), See_sin_1_aux(linf_aux, :));
                a1bar_fullcircle_aux = trapz(adcpL2.frequency(linf_aux), See_cos_1_aux(linf_aux, :));
        
                % ----------------------------------------
                % Compute significant wave height
                adcpL2.(dspec_method{i2}).(bulkbands(i4).ID).Hs(:, i3) = 4*sqrt(See_int_2_aux(:));
                
                % ----------------------------------------
                % Compute mean direction
                meandir_aux = atan2(b1bar_fullcircle_aux, a1bar_fullcircle_aux);
                adcpL2.(dspec_method{i2}).(bulkbands(i4).ID).meandir(:, i3) = meandir_aux(:);
        
                % ----------------------------------------
                % Compute directional spread

                %
                meandir_array_aux = reshape(meandir_aux, [1, 1, length(adcpL2.dtime)]);
                meandir_array_aux = repmat(meandir_array_aux, [length(adcpL2.frequency), length(adcpL2.EMEM.direction), 1]);

                %
                sin2_aux = sin((direction_radians_array(:, lin_aux, :) - meandir_array_aux(:, lin_aux, :))).^2;
                %
                See_spreadweight_aux = (dtheta*pi/180) * squeeze(sum(sin2_aux .* adcpL2.(dspec_method(i2)).See(:, lin_aux, :), 2));

                %
                dirspread_aux = (1./See_int_2_aux) .* trapz(adcpL2.frequency(linf_aux), See_spreadweight_aux(linf_aux, :), 1);
                dirspread_aux = sqrt(dirspread_aux);

                %
                adcpL2.(dspec_method{i2}).(bulkbands(i4).ID).dirspread(:, i3) = dirspread_aux(:);

                % ----------------------------------------
                % Integrate flux
                adcpL2.(dspec_method{i2}).(bulkbands(i4).ID).Fx(:, i3) = 1025 * 9.8 * trapz(adcpL2.frequency(linf_aux), (See_cos_1_aux(linf_aux, :).' .* adcpL2.cg(:, linf_aux)), 2);
                adcpL2.(dspec_method{i2}).(bulkbands(i4).ID).Fy(:, i3) = 1025 * 9.8 * trapz(adcpL2.frequency(linf_aux), (See_sin_1_aux(linf_aux, :).' .* adcpL2.cg(:, linf_aux)), 2);

            end
        end
    end
    
    
    %% Save results
    %
    disp('--- Saving data structure with results ---')
    %
    save(fullfile(dir_outputdata, ...
                  ['roxsi_adcpL2_' list_moorings{i1} '.mat']), ...
         'adcpL2', '-v7.3')
    

    %% Save a version without directional spectra to see
    % how smaller the data structure is

    %
    for i2 = 1:length(dspec_method)
        %
        adcpL2.(dspec_method{i2}) = rmfield(adcpL2.(dspec_method{i2}), {'See', 'D', 'Flux'});
    end

    %
    disp('--- Saving another structure with results (no directional spectra) ---')
    %
    save(fullfile(dir_outputdata, ...
                  ['roxsi_adcpL2_nodirspec_' list_moorings{i1} '.mat']), ...
         'adcpL2', '-v7.3')
    

    %%
    disp(['--- Done with calculations for mooring ' list_moorings{i1} ' ---'])
    disp(' ')

    %
    clear adcpL2

    %
    toc(totalRunTime)
end


%% Convert mean direction and directional spread from radians to degrees

%% Convert direction dimension to whatever is physically useful

%%
disp('--- DONE WITH CALCULATIONS FOR ALL MOORINGS ---')
%
toc(totalRunTime)

%%

%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------
return

%% Pcolor of elevation spectra

%
figure
    %
    pcolor(adcpL2.dtime, adcpL2.frequency, log10(adcpL2.See.'))
    shading flat

    %
    hcb = colorbar;

    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'XGrid', 'on', 'YGrid', 'on', ...
             'YScale', 'log')


%% Compare simple frequency spectrum with that
% averaged from the integrated directional spectum

%
See_intheta = squeeze((pi/180)*sum(adcpL2.EMEM.See, 2));
See_intheta = mean(See_intheta, 2);

%
figure
    %
    plot(adcpL2.frequency, mean(adcpL2.See, 1), '.-')
    hold on
    plot(adcpL2.frequency, mean(See_intheta, 2), '.-')


%%
figure
    %
    plot(adcpL2.frequency, adcpL2.See.' ./ See_intheta, '--k')

    %
    grid on
    set(gca, 'FontSize', 16)
    %
    overlayline('h', 1, '-k')

    %
    xlabel('Frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 20)
    title('Ratio between hourly spectra', 'Interpreter', 'Latex', 'FontSize', 20)


%% Plot timeseries of mean direction

figure
    hold on
    %
    plot(adcpL2.dtime, (180/pi)*adcpL2.EMEM.swell.meandir(:, 1), '-b')
    plot(adcpL2.dtime, (180/pi)*adcpL2.moments.swell.meandir, '--b')
    %
    plot(adcpL2.dtime, (180/pi)*adcpL2.EMEM.sea.meandir(:, 1), '-r')
    plot(adcpL2.dtime, (180/pi)*adcpL2.moments.sea.meandir, '--r')

    %
    title('Mean direction ($\overline{\theta}$) at sea and swell', 'Interpreter', 'Latex', 'FontSize', 18)


%% Plot timeseries of directional spread

figure
    %
    hold on
    %
    plot(adcpL2.dtime, (180/pi)*adcpL2.EMEM.swell.dirspread(:, 1), '-b')
    plot(adcpL2.dtime, (180/pi)*adcpL2.moments.swell.dirspread, '--b')
    %
    plot(adcpL2.dtime, (180/pi)*adcpL2.EMEM.sea.dirspread(:, 1), '-r')
    plot(adcpL2.dtime, (180/pi)*adcpL2.moments.sea.dirspread, '--r')

    %
    title('Directional spread ($\sigma_{\theta}$) at sea and swell', 'Interpreter', 'Latex', 'FontSize', 18)


%% Plot fluxes - from moments (should be the net fluxes)

figure
    %
    hold on
    %
    plot(adcpL2.dtime, adcpL2.moments.swell.Fx, '-b')
    plot(adcpL2.dtime, adcpL2.moments.swell.Fy, '--b')
    %
    plot(adcpL2.dtime, adcpL2.moments.sea.Fx, '-r')
    plot(adcpL2.dtime, adcpL2.moments.sea.Fy, '--r')

    %
    grid on
    set(gca, 'FontSize', 16)
    %
    overlayline('h', 0, '-k')

    %
    title('Flux components from moments at sea and swell', 'Interpreter', 'Latex', 'FontSize', 18)


%% Cross-shore fluxes -- from moments and spectra


figure
    %
    haxs = makeSubPlots(0.1, 0.1, 0.1, ...
                        0.1, 0.1, 0.1, 1, 2);
    holdAll

        %
        plot(haxs(1), adcpL2.dtime, adcpL2.moments.swell.Fx, '--b')
        plot(haxs(1), adcpL2.dtime, adcpL2.EMEM.swell.Fx(:, 1), '-b')
        plot(haxs(1), adcpL2.dtime, adcpL2.EMEM.swell.Fx(:, 2), '-r')
        plot(haxs(1), adcpL2.dtime, adcpL2.EMEM.swell.Fx(:, 3), '-g')
    
        %
        plot(haxs(2), adcpL2.dtime, adcpL2.moments.sea.Fx, '--b')
        plot(haxs(2), adcpL2.dtime, adcpL2.EMEM.sea.Fx(:, 1), '-b')
        plot(haxs(2), adcpL2.dtime, adcpL2.EMEM.sea.Fx(:, 2), '-r')
        plot(haxs(2), adcpL2.dtime, adcpL2.EMEM.sea.Fx(:, 3), '-g')


    %
    set(haxs, 'FontSize', 16, 'Box', 'on', ...
              'XGrid', 'on', 'YGrid', 'on')
    %
    overlayline(haxs(1), 'h', 0, '-k')
    overlayline(haxs(2), 'h', 0, '-k')

    %
    title(haxs(1), ['Flux components from moments at ' ...
                    'sea and swell'], ...
                   'Interpreter', 'Latex', 'FontSize', 18)
