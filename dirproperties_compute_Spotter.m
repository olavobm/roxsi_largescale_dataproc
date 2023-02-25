%% Compute directional properties, consistent with the ADCP data

clear
close all

%
addpath(genpath(['/home/omarques/Documents/MATLAB/wafo/']))

%%

xylocalgrids = load('/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc/ROXSI_xygrids.mat');
xylocalgrids = xylocalgrids.roxsigrid;


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
% % dir_outputdata = dir_dataaux;
dir_outputdata = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results_3/data_with_dirspectra/';


%% Moorings to use for calculation

%
% list_moorings = {'B01', 'B03', 'B05'};
% list_moorings = {'B03'};
% list_moorings = {'E08'};
%
% list_moorings = {'E08', 'B01', 'B03', 'B05'};


list_moorings = {'E09', ...
                 'B01', 'B03', 'B05', ...
                 'X01', 'X03', 'X04', ...
                 'E01', 'E02', 'E05', 'E07', 'E08', 'E10', 'E11', 'E13'};



%% Define time grid for computing spectra

%
dt_spec = hours(1);

% %
% dtimegridlims = [datetime(2022, 06, 25, 12, 00, 00), ...
%                  datetime(2022, 07, 03, 12, 00, 00)];


% Full timeseries (excluding remaining data after
% recovery started, and a bit on the first day)
dtimegridlims = [datetime(2022, 06, 15, 20, 00, 00), ...
                 datetime(2022, 07, 20, 05, 00, 00)];
             
             
% % % Most of the timeseries
% % dtimegridlims = [datetime(2022, 06, 24, 18, 00, 00), ...
% %                  datetime(2022, 07, 20, 04, 00, 00)];
             
% % % Test beginning edge
% % dtimegridlims = [datetime(2022, 06, 14, 21, 00, 00), ...
% %                  datetime(2022, 06, 21, 00, 00, 00)];
% % % Test end edge
% % dtimegridlims = [datetime(2022, 07, 19, 00, 00, 00), ...
% %                  datetime(2022, 07, 21, 00, 00, 00)];

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
    spotterL1 = load(fullfile(dir_dataaux, ['spotterdata_' list_moorings{i1} '.mat']));
    spotterL1 = spotterL1.spotterL1;

    %
    spotterL2.mooringID = list_moorings{i1};

    %% --------------------------------------
    % Set FFT parameters based on the sampling
    % rate of the i1'th instrument
    
    %
    dtsampling = seconds(diff(spotterL1.displacement.dtime(1:2)));
    nfft = (1/dtsampling) * dt_fftsegment;
    
    
    %% --------------------------------------
    % Add fields to output data structure

    %
% %     spotterL2.instrument = spotterL1.instrument;
    spotterL2.latitude = spotterL1.latitude;
    spotterL2.longitude = spotterL1.longitude;
    spotterL2.site = spotterL1.site;
    spotterL2.X = spotterL1.X;
    spotterL2.Y = spotterL1.Y;

    %
    spotterL2.fftwindow = dt_fftsegment;
    spotterL2.fftavgwindow = dt_avgfft;
    %
    spotterL2.freqcutoff = freq_cutoff;

    %
    spotterL2.dt = dt_spec;
    spotterL2.dtime = dtimegrid(:);

    
    %% Rotate quantities so that they are in 
    % the local (x, y) coordinate system insteam
    % of (east, north)
    
    % 2xN matrix
    xy_array = [spotterL1.displacement.x, spotterL1.displacement.y].';
    
    %
    anglerot = xylocalgrids.ChinaRock.angleref - 270;
    %
    rot_matrix = [cosd(anglerot), -sind(anglerot); ...
                  sind(anglerot),  cosd(anglerot)];
    
	%
    xy_array = rot_matrix * xy_array;
    
    %
    spotterL1.displacement.x = xy_array(1, :).';
    spotterL1.displacement.y = xy_array(2, :).';
    
    %
    disp('----- ROTATED x, y Spotter DISPLACEMENTS INTO CHINA ROCK (x, y) !!!! -----')


    %% --------------------------------------
    % Compute vertical elevation spectra

    %
    dtime_lims = dtimegrid([1, end]);

    %
    disp('--- Computing elevation spectra ---')
    %
    [See_aux, dtimespec, frequency_aux, dof, avgpres] = ...
                spectra_scalar_reg(spotterL1.displacement.dtime, ...
                                   spotterL1.displacement.z, ...
                                   spotterL2.fftwindow, spotterL2.fftavgwindow, ...
                                   dtime_lims);
    %
    See_aux = See_aux.';    % put time in the first dimension
    disp('--- Done computing elevation spectra ---')

   
    % with dtime_lims in the input, output dtimespec should
    % be the same as the time grid
    
    %
    lbelowcutoff = (frequency_aux <= freq_cutoff);

    %
    spotterL2.frequency = frequency_aux(lbelowcutoff);
    spotterL2.frequency = spotterL2.frequency(:);

    %
    spotterL2.See = See_aux(:, lbelowcutoff);

    
    %% Get the frequency limits for each band for bulk statistics
    
    for i2 = 1:length(bulkbands)
        %
        bulkbands(i2).linlims = (spotterL2.frequency >= bulkbands(i2).freqlims(1)) & ...
                                (spotterL2.frequency <= bulkbands(i2).freqlims(2));
    end
    
    
    %% --------------------------------------
    % Compute mean water depth corresponding
    % to the timestamps of the spectra
    
    % using inefficient code...

    %
    bottomdepthavg = NaN(length(dtimespec), 1);

    %
    for i2 = 1:length(dtimespec)
        %
        linlims_aux = (spotterL1.location.dtime >= (dtimespec(i2) - seconds(spotterL2.fftavgwindow/2))) & ...
                      (spotterL1.location.dtime < (dtimespec(i2) + seconds(spotterL2.fftavgwindow/2)));
        %
        bottomdepthavg(i2) = mean(spotterL1.location.z_msl(linlims_aux), 'omitnan');

    end
    
    %
    spotterL2.bottomdepth = abs(bottomdepthavg);


    %% Compute wavenumber

    %
    freq_bounds = [(1/5000), (1/0.5)];

    %
    k_aux = wave_freqtok(spotterL2.frequency, spotterL2.bottomdepth, [], freq_bounds);
    
    %
    spotterL2.k = k_aux(:, lbelowcutoff);

    
    %% --------------------------------------
% %     % Compute a1, b1, a2, b2
% %     
% %     %
% %     spotterL2.a1 = NaN(length(spotterL2.dtime), length(spotterL2.frequency));
% %     spotterL2.b1 = spotterL2.a1;
% %     spotterL2.a2 = spotterL2.a1;
% %     spotterL2.b2 = spotterL2.a1;
% %     
% %     %
% %     disp('--- Computing first and second directional moments ---')
% %     %
% %     for i2 = 1:length(spotterL2.dtime)
% %         
% %         %
% %         linsubchunk_data = (spotterL1.displacement.dtime >= (dtimegrid(i2) - (dt_spec/2))) & ...
% %                            (spotterL1.displacement.dtime < (dtimegrid(i2) + (dt_spec/2)));
% % % % %         if i2>50
% % % % %             keyboard
% % % % %         end
% %         % -----------------------------------
% %         %
% %         if ~any(linsubchunk_data)
% %             continue
% %         end
% %         %
% %         if any(isnan(spotterL1.displacement.x(linsubchunk_data)))
% %             continue
% %         end
% %         %
% %         timediff_aux = spotterL1.displacement.dtime(linsubchunk_data);
% %         
% %         if (0.4*length(timediff_aux))<seconds(dt_spec)
% %             continue
% %         end
% %         
% %         % -----------------------------------
% %         
% %         % MINUS SIGNS ON X AND Y TO GET THE RESULTS CORRECTLY!!!
% %         % MAY WANT TO CHANGE THIS IN THE FUNCTION!!!
% %         [freq_moments, a1, a2, b1, b2] = ...
% %                     wave_dirmoments_xyz(-double(spotterL1.displacement.x(linsubchunk_data)), ...   % if data is not double, the frequency is single and approximation makes it a slightly different vector than what it should be
% %                                         -double(spotterL1.displacement.y(linsubchunk_data)), ...
% %                                         double(spotterL1.displacement.z(linsubchunk_data)), ...
% %                                         dtsampling, ...
% %                                         nfft, nfft/2, ...
% %                                         hanning(nfft));
% %         if i2==50
% %             disp('it is working!')
% %         end
% %         
% %         %
% %         if i2==1
% %             lbelowcutoff = (freq_moments<=freq_cutoff);
% %         end
% %         %
% %         spotterL2.a1(i2, :) = a1(lbelowcutoff);
% %         spotterL2.b1(i2, :) = b1(lbelowcutoff);
% %         spotterL2.a2(i2, :) = a2(lbelowcutoff);
% %         spotterL2.b2(i2, :) = b2(lbelowcutoff);
% %     end
% %     disp('--- Done computing a1, b1, a2, and b2 ---')
% %     toc(totalRunTime)
        
    %% --------------------------------------
    % Compute a1, b2, a2, b2 -- similar as above
    % taking advantage that the data is gridded
    % to do the calculation efficiently
    
    %
    spotterL2.a1 = NaN(length(spotterL2.dtime), length(spotterL2.frequency));
    spotterL2.b1 = spotterL2.a1;
    spotterL2.a2 = spotterL2.a1;
    spotterL2.b2 = spotterL2.a1;
    
    %
    disp('--- Computing first and second directional moments ---')

    %
    [indsub, reshapeNdims, indgridlims] = reshapeintoWindows(spotterL1.displacement.dtime, spotterL2.dtime);
    
    %
    indfill = indgridlims(1) : 1 : indgridlims(2);
    
    %
    x_array = reshape(spotterL1.displacement.x(indsub), reshapeNdims);
    y_array = reshape(spotterL1.displacement.y(indsub), reshapeNdims);
    z_array = reshape(spotterL1.displacement.z(indsub), reshapeNdims);


    % The function wave_dirmoments_uvp.m currently uses cpsd, which
    % does not accept NaNs. For now, I'll ignore the whole hour that
    % has even one NaN
    lgoodcolumns = ~any(isnan(x_array), 1);
    

    % MINUS SIGNS ON X AND Y TO GET THE RESULTS CORRECTLY!!!
    % MAY WANT TO CHANGE THIS IN THE FUNCTION!!!
    [freq_moments, a1, a2, b1, b2] = ...
                wave_dirmoments_xyz(-double(x_array(:, lgoodcolumns)), ...   % if data is not double, the frequency is single and approximation makes it a slightly different vector than what it should be
                                    -double(y_array(:, lgoodcolumns)), ...
                                    double(z_array(:, lgoodcolumns)), ...
                                    dtsampling, ...
                                    nfft, nfft/2, ...
                                    hanning(nfft));

    %
    lbelowcutoff = (freq_moments<=freq_cutoff);

    %
    spotterL2.a1(indfill(lgoodcolumns), :) = a1(lbelowcutoff, :).';
    spotterL2.b1(indfill(lgoodcolumns), :) = b1(lbelowcutoff, :).';
    spotterL2.a2(indfill(lgoodcolumns), :) = a2(lbelowcutoff, :).';
    spotterL2.b2(indfill(lgoodcolumns), :) = b2(lbelowcutoff, :).';
    %
    disp('--- Done computing a1, b1, a2, and b2 ---')
    toc(totalRunTime)
    

    %%


% %     %
% %     for i2 = 1:length(spotterL2.dtime)
% %         
% %         %
% % % %         linsubchunk_data = (spotterL1.displacement.dtime >= (dtimegrid(i2) - (dt_spec/2))) & ...
% % % %                            (spotterL1.displacement.dtime < (dtimegrid(i2) + (dt_spec/2)));
% % % % %         if i2>50
% % % % %             keyboard
% % % % %         end
% %         % -----------------------------------
% %         %
% %         if ~any(linsubchunk_data)
% %             continue
% %         end
% %         %
% %         if any(isnan(spotterL1.displacement.x(linsubchunk_data)))
% %             continue
% %         end
% %         %
% %         timediff_aux = spotterL1.displacement.dtime(linsubchunk_data);
% %         
% %         if (0.4*length(timediff_aux))<seconds(dt_spec) %????
% %             continue
% %         end
% %         
% % 
% %     end
% % % %     disp('--- Done computing a1, b1, a2, and b2 ---')
% % % %     toc(totalRunTime)


    %% --------------------------------------
    % Add kH, cp, and cg do data structure
    
	%
    H_array_aux = repmat(spotterL2.bottomdepth, 1, length(spotterL2.frequency));
    
	%
    spotterL2.kH = spotterL2.k .* H_array_aux;
    
	%
    spotterL2.cp = wave_cp(spotterL2.k, H_array_aux);
    spotterL2.cg = wave_cg(spotterL2.k, H_array_aux);
    
    

    %% --------------------------------------
    % Integrate flux with no direction in frequency bands

    % Flux without taking into account direction
    spotterL2.nodirection.Flux = (1025 * 9.8) * spotterL2.See .* spotterL2.cg;

    %
    for i2 = 1:length(bulkbands)
        %
        lin_aux = bulkbands(i2).linlims;
        %
        spotterL2.nodirection.(bulkbands(i2).ID).Flux = ...
                            trapz(spotterL2.frequency(lin_aux), ...
                                  spotterL2.nodirection.Flux(:, lin_aux), 2);
    end
    
    
    %% --------------------------------------
    % Compute directional spectra:
    %
    % The direction in the output is in degrees,
    % it follows the trigonometric convention,
    % and says the direction where waves are going to.

    %
    spotterL1.displacement.x = double(spotterL1.displacement.x);
    spotterL1.displacement.y = double(spotterL1.displacement.y);  % double for x and y ???
    spotterL1.displacement.z = single(spotterL1.displacement.z);

    % Although only the hourly averaged bottom depth is used,
    % interpolate to the displacement time grid to run the
    % function in a simple way
    bottomdepth_aux = interp1(spotterL1.location.dtime, ...
                              spotterL1.location.z_msl, ...
                              spotterL1.displacement.dtime);
	bottomdepth_aux = abs(bottomdepth_aux);

    %
    disp('--- Computing directional spectra ---')
    %
    directional_spec = dirspectra_xyz(dt_fftsegment, freq_cutoff, ...
                                           dtimegrid, seconds(dt_spec), ...
                                           spotterL1.displacement.dtime, ...
                                           spotterL1.displacement.x, ...
                                           spotterL1.displacement.y, ...
                                           spotterL1.displacement.z, ...
                                           bottomdepth_aux, ...
                                           dspec_method);
    disp('--- Done Computing directional spectra ---')
    toc(totalRunTime)
        
    % good to check that adcpL2.frequency = directional_spec.frequency,
    % which I did and it is right.

    % Add to the data structure
    for i2 = 1:length(dspec_method)
        %
        spotterL2.(dspec_method{i2}).direction = directional_spec.direction;
        %
        spotterL2.(dspec_method{i2}).See = directional_spec.(dspec_method{i2}).See;
        spotterL2.(dspec_method{i2}).D = directional_spec.(dspec_method{i2}).D;

        % Compute directional spectra of flux
        spotterL2.(dspec_method{i2}).Flux = 1025 * 9.8 * spotterL2.(dspec_method{i2}).See;

    end
    

    %% --------------------------------------
    % Spectrally compute directional properties
    % from moments, and net flux components

    %
    spotterL2.moments.meandir = atan2(spotterL2.b1, spotterL2.a1);

    %
    spotterL2.moments.dirspread = sqrt(0.5.*(1 - sqrt(spotterL2.a2.^2 + spotterL2.b2.^2)));


    %
    spotterL2.moments.Fx = 1025 * 9.8 * spotterL2.See .* spotterL2.cg .* spotterL2.a1;
    spotterL2.moments.Fy = 1025 * 9.8 * spotterL2.See .* spotterL2.cg .* spotterL2.b1;


    %% --------------------------------------
    % Compute mean direction and spread from moments
    % at different frequency bands
    
    % Then loop over bands and average a1, b1, a2, and b2.
    % Then use the standard formulas to compute mean direction
    % and directional spread
    for i2 = 1:length(bulkbands)
        %
        spotterL2.moments.(bulkbands(i2).ID).freqlims = bulkbands(i2).freqlims;
        
        %
        See_int_aux = trapz(spotterL2.frequency(bulkbands(i2).linlims), ...
                            spotterL2.See(:, bulkbands(i2).linlims), 2);

        %
        for i3 = 1:length(list_moments)
            %
            spotterL2.moments.(bulkbands(i2).ID).(list_moments{i3}) = ...
                       trapz(spotterL2.frequency(bulkbands(i2).linlims), ...
                             (spotterL2.(list_moments{i3})(:, bulkbands(i2).linlims) .* ...
                              spotterL2.See(:, bulkbands(i2).linlims)), 2);
            %
            spotterL2.moments.(bulkbands(i2).ID).(list_moments{i3}) = ...
                        (1./See_int_aux) .* spotterL2.moments.(bulkbands(i2).ID).(list_moments{i3});
        end

        % Compute mean direction and directional
        % spread for the i2'th frequency band
        %
        spotterL2.moments.(bulkbands(i2).ID).meandir = atan2(spotterL2.moments.(bulkbands(i2).ID).b1, ...
                                                             spotterL2.moments.(bulkbands(i2).ID).a1);
        %
        spotterL2.moments.(bulkbands(i2).ID).dirspread = ...
                        sqrt(0.5.*(1 - sqrt(spotterL2.moments.(bulkbands(i2).ID).a2.^2 + ...
                                            spotterL2.moments.(bulkbands(i2).ID).b2.^2)));
    end

    %% --------------------------------------
    % Compute net flux components at different
    % frequency bands

    %
    for i2 = 1:length(bulkbands)
        %
        lin_aux = bulkbands(i2).linlims;
        %
        spotterL2.moments.(bulkbands(i2).ID).Fx = trapz(spotterL2.frequency(lin_aux), spotterL2.moments.Fx(:, lin_aux), 2);
        spotterL2.moments.(bulkbands(i2).ID).Fy = trapz(spotterL2.frequency(lin_aux), spotterL2.moments.Fy(:, lin_aux), 2);
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
    alloc_selectsector = false(length(spotterL2.EMEM.direction), 1);
    
    %
    sector.fullcircle.lindir = alloc_selectsector;
% %     sector.onshore.lindir = alloc_selectsector;
% %     sector.offshore.lindir = alloc_selectsector;

    %
    sector.fullcircle.lindir = ~sector.fullcircle.lindir;
    %
    sector.onshore.lindir((spotterL2.EMEM.direction >= (-90)) & (spotterL2.EMEM.direction < (90))) = true;
    %
    sector.offshore.lindir((spotterL2.EMEM.direction < (-90)) | (spotterL2.EMEM.direction >= (90))) = true;

    %
    list_sectors = fieldnames(sector);
    Nsectors = length(list_sectors);


    %% --------------------------------------
    % Compute mean direction and directional spread
    % by integrating the directional spectrum

    %
    direction_radians_array = (pi/180) * spotterL2.EMEM.direction;
    direction_radians_array = direction_radians_array.';
    %
    direction_radians_array = repmat(direction_radians_array, [length(spotterL2.frequency), 1, length(spotterL2.dtime)]);

    % Loop over directional spectrum methods
    for i2 = 1:length(dspec_method)

        % Pre-allocate space
        prealloc_aux = NaN(length(spotterL2.dtime), Nsectors);
        %
        for i3 = 1:length(bulkbands)
            %
            spotterL2.EMEM.(bulkbands(i3).ID).freqlims = bulkbands(i3).freqlims;
            %
            spotterL2.EMEM.(bulkbands(i3).ID).Hs = prealloc_aux;
            spotterL2.EMEM.(bulkbands(i3).ID).meandir = prealloc_aux;
            spotterL2.EMEM.(bulkbands(i3).ID).dirspread = prealloc_aux;
            %
            spotterL2.EMEM.(bulkbands(i3).ID).Fx = prealloc_aux;
            spotterL2.EMEM.(bulkbands(i3).ID).Fy = prealloc_aux;
        end

        %
        for i3 = 1:Nsectors
            %
            lin_aux = sector.(list_sectors{i3}).lindir;
            %
            See_cos_1_aux = (dtheta*pi/180) * squeeze(sum(cos(direction_radians_array(:, lin_aux, :)) .* spotterL2.(dspec_method(i2)).See(:, lin_aux, :), 2));    % 2 == sum over direction
            See_sin_1_aux = (dtheta*pi/180) * squeeze(sum(sin(direction_radians_array(:, lin_aux, :)) .* spotterL2.(dspec_method(i2)).See(:, lin_aux, :), 2));

            % Integrate within frequency band
            for i4 = 1:length(bulkbands)
                %
                linf_aux = bulkbands(i4).linlims;

                %
                See_int_1_aux = (dtheta*pi/180) * squeeze(sum(spotterL2.(dspec_method(i2)).See(:, lin_aux, :), 2));
                See_int_2_aux = trapz(spotterL2.frequency(linf_aux), See_int_1_aux(linf_aux, :), 1);

                % Integrate over frequency frequency
                b1bar_fullcircle_aux = trapz(spotterL2.frequency(linf_aux), See_sin_1_aux(linf_aux, :));
                a1bar_fullcircle_aux = trapz(spotterL2.frequency(linf_aux), See_cos_1_aux(linf_aux, :));
        
                % ----------------------------------------
                % Compute significant wave height
                spotterL2.(dspec_method{i2}).(bulkbands(i4).ID).Hs(:, i3) = 4*sqrt(See_int_2_aux(:));
                
                % ----------------------------------------
                % Compute mean direction
                meandir_aux = atan2(b1bar_fullcircle_aux, a1bar_fullcircle_aux);
                spotterL2.(dspec_method{i2}).(bulkbands(i4).ID).meandir(:, i3) = meandir_aux(:);
        
                % ----------------------------------------
                % Compute directional spread

                %
                meandir_array_aux = reshape(meandir_aux, [1, 1, length(spotterL2.dtime)]);
                meandir_array_aux = repmat(meandir_array_aux, [length(spotterL2.frequency), length(spotterL2.EMEM.direction), 1]);

                %
                sin2_aux = sin((direction_radians_array(:, lin_aux, :) - meandir_array_aux(:, lin_aux, :))).^2;
                %
                See_spreadweight_aux = (dtheta*pi/180) * squeeze(sum(sin2_aux .* spotterL2.(dspec_method(i2)).See(:, lin_aux, :), 2));

                %
                dirspread_aux = (1./See_int_2_aux) .* trapz(spotterL2.frequency(linf_aux), See_spreadweight_aux(linf_aux, :), 1);
                dirspread_aux = sqrt(dirspread_aux);

                %
                spotterL2.(dspec_method{i2}).(bulkbands(i4).ID).dirspread(:, i3) = dirspread_aux(:);

                % ----------------------------------------
                % Integrate flux
                spotterL2.(dspec_method{i2}).(bulkbands(i4).ID).Fx(:, i3) = 1025 * 9.8 * trapz(spotterL2.frequency(linf_aux), (See_cos_1_aux(linf_aux, :).' .* spotterL2.cg(:, linf_aux)), 2);
                spotterL2.(dspec_method{i2}).(bulkbands(i4).ID).Fy(:, i3) = 1025 * 9.8 * trapz(spotterL2.frequency(linf_aux), (See_sin_1_aux(linf_aux, :).' .* spotterL2.cg(:, linf_aux)), 2);

            end
        end
    end


    %% Save results
    %
    disp('--- Saving data structure with results ---')
    %
    save(fullfile(dir_outputdata, ...
                  ['roxsi_spotterL2_' list_moorings{i1} '.mat']), ...
         'spotterL2', '-v7.3')
    

    %% Save a version without directional spectra to see
    % how smaller the data structure is

    %
    for i2 = 1:length(dspec_method)
        %
        spotterL2.(dspec_method{i2}) = rmfield(spotterL2.(dspec_method{i2}), {'See', 'D', 'Flux'});
    end

    %
    disp('--- Saving another structure with results (no directional spectra) ---')
    %
    save(fullfile(dir_outputdata, ...
                  ['roxsi_spotterL2_nodirspec_' list_moorings{i1} '.mat']), ...
         'spotterL2', '-v7.3')


    %%
    disp(['--- Done with calculations for mooring ' list_moorings{i1} ' ---'])
    disp(' ')

    %
    clear spotterL2

    %
    toc(totalRunTime)
end


%%
disp('--- DONE WITH CALCULATIONS FOR ALL MOORINGS ---')
%
toc(totalRunTime)
    