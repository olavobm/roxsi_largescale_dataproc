%% Script that process Aquadopp data from RAW to level 1.

clear
close all

%%
% --------------------------------------
% --------- PRELIMINARY STUFF ----------
% --------------------------------------

%%

%
% % dirparent_data = data_dirpath();
dirparent_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
%
dir_rawdata_parent = fullfile(dirparent_data, 'RAW', 'Aquadopp');


%%

%
% % dir_output_parent = data_dirpath();
% % dir_output_parent = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/';
% dir_output_parent = '/home/omarques/Documents/MATLAB/ROXSIproc_output/';
% % %
% % %
% % dir_output_data_L1 = pwd;
% % dir_output_figs_L1 = pwd;

%
dir_output_parent = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/';
dir_output_data_L1 = fullfile(dir_output_parent, 'Level1_Data', 'Aquadopp_Level1');
dir_output_figs_L1 = fullfile(dir_output_data_L1, 'figs_QC');

% Logical switches to save or not save data and figures
lsave_file = true;
lsave_fig = true;


%% Load ADCP deployment information

% Just to be clear: file and variable have
% the same name (though not a requirement)
% dir_coderepo = repo_dirpath();
dir_coderepo = '/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc';
%
load(fullfile(dir_coderepo, 'deploymentInfo_ROXSI2022.mat'), 'deploymentInfo_ROXSI2022')


%% List of Aquadopps that will be processed

% All Aquadopps
list_Aquadopp = {'A03_5380', ...
                 'B02_12507', 'B04_2147', 'B07_2141', 'B08_13288', 'B11_12280', ...
                 'C03_0709', ...
                 'D01_12346', 'D02_0653', ...
                 'E03_13300', 'E04_13172', 'E06_9736', 'E12_11150', ...
                 'F01_9995', 'F02_5838', 'F03_5384', 'F04_5401', 'F05_14032', ...
                 'X06_13290', 'X13_9945'};

% Just a test
% list_Aquadopp = {'B02_12507', 'X13_9945'};
% list_Aquadopp = {'X13_9945'};
% list_Aquadopp = {'F03_5384'};
% list_Aquadopp = {'A03_5380'};

% % % A subset with a few
% % list_Aquadopp = {'A03_5380', ...
% %                  'B02_12507', 'B11_12280', ...
% %                  'C03_0709', ...
% %                  'E03_13300', ...
% %                  'F01_9995', 'F02_5838', 'F03_5384'};

%
% % list_Aquadopp = {'B02_12507', 'E03_13300'};

%
Naquadopps = length(list_Aquadopp);


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------------- DO DATA PROCESSING ----------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------



%% Load atmospheric pressure

%
atm_pressure = load(fullfile(dirparent_data, 'RAW', ...
                      'noaa_mry_barometric_pressure', 'atm_pressure.mat'));
% Make correction for the spatial variability
atm_pressure.atm_pres = atm_pressure.atm_pres - (100*0.032);    % the correction is 0.032 dbar, and atmospheric pressure is in milibar

% Time limits when pressure offset was set for Aquadopps
time_lims = [datetime(2022, 06, 13, 08, 0, 0), ...
             datetime(2022, 06, 14, 15, 0, 0)];
time_lims.TimeZone = atm_pressure.time_vec.TimeZone;
%
lintimelims = (atm_pressure.time_vec >= time_lims(1)) & ...
              (atm_pressure.time_vec <= time_lims(2));

%
atmpresanomaly.dtime = atm_pressure.time_vec;
atmpresanomaly.atm_anomaly = (atm_pressure.atm_pres - mean(atm_pressure.atm_pres(lintimelims)))./100;
atmpresanomaly.units = 'dbar';

% % % % Check the atmospheric pressure correction
% % %
% % figure
% %     plot(atm_pressure.time_vec, atm_pressure.atm_pres)
% %     hold on
% %     overlayline('v', time_lims, '--k')
% % 
% % %
% % figure
% %     plot(atmpresanomaly.dtime, atmpresanomaly.atm_anomaly)
% %     hold on
% %     overlayline('v', time_lims, '--k')


%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_Aquadopp_procL1_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_output_data_L1, log_file_name))

%
totalRunTime = tic;

%% Display on the screen:

%
disp(' '), disp(' ')
disp('------------------------------ Processing data from Aquadopps ------------------------------')
disp('List of Aquadopps being processed:')
%
for i = 1:Naquadopps
    disp([num2str(i) ' - ' list_Aquadopp{i}])
end



%% Process RAW Aquadopp data

%
list_senfile_vars = ["time", "heading", "pitch", "roll", "pressure", "temperature"];

%
tic
% Loop over Aquadopps in the list
for i = 1:Naquadopps

    % ----------------------------------------------------
    %
    disp(' '), disp(' ')
    disp(['----- Start processing raw Aquadopp data: ' list_Aquadopp{i} ' -----'])

    % ----------------------------------------------------
    % Get header
    file_header_aux = dir(fullfile(dir_rawdata_parent, list_Aquadopp{i}, '*.hdr'));

    % E12 has extra files (E1202.hdr, E1202.PRF, E1202.ssl). From header
    % file, it seems these belong to this deployment, but these were not
    % used (some kind of error?). But all the data seems to be in E1203.*
    if strcmp(list_Aquadopp{i}, 'E12_11150')
        file_header_aux = dir(fullfile(dir_rawdata_parent, list_Aquadopp{i}, 'E1203.hdr'));
    end

    %
    disp('----- Loading header file -----')

    %
    header_aux = Aquadopp_read_header(fullfile(file_header_aux.folder, file_header_aux.name));

    % ----------------------------------------------------
    %
    file_sen_aux = dir(fullfile(dir_rawdata_parent, list_Aquadopp{i}, '*.sen'));

    %
    senAQDP_aux = Aquadopp_read_senfile(fullfile(file_sen_aux.folder, file_sen_aux.name), list_senfile_vars);

    %
    disp('----- Correcting for clock drift -----')

    % Now correct for clock drift
    time_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, ...
                                            list_Aquadopp{i}(5:end), ...
                                            datenum(senAQDP_aux.time, "yyyy/mm/dd HH:MM:SS"));
    %
    senAQDP_aux.time = datetime(time_aux, ...
                                    'ConvertFrom', 'datenum', ...
                                    'TimeZone', 'America/Los_Angeles');

    % ----------------------------------------------------
    % Load velocity and backscatter data

    disp('----- Loading velocity data -----')

    %
    file_noextension_aux = file_sen_aux.name(1:(strfind(file_sen_aux.name, '.') - 1));

    %
    beamAQDP_aux = Aquadopp_read_beamdata(fullfile(dir_rawdata_parent, list_Aquadopp{i}, file_noextension_aux), ...
                                          ["a1", "a2", "a3", "v1", "v2", "v3"]);


    % ----------------------------------------------------
    % Trim for when instrument was in the water

    %
    ind_row_match = find(strcmp(deploymentInfo_ROXSI2022.SN, list_Aquadopp{i}(5:end)));
    
    %
    time_1 = deploymentInfo_ROXSI2022.time_begin_trim(ind_row_match);
    time_2 = deploymentInfo_ROXSI2022.time_end_trim(ind_row_match);

    %
    time_1 = datetime(datenum(time_1, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', senAQDP_aux.time.TimeZone);
    time_2 = datetime(datenum(time_2, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', senAQDP_aux.time.TimeZone);

    %
    lin_deployment = (senAQDP_aux.time >= time_1) & (senAQDP_aux.time <= time_2);

    %
    disp(['----- Trimmed data between deployment and recovery: between ' datestr(time_1) ' and ' datestr(time_2) ' -----'])


    % ----------------------------------------------------
    % Start adding data into level 1 data structure
    
    %
    lmatch_ontable = strcmp(deploymentInfo_ROXSI2022.SN, list_Aquadopp{i}(5:end));

    %
    aquadoppL1.SN = convertCharsToStrings(list_Aquadopp{i}(5:end));
    aquadoppL1.mooringID = deploymentInfo_ROXSI2022.mooringID(lmatch_ontable);

    % Latitude/longitude
    info_mooringtable = ROXSI_mooringlocation(aquadoppL1.mooringID, "ADCP");
    %
    aquadoppL1.latitude = info_mooringtable.latitude;
    aquadoppL1.longitude = info_mooringtable.longitude;
    %
    aquadoppL1.site = info_mooringtable.roxsiarray;

    % Get (x, y) location in local coordinate system
    [aquadoppL1.X, ...
     aquadoppL1.Y] = ROXSI_lltoxy(aquadoppL1.latitude, ...
                                  aquadoppL1.longitude, aquadoppL1.site);

    % ----------------------------------------------------
    % Metadata of the ADCP settings

    %
    aquadoppL1.header = header_aux.header;

    % Both of these are in seconds
    aquadoppL1.samplingtime = header_aux.profileInterval;
    aquadoppL1.averagingtime = header_aux.averageInterval;

    % In meters
    aquadoppL1.binsize = header_aux.binsize;

    % In meters
    aquadoppL1.cellcenter = header_aux.cellcenter;

    % Height of the transducers above the bottom (in meters)
    if strcmp(aquadoppL1.mooringID, "B02at")
        %
% %         aquadoppL1.transducerHAB = 1;   % PLACEHOLDER VALUE
        %
        aquadoppL1.transducerHAB = (80)/100;
        % approximate value (taking into account that the Aquadopp was
        % positioned such that the end of the dummy plug was just above the
        % ground, when the SeaSpider did not have the lead feet). The
        % height of the lead feet is 5.5 cm and the length between
        % transducer and end of the ADCP (without dummy plug) is about
        % 63 cm).
    %
    else
        % Based on top and bottom of the transducers
        % from the ADCP mount schematics on Solidworks:
        aquadoppL1.transducerHAB = ((12.44 + 10.53)/2)/100;
    end

    % ----------------------------------------------------
    % Create (time-constant) height above the bottom vector
    aquadoppL1.zhab = aquadoppL1.transducerHAB + aquadoppL1.cellcenter;

    % ----------------------------------------------------
    % Put all the raw data in a single structure that will be saved

    %
% % %     aquadoppL1.timezone = senAQDP_aux.time.TimeZone;
    %
    aquadoppL1.dtime = senAQDP_aux.time(lin_deployment).';
    aquadoppL1.heading = senAQDP_aux.heading(lin_deployment).';
    aquadoppL1.pitch = senAQDP_aux.pitch(lin_deployment).';
    aquadoppL1.roll = senAQDP_aux.roll(lin_deployment).';

    %
    aquadoppL1.pressureoffset = deploymentInfo_ROXSI2022.pressure_offset(ind_row_match);
    
    % Applies pressure offset set when programming
    % the Aquadopp (should be 1, but there are a
    % couple that seem to have been set to 0).
    aquadoppL1.pressure = senAQDP_aux.pressure(lin_deployment).' - aquadoppL1.pressureoffset;

    % ------------------------------
    % Correct for the atmospheric pressure variability

    %
    disp('----- Making small correction to pressure due to atmospheric pressure variability -----')

    %
    atm_anomaly_interp = interp1(atmpresanomaly.dtime, ...
                                 atmpresanomaly.atm_anomaly, ...
                                 aquadoppL1.dtime);
    
    %
    aquadoppL1.pressure = aquadoppL1.pressure - atm_anomaly_interp;

    % ------------------------------

    %
    aquadoppL1.temperature = senAQDP_aux.temperature(lin_deployment).';

    % ----------------------------------------------------
    % Trim the data in the vertical using the bottom depth estimate

    %
    disp('----- Computing SSH from pressure and NaNing data above the surface -----')

% % 
% %     % One very simple and quick selection
% %     aquadoppL1.ptrimTH = max(aquadoppL1.pressure) + (2 * aquadoppL1.binsize);

    % Compute depth from hydrostatic pressure 
    aquadoppL1.bottomdepthfrompres = 1e4*aquadoppL1.pressure ./ (1030*9.8);

    % Find all bins that are below the maximum bottom depth
    lin_verticalrange = ((aquadoppL1.zhab + (aquadoppL1.binsize/2)) < ...
                         max(aquadoppL1.bottomdepthfrompres));
    %
    aquadoppL1.zhab = aquadoppL1.zhab(lin_verticalrange);

    % ----------------------------------------------------
    %
%     aquadoppL1.coordsystem = 'ENU';
    aquadoppL1.coordsystem = 'localXY';

    % In clockwise degrees from the true north
    aquadoppL1.magdec = 12.86;

    % At this point, Ue and Vn are relative to the MAGNETIC north
    aquadoppL1.Ue = beamAQDP_aux.v1(lin_verticalrange, lin_deployment);
    aquadoppL1.Vn = beamAQDP_aux.v2(lin_verticalrange, lin_deployment);
    aquadoppL1.Wup = beamAQDP_aux.v3(lin_verticalrange, lin_deployment);
    %
    aquadoppL1.a1 = beamAQDP_aux.a1(lin_verticalrange, lin_deployment);
    aquadoppL1.a2 = beamAQDP_aux.a2(lin_verticalrange, lin_deployment);
    aquadoppL1.a3 = beamAQDP_aux.a3(lin_verticalrange, lin_deployment);

    % ----------------------------------------------------
    % NaN values above the surface at every time
    %
    zhab_halfstep = aquadoppL1.zhab + (aquadoppL1.binsize/2);
    Nbins = length(aquadoppL1.zhab);
    %
    ind_abovesurface = 1:(size(aquadoppL1.Ue, 1) * size(aquadoppL1.Ue, 2));
    ind_abovesurface = reshape(ind_abovesurface, size(aquadoppL1.Ue, 1), size(aquadoppL1.Ue, 2));
    %
    for i2 = 1:length(aquadoppL1.pressure)
        %
        ind_above_aux = find((zhab_halfstep < aquadoppL1.bottomdepthfrompres(i2)), 1, 'last');
        %
        ind_abovesurface(1:ind_above_aux, i2) = NaN;
    end
    %
    ind_abovesurface = ind_abovesurface(:);
    ind_abovesurface = ind_abovesurface(~isnan(ind_abovesurface));
    %
    aquadoppL1.Ue(ind_abovesurface) = NaN;
    aquadoppL1.Vn(ind_abovesurface) = NaN;
    aquadoppL1.Wup(ind_abovesurface) = NaN;
    %
    aquadoppL1.a1(ind_abovesurface) = NaN;
    aquadoppL1.a2(ind_abovesurface) = NaN;
    aquadoppL1.a3(ind_abovesurface) = NaN;


% %     % ----------------------------------------------------
% %     % Rotate horizontal velocity from magnetic to true north
% % 
% %     %
% %     disp(['----- Rotating horizontal velocity to ' ...
% %           aquadoppL1.coordsystem ', with mag. declination of ' ...
% %           num2str(aquadoppL1.magdec, '%.2f') ' degrees -----'])
% % 
% %     %
% %     rotMatrix = [cosd(aquadoppL1.magdec), sind(aquadoppL1.magdec); ...
% %                  -sind(aquadoppL1.magdec), cosd(aquadoppL1.magdec)];
% % 
% % % %     % Check the rotation (i.e. velocity aligned with magnetic north
% % % %     % should have a small zonal component and large meridional
% % % %     % component in a geographical north coordinate system)
% % % %     rotMatrix * [0; 1]
% % 
% %     %
% %     for i2 = 1:size(aquadoppL1.Ue, 1)
% %         %
% %         uv_aux = [aquadoppL1.Ue(i2, :); aquadoppL1.Vn(i2, :)];
% %         %
% %         uv_rot_aux = rotMatrix * uv_aux;
% % 
% %         %
% %         aquadoppL1.Ue(i2, :) = uv_aux(1, :);
% %         aquadoppL1.Vn(i2, :) = uv_aux(2, :);
% %     end


    % ----------------------------------------------------
    % Rotate horizontal velocity from magnetic north to local X/Y

    %
    disp(['--- Rotating horizontal velocity to ' ...
          aquadoppL1.coordsystem ', with mag. declination of ' ...
          num2str(aquadoppL1.magdec, '%.2f') ' degrees ---'])

    %
    [aquadoppL1.u, aquadoppL1.v] = ROXSI_uv_ENtoXY(aquadoppL1.Ue, aquadoppL1.Vn, aquadoppL1.site, true);

    % Rename w
    aquadoppL1.w = aquadoppL1.Wup;


    %%
    % ----------------------------------------------------
    % Time-grid all the variables to the resolution of the sampling rate

    %
    disp('--- Gridding data in (date)time ---')

    % As time edges of the grid, round to the whole minute (done in a more
    % complicated way so that it doesn't rely on Matlab versions newer than
    % at least 2021b)
    dtime_edge_1 = datetime(aquadoppL1.dtime(1).Year, aquadoppL1.dtime(1).Month, aquadoppL1.dtime(1).Day, ...
                            aquadoppL1.dtime(1).Hour, (aquadoppL1.dtime(1).Minute + 1), 00);
    dtime_edge_2 = datetime(aquadoppL1.dtime(end).Year, aquadoppL1.dtime(end).Month, aquadoppL1.dtime(end).Day, ...
                            aquadoppL1.dtime(end).Hour, aquadoppL1.dtime(end).Minute, 00);
    % In seconds
    dt_sampling = aquadoppL1.samplingtime;

    %
    dtime_grid = dtime_edge_1 : seconds(dt_sampling) : dtime_edge_2;
    dtime_grid.TimeZone = aquadoppL1.dtime.TimeZone;

    %
    Nlengthtimeseries = length(aquadoppL1.dtime);
    list_time_vars = {'dtime'};
    %
    list_variables_aux = fieldnames(aquadoppL1);

    %
    for i2 = 1:length(list_variables_aux)

        % Only interpolate fields that are NOT contained in list_time_vars
        % and are vectors of the correct length
        if ~any(contains(list_time_vars, list_variables_aux{i2})) && ...
           isvector(aquadoppL1.(list_variables_aux{i2})) && ...
           (length(aquadoppL1.(list_variables_aux{i2})) == Nlengthtimeseries)
           
            %
            aquadoppL1.(list_variables_aux{i2}) = ...
                            interp1(aquadoppL1.dtime, ...              
                                    aquadoppL1.(list_variables_aux{i2}), ...
                                    dtime_grid);

            % Turn to column vector
            aquadoppL1.(list_variables_aux{i2}) = aquadoppL1.(list_variables_aux{i2})(:);

        end

        % Only interpolate fields that are NOT contained in list_time_vars
        % and are matrices with correct number of columns
        if ~any(contains(list_time_vars, list_variables_aux{i2})) && ...
           ~isvector(aquadoppL1.(list_variables_aux{i2})) && ...
           (size(aquadoppL1.(list_variables_aux{i2}), 2) == Nlengthtimeseries)
           
            %
            var_aux = NaN(size(aquadoppL1.(list_variables_aux{i2}), 1), length(dtime_grid));

            %
            for i3 = 1:size(aquadoppL1.(list_variables_aux{i2}), 1)
                %
                var_aux(i3, :) = interp1(aquadoppL1.dtime, ...
                                         aquadoppL1.(list_variables_aux{i2})(i3, :), ...
                                         dtime_grid);
            end

            % Replace variable in data structure
            aquadoppL1.(list_variables_aux{i2}) = var_aux;

        end
    end

    % Replace measured time stamps by time grid
    aquadoppL1.dtime = dtime_grid(:);

    %
    disp('--- Done with time gridding ---')
    
    % ----------------------------------------------------
    % Filter out velocity where amplitude is below a threshold value

    %
    disp('----- Removing velocity where backscatter is smaller than threshold value -----')

    %
    aquadoppL1.backscatterTH = 30;

% %     %
% %     l_belowTH = (aquadoppL1.a1 <= aquadoppL1.backscatterTH) | ...
% %                 (aquadoppL1.a2 <= aquadoppL1.backscatterTH) | ...
% %                 (aquadoppL1.a3 <= aquadoppL1.backscatterTH);
% %     %
% %     aquadoppL1.Ue(l_belowTH) = NaN;
% %     aquadoppL1.Vn(l_belowTH) = NaN;
% %     aquadoppL1.Wup(l_belowTH) = NaN;

    %
    aquadoppL1.la1belowTH = (aquadoppL1.a1 <= aquadoppL1.backscatterTH);
    aquadoppL1.la2belowTH = (aquadoppL1.a2 <= aquadoppL1.backscatterTH);
    aquadoppL1.la3belowTH = (aquadoppL1.a3 <= aquadoppL1.backscatterTH);
    %
    aquadoppL1.lanyabelowTH = (aquadoppL1.la1belowTH | ...
                               aquadoppL1.la2belowTH | ...
                               aquadoppL1.la3belowTH);
                

    % ----------------------------------------------------
    % Turn all vectors into column vectors so that Matlab
    % can quickly displace the structure variable in the 
    % command window (Matlab displays first elements of
    % row vectors, and it takes longer)
    list_fields_aux = fieldnames(aquadoppL1);
    %
    for i2 = 1:length(list_fields_aux)
        %
        if isvector(aquadoppL1.(list_fields_aux{i2})) && ...
           ~isstruct(aquadoppL1.(list_fields_aux{i2})) && ...
           ~ischar(aquadoppL1.(list_fields_aux{i2}))
            %
            aquadoppL1.(list_fields_aux{i2}) = aquadoppL1.(list_fields_aux{i2})(:);
        end
    end


    % ----------------------------------------------------
    % Compute averaged fields (compute 10 or 20 min averages
    % because B02 gives data every 10 min, so only multiples
    % of 10 min will work for all ADCPs).

    %
    disp('----- Computing averaged fields -----')

    %
    aquadoppL1.averaged.dt = 10*60;

    %
    npts_avg = aquadoppL1.averaged.dt / aquadoppL1.samplingtime;

    % First average time
    aquadoppL1.averaged.dtime = regRunMean(npts_avg, datenum(aquadoppL1.dtime), npts_avg, @rectwin);
    aquadoppL1.averaged.dtime = datetime(aquadoppL1.averaged.dtime, 'ConvertFrom', 'datenum');
    aquadoppL1.averaged.dtime.TimeZone = aquadoppL1.dtime.TimeZone;
    %
    aquadoppL1.averaged.dtime = aquadoppL1.averaged.dtime(:);

    %
    list_uvw = {'u', 'v', 'w'};

    % Now average the beam data
    prealloc_aux = NaN(size(aquadoppL1.Ue, 1), length(aquadoppL1.averaged.dtime));
    %
    for i2 = 1:length(list_uvw)
        aquadoppL1.averaged.(list_uvw{i2}) = prealloc_aux;
    end
    %
    aquadoppL1.averaged.a1 = prealloc_aux;
    aquadoppL1.averaged.a2 = prealloc_aux;
    aquadoppL1.averaged.a3 = prealloc_aux;

    % Loop over rows
    for i2 = 1:length(aquadoppL1.zhab)
        %
        for i3 = 1:length(list_uvw)
            aquadoppL1.averaged.(list_uvw{i3})(i2, :) = regRunMean(npts_avg, aquadoppL1.(list_uvw{i3})(i2, :), npts_avg, @rectwin);
        end
        %
        aquadoppL1.averaged.a1(i2, :) = regRunMean(npts_avg, aquadoppL1.a1(i2, :), npts_avg, @rectwin);
        aquadoppL1.averaged.a2(i2, :) = regRunMean(npts_avg, aquadoppL1.a2(i2, :), npts_avg, @rectwin);
        aquadoppL1.averaged.a3(i2, :) = regRunMean(npts_avg, aquadoppL1.a3(i2, :), npts_avg, @rectwin);
    end

    % ----------------------------------------------------
    % Remove ENU velocities
    aquadoppL1 = rmfield(aquadoppL1, {'Ue', 'Vn', 'Wup'});

    % ----------------------------------------------------
    % Add README
    time_dataproc = datetime('now', 'TimeZone', 'Local');
    time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');
    % Add a general README
    aquadoppL1.README = ['Level 1 Aquadopp data from ROXSI 2022. The data is from Aquadopp ' ...
                         'with serial number SN and deployed at mooring site mooringID. ' ...
                         'Data processed by script ' mfilename() '.m on ' time_dataproc_char ' (TimeZone ' time_dataproc.TimeZone '). ' ...
                         'Velocity components (u, v, w) are relative to the local cartesian ' ...
                         'coordinate system of the ROXSI site (i.e. where positive u is onshore ' ...
                         'and positive v is 90 degrees counterclockwise). w is positive for ' ...
                         'upward velocity. Magnetic declination was taken from www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml . ' ...
                         'Data in the Level 1 structure has been trimmed for the deployment ' ...
                         'period, as defined in the table deploymentInfo_ROXSI2022.mat. ' ...
                         'Pressure is in dbar, where atmospheric pressure has been removed.'];


    %% Make QC plots, save figues and data structure

    %
    disp('--- Done with processing. Making QC figures ---')

    % ----------------------------------------------------
    % Make a diagnostic plot of pressure, heading, pitch and roll
    fig_L1_QC_tilt = Aquadopp_scalars_QCplot(aquadoppL1, time_1, time_2);

    % ----------------------------------------------------
    % Make a diagnostic plot with pcolor
    %
    fig_L1_QC = Aquadopp_pcolor_lvl_1(aquadoppL1);

    
    % ----------------------------------------------------
    % Save level 1 data and QC figures

    %
    if lsave_file
        %
        disp('----- Saving level 1 data -----')

        %
        str_filename = ['roxsi_aquadopp_L1_' list_Aquadopp{i}(1:3) '_' char(aquadoppL1.SN)];
        %
        save(fullfile(dir_output_data_L1, [str_filename '.mat']), 'aquadoppL1', '-v7.3')
    end

    %
    if lsave_fig
        %
        disp('----- Saving level 1 QC plot figures -----')

        %
        str_filename = ['roxsi_aquadopp_L1_' list_Aquadopp{i}(1:3) '_' char(aquadoppL1.SN) '_QC_1'];
        % Save figure as *.png
        exportgraphics(fig_L1_QC_tilt, fullfile(dir_output_figs_L1, [str_filename '.png']), 'Resolution', 300)

        %
        str_filename = ['roxsi_aquadopp_L1_' list_Aquadopp{i}(1:3) '_' char(aquadoppL1.SN) '_QC_2'];
        % Save figure as *.png
        exportgraphics(fig_L1_QC, fullfile(dir_output_figs_L1, [str_filename '.png']), 'Resolution', 300)
    
        % Save figure as *.fig -- may not work
        % % savefig(fig_L1_QC, fullfile(dir_output_Aquadopp, [str_filename '.fig']), 'compact')
    end

    
    % Clear some variables to avoid issues in the next loop iteration
    close(fig_L1_QC_tilt), close(fig_L1_QC)
    clear aquadoppL1 header_aux beamAQDP_aux senAQDP_aux

    
    % ----------------------------------------------------
    %
    disp(['----- DONE with raw Aquadopp data proc: ' list_Aquadopp{i} ' -----'])
    %
    toc
end



%%

%
disp('###################### Done with data processing for all Aquadopps ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');

