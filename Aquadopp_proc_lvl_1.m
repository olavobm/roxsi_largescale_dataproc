clear
close all



%% Get 

% This creates the table ROXSI2022_deploymentInfo
run(['/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI' ...
     '/other_matlab_functions/codebyOBM/obm_ROXSI_toolbox' ...
     '/instruments_toolbox/define_ROXSI2022_deploymentinfo.m'])


%%

%
dir_data_parent = '/Volumes/LaCie/RAW/Aquadopp/';

% % %
% % list_Aquadopp = {'A03_5380', ...
% %                  'B02_12507', 'B04_2147', 'B07_2141', 'B08_13288', 'B11_12280', ...
% %                  'C03_0709', ...
% %                  'D01_12346', 'D02_0653', ...
% %                  'E03_13300', 'E04_13172', 'E06_9736', 'E12_11150', ...
% %                  'F01_9995', 'F02_5838', 'F03_5384', 'F04_5401', 'F05_14032', ...
% %                  'X06_13290', 'X13_9945'};

%
list_Aquadopp = {'B02_12507'};
% list_Aquadopp = {'B07_2141'};

%
Naquadopps = length(list_Aquadopp);

%
dir_output_Aquadopp = ['/Volumes/OBM-HD/docs/researchPostdoc/datasets' ...
                       '/ROXSI/fieldworks/experiment_2022/ADCPs/status/Aquadopp/'];

%% Load atmospheric pressure

%
atm_pressure = load(['/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch' ...
                     '/figures_bydate/2022_08_17/obm_edited_noaa_mry_barometric_pressure/' ...
                     'atm_pressure.mat']);




%%
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------


%%

%
list_senfile_vars = ["time", "heading", "pitch", "roll", "pressure", "temperature"];

%
for i = 1:Naquadopps

    % ----------------------------------------------------
    % Get header (in particular, cell centers)
    file_header_aux = dir(fullfile(dir_data_parent, list_Aquadopp{i}, '*.hdr'));
    %
    header_aux = Aquadopp_read_header(fullfile(file_header_aux.folder, file_header_aux.name));
    
    % ----------------------------------------------------
    %
    file_sen_aux = dir(fullfile(dir_data_parent, list_Aquadopp{i}, '*.sen'));

    %
    senAQDP_aux = Aquadopp_read_senfile(fullfile(file_sen_aux.folder, file_sen_aux.name), list_senfile_vars);
    

    
% % %     %
% % %     senAQDP_aux.time = datetime(datenum(senAQDP_aux.time, 'yyyy/mm/dd HH:MM:SS'), ...
% % %                                         'ConvertFrom', 'datenum', ...
% % %                                         'TimeZone', 'America/Los_Angeles');


    % Now correct for clock drift
    time_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, list_Aquadopp{i}(5:end), datenum(senAQDP_aux.time));
    %
    senAQDP_aux.time = datetime(time_aux, ...
                                    'ConvertFrom', 'datenum', ...
                                    'TimeZone', 'America/Los_Angeles');


    % ----------------------------------------------------
    % Load velocity and backscatter data

    %
    file_noextension_aux = file_sen_aux.name(1:(strfind(file_sen_aux.name, '.') - 1));

    %
    beamAQDP_aux = Aquadopp_read_beamdata(fullfile(dir_data_parent, list_Aquadopp{i}, file_noextension_aux), ...
                                          ["a1", "a2", "a3", "v1", "v2", "v3"]);

   
    

    % ----------------------------------------------------
    % Trim for when instrument was in the water
    ind_row_match = find(strcmp(deploymentInfo_ROXSI2022.SN, list_Aquadopp{i}(5:end)));
    
    %
    time_1 = deploymentInfo_ROXSI2022.time_begin_trim(ind_row_match);
    time_2 = deploymentInfo_ROXSI2022.time_end_trim(ind_row_match);

    %
    time_1 = datetime(datenum(time_1, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    time_2 = datetime(datenum(time_2, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');

    %
    lin_water = (senAQDP_aux.time >= time_1) & (senAQDP_aux.time <= time_2);

    %
    statusAquadopp.all = senAQDP_aux;
    %
    for i2 = 1:length(list_senfile_vars)
        %
        statusAquadopp.trimmed.(list_senfile_vars(i2)) = statusAquadopp.all.(list_senfile_vars(i2))(lin_water);
    end

    % ----------------------------------------------------
    % Start adding data into level 1 data structure
    
    %
    lmatch_ontable = strcmp(deploymentInfo_ROXSI2022.SN, list_Aquadopp{i}(5:end));

    %
    aquadoppL1.SN = convertCharsToStrings(list_Aquadopp{i}(5:end));
    %
    aquadoppL1.mooringID = deploymentInfo_ROXSI2022.mooringID(lmatch_ontable);

    % Latitude/longitude

    % x/y

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
        aquadoppL1.transducerHAB = 0.5;   % PLACEHOLDER VALUE (0.5 m seems like a reasonable guess)
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
    aquadoppL1.timezone = 'PDT';
    %
    aquadoppL1.dtime = senAQDP_aux.time(lin_water).';
    aquadoppL1.heading = senAQDP_aux.heading(lin_water).';
    aquadoppL1.pitch = senAQDP_aux.pitch(lin_water).';
    aquadoppL1.roll = senAQDP_aux.roll(lin_water).';

    %
    aquadoppL1.pressureoffset = deploymentInfo_ROXSI2022.pressure_offset(ind_row_match);
    
    %
    aquadoppL1.pressure = senAQDP_aux.pressure(lin_water).' - aquadoppL1.pressureoffset;

    % ------------------------------
    % Correct for the atmospheric pressure variability
    %
    keyboard


    % ------------------------------

    %
    aquadoppL1.temperature = senAQDP_aux.temperature(lin_water).';

    % ----------------------------------------------------
    %
    aquadoppL1.coordsystem = 'ENU';    % sort of because that's not true yet

    % In clockwise degrees from the true north
    aquadoppL1.magdec = 12.86;

    % At this point, Ue and Vn are relative to the MAGNETIC north
    aquadoppL1.Ue = beamAQDP_aux.v1(:, lin_water);
    aquadoppL1.Vn = beamAQDP_aux.v2(:, lin_water);
    aquadoppL1.Wup = beamAQDP_aux.v3(:, lin_water);
    %
    aquadoppL1.a1 = beamAQDP_aux.a1(:, lin_water);
    aquadoppL1.a2 = beamAQDP_aux.a2(:, lin_water);
    aquadoppL1.a3 = beamAQDP_aux.a3(:, lin_water);

    % ----------------------------------------------------
    % Rotate horizontal velocity from magnetic to true north

    %
    rotMatrix = [cosd(aquadoppL1.magdec), sind(aquadoppL1.magdec); ...
                 -sind(aquadoppL1.magdec), cosd(aquadoppL1.magdec)];

% %     % Check the rotation (i.e. velocity aligned with magnetic north
% %     % should have a small zonal component and large meridional
% %     % component)
% %     rotMatrix * [0; 1]

    %
    for i2 = 1:size(aquadoppL1.Ue, 1)
        %
        uv_aux = [aquadoppL1.Ue(i2, :); aquadoppL1.Vn(i2, :)];
        %
        uv_rot_aux = rotMatrix * uv_aux;

        %
        aquadoppL1.Ue(i2, :) = uv_aux(1, :);
        aquadoppL1.Vn(i2, :) = uv_aux(2, :);
    end
    
    % ----------------------------------------------------
    % Compute averaged fields (compute 10 and 20 min averages
    % are whole number choices because B02 gives data
    % every 10 min).

    %
    aquadoppL1.averaged.dt = 10*60;

    %
    npts_avg = aquadoppL1.averaged.dt / aquadoppL1.samplingtime;

    % First average time (also to get its length)
    aquadoppL1.averaged.dtime = regRunMean(npts_avg, datenum(aquadoppL1.dtime), npts_avg, @rectwin);
    aquadoppL1.averaged.dtime = datetime(aquadoppL1.averaged.dtime, 'ConvertFrom', 'datenum');
    aquadoppL1.averaged.dtime.TimeZone = aquadoppL1.dtime.TimeZone;

    % Now average the beam data
    prealloc_aux = NaN(size(aquadoppL1.Ue, 1), length(aquadoppL1.averaged.dtime));
    %
    aquadoppL1.averaged.Ue = prealloc_aux;
    aquadoppL1.averaged.Vn = prealloc_aux;
    aquadoppL1.averaged.Wup = prealloc_aux;
    %
    aquadoppL1.averaged.a1 = prealloc_aux;
    aquadoppL1.averaged.a2 = prealloc_aux;
    aquadoppL1.averaged.a3 = prealloc_aux;

    % Loop over rows
    for i2 = 1:length(aquadoppL1.cellcenter)
        %
        aquadoppL1.averaged.Ue(i2, :) = regRunMean(npts_avg, aquadoppL1.Ue(i2, :), npts_avg, @rectwin);
        aquadoppL1.averaged.Vn(i2, :) = regRunMean(npts_avg, aquadoppL1.Vn(i2, :), npts_avg, @rectwin);
        aquadoppL1.averaged.Wup(i2, :) = regRunMean(npts_avg, aquadoppL1.Wup(i2, :), npts_avg, @rectwin);
        %
        aquadoppL1.averaged.a1(i2, :) = regRunMean(npts_avg, aquadoppL1.a1(i2, :), npts_avg, @rectwin);
        aquadoppL1.averaged.a2(i2, :) = regRunMean(npts_avg, aquadoppL1.a2(i2, :), npts_avg, @rectwin);
        aquadoppL1.averaged.a3(i2, :) = regRunMean(npts_avg, aquadoppL1.a3(i2, :), npts_avg, @rectwin);
    end

    % ----------------------------------------------------
    % Make a diagnostic plot
    fig_L1_QC = Aquadopp_pcolor_lvl_1(aquadoppL1);


    % ----------------------------------------------------

% %     % Add a general README
% %     aquadoppL1.README = ['Aquadopp level 1 data structure. mooringID indicates the ' ...
% %                           mooring location according to the deployment plan in the ' ...
% %                           ROXSI 2022 experiment. Pressure is in dbar and ???: ' ...
% %                           hydrostatic or total? OFFSET? All time quantities are in seconds. ' ...
% %                           Magnetic declination taken from
% %                           www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml.'];


    % ----------------------------------------------------
    % Save level 1 data

% % %     %
% % %     str_filename = ['roxsi_aquadopp_L1_' char(aquadoppL1.SN) '_' list_Aquadopp{i}(1:3)];
% % % 
% % %     %
% % %     save(fullfile(dir_output_Aquadopp, [str_filename '.mat']), 'aquadopplvl1')
% % % 
% % %     % Save figure as *.png
% % %     exportgraphics(fig_L1_QC, fullfile(dir_output_Aquadopp, [str_filename '.png']), 'Resolution', 300)
% % % 
% % %     % Save figure as *.fig -- may not work
% % %     % % savefig(fig_L1_QC, fullfile(dir_output_Aquadopp, [str_filename '.fig']), 'compact')



end