%% Bathymetry around Spotters
%
% (ideally for Smart Moorings as well, although
% we don't have high-resolution bathymetry around
% Smart Moorings yet).

clear
close all


%% Load (CSUMB) bathymetry

%
dir_CSUMB = '/Users/olavobm/Library/CloudStorage/Box-Box/olavo_jamie/ROXSI_experiments/';

%
bathyCSUMB.CyPt_PtPn = load(fullfile(dir_CSUMB, 'bathy_CyPt_PtPn2m_xyz.mat'));

% Trim for both China Rock Asilomar
northing_box_lims = [4051000, 4054500];

%
lsub_aux = (bathyCSUMB.CyPt_PtPn.utm_north >= northing_box_lims(1)) & ...
           (bathyCSUMB.CyPt_PtPn.utm_north <= northing_box_lims(2));

%
bathyCSUMB.CyPt_PtPn.utm_east = bathyCSUMB.CyPt_PtPn.utm_east(lsub_aux);
bathyCSUMB.CyPt_PtPn.utm_north = bathyCSUMB.CyPt_PtPn.utm_north(lsub_aux);
bathyCSUMB.CyPt_PtPn.z_msl = bathyCSUMB.CyPt_PtPn.z_msl(lsub_aux);


%% Mooring table with reference locations

%
mooringLocs = load(['/Users/olavobm/Documents/ROXSI_Postdoc' ...
                    '/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022' ...
                    '/code_proc/ROXSI2022_mooringtable.mat']);
mooringLocs = mooringLocs.mooringtable;


%
[mooringLocs.easting, ...
 mooringLocs.northing] = ll2utm(mooringLocs.latitude, mooringLocs.longitude);


%% Load Spotter deployment info table

dplySpotters = load(['/Users/olavobm/Documents/ROXSI_Postdoc' ...
                     '/MyResearch/ROXSI/Common_Code' ...
                     '/LargeScale_Data_2022/code_proc' ...
                     '/deploymentInfo_Spotters_ROXSI2022.mat']);
dplySpotters = dplySpotters.dployInfo_Spotters;


%% Spotters that will be loaded

% % %
% % dir_data = ['/Volumes/GoogleDrive/Shared drives/ROXSI' ...
% %             '/LargeScale_Data_2022/Level1_Data/Spotter_Level1/'];

%
dir_data = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';

%
list_files = {'B01_spot1150.mat', 'B01_spot1158.mat', 'B03_spot1152.mat', ...
              'B05_spot1153.mat', 'X01_spot1151.mat', 'X03_spot1157.mat', ...
              'X04_spot1155.mat'};


%% Load tides (I could also compute a tide product from the
% pressure sensors instead of using Monterey harbor tides.
% I should probably do that, but the difference will probably
% be minor)

%
tidal_elevation = load(['/Volumes/ROXSI_Data/LargeScale_Data_2022', ...
                        '/RAW/noaa_mry_tides/tides_NOAA_MRY.mat']);
tidal_elevation = tidal_elevation.noaaTides;


%% Load data and get trimming variable

%
for i = 1:length(list_files)
    
    %
    data_aux = load(fullfile(dir_data, list_files{i}));
    data_aux = data_aux.s;
    
    %
    data_aux.location.time.TimeZone = 'America/Los_Angeles';
    
    %
    spotterAll(i).dataID = list_files{i}(1:12);
    %
    spotterAll(i).location = data_aux.location;


    % Trim Spotter location data
    %
    lmatch = strcmp(dplySpotters.SN, list_files{i}(9:12));
    
    %
    time_trim_1 = dplySpotters.time_begin_trim(lmatch);
    time_trim_2 = dplySpotters.time_end_trim(lmatch);
    %
    time_trim_1 = datetime(datenum(time_trim_1, "yyyy/mm/dd HH:MM:SS"), 'ConvertFrom', 'datenum');
    time_trim_2 = datetime(datenum(time_trim_2, "yyyy/mm/dd HH:MM:SS"), 'ConvertFrom', 'datenum');
    %
    time_trim_1.TimeZone = 'America/Los_Angeles';
    time_trim_2.TimeZone = 'America/Los_Angeles';
    
    %
    spotterAll(i).timetrimedges = [time_trim_1, time_trim_2];

    %
    spotterAll(i).ltrimedges = (spotterAll(i).location.time >= time_trim_1) & ...
                                  (spotterAll(i).location.time <= time_trim_2);

    % Compute easting/northing coordinates
    [spotterAll(i).location.easting, ...
     spotterAll(i).location.northing] = ll2utm(spotterAll(i).location.("latitude (decimal degrees)"), ...
                                               spotterAll(i).location.("longitude (decimal degrees)"));

    %
    clear data_aux
end


%% Get bathymetry around each Spotter


%
for i = 1:length(list_files)

    %
    mean_easting = mean(spotterAll(i).location.easting(spotterAll(i).ltrimedges));
    mean_northing = mean(spotterAll(i).location.northing(spotterAll(i).ltrimedges));

    %
    minmax_easting = [min(spotterAll(i).location.easting(spotterAll(i).ltrimedges)), max(spotterAll(i).location.easting(spotterAll(i).ltrimedges))];
    minmax_northing = [min(spotterAll(i).location.northing(spotterAll(i).ltrimedges)), max(spotterAll(i).location.northing(spotterAll(i).ltrimedges))];
    
    %
    max_easting_dist = max(abs(minmax_easting - mean_easting));
    max_northing_dist = max(abs(minmax_northing - mean_northing));

    %
    max_dist = max([max_easting_dist, max_northing_dist]);

    %
    dist_get_bathy = 1.35*max_dist;

    %
    lbathyinlims = (bathyCSUMB.CyPt_PtPn.utm_east >= (mean_easting - dist_get_bathy)) & ...
                   (bathyCSUMB.CyPt_PtPn.utm_east <= (mean_easting + dist_get_bathy)) & ...
                   (bathyCSUMB.CyPt_PtPn.utm_north >= (mean_northing - dist_get_bathy)) & ...
                   (bathyCSUMB.CyPt_PtPn.utm_north <= (mean_northing + dist_get_bathy));
    %
    spotterAll(i).bathy.easting = bathyCSUMB.CyPt_PtPn.utm_east(lbathyinlims);
    spotterAll(i).bathy.northing = bathyCSUMB.CyPt_PtPn.utm_north(lbathyinlims);
    spotterAll(i).bathy.z_msl = bathyCSUMB.CyPt_PtPn.z_msl(lbathyinlims);
    %
    spotterAll(i).bathy.bathyInterpolant = ...
                            scatteredInterpolant(spotterAll(i).bathy.easting, ...
                                                 spotterAll(i).bathy.northing, ...
                                                 spotterAll(i).bathy.z_msl);

end


%% Interpolate bathymetry in a box around Spotters,
% define smoothing scale, interpolate smoothed
% bathymetry around Spotter, and interpolate bathymetry
% at every location of the spotter

% In meters
dx_grid = 1;
% % dx_smoothing = 25;
dx_smoothing = 11;   
%
disk_filter_weights = fspecial("disk", dx_smoothing);  % THIS ASSUMES GRID
                                                       % SPACING IS 1 M,
                                                       % BECAUSE THE INPUT
                                                       % IS IN TERMS OF
                                                       % NUMBER OF POINTS
                                                       % (NOT DISTANCE)

tic
%
for i = 1:length(spotterAll)

    %
    spotterAll(i).bathy.gridded.easting = ceil(min(spotterAll(i).bathy.easting)) : dx_grid : floor(max(spotterAll(i).bathy.easting));
    spotterAll(i).bathy.gridded.northing = ceil(min(spotterAll(i).bathy.northing)) : dx_grid : floor(max(spotterAll(i).bathy.northing));

    %
    [x_mesh, y_mesh] = meshgrid(spotterAll(i).bathy.gridded.easting, ...
                                spotterAll(i).bathy.gridded.northing);
    %
    spotterAll(i).bathy.gridded.z_msl = ...
                    spotterAll(i).bathy.bathyInterpolant(x_mesh, ...
                                                         y_mesh);
    %
    spotterAll(i).bathy.gridded.smoothing_radius = dx_smoothing;

    %
    spotterAll(i).bathy.gridded.z_msl_smoothed = ...
                        imfilter(spotterAll(i).bathy.gridded.z_msl, ...
                                 disk_filter_weights);
    
    % NaN edges of smoothed product
    spotterAll(i).bathy.gridded.z_msl_smoothed(1:dx_smoothing, :) = NaN;
    spotterAll(i).bathy.gridded.z_msl_smoothed(end-dx_smoothing:end, :) = NaN;
    spotterAll(i).bathy.gridded.z_msl_smoothed(:, 1:dx_smoothing) = NaN;
    spotterAll(i).bathy.gridded.z_msl_smoothed(:, end-dx_smoothing:end) = NaN;


    %
    spotterAll(i).bathy.smoothedbathyInterpolant = scatteredInterpolant(x_mesh(:), ...
                                                                        y_mesh(:), ...
                                                                        spotterAll(i).bathy.gridded.z_msl_smoothed(:));

    %
    spotterAll(i).location.z_msl = ...
               spotterAll(i).bathy.smoothedbathyInterpolant(spotterAll(i).location.easting, ...
                                                            spotterAll(i).location.northing);

end
toc


%% Now take into account tidal elevation to compute the bottom depth
%
% Note the sign definition of each variable

%
for i = 1:length(spotterAll)

    %
    interp_tidal_elevation = interp1(tidal_elevation.dtime, tidal_elevation.zMSL, ...
                                     spotterAll(i).location.time);

    % Minus tidal elevation means that z_msl gets bigger
    % bigger magnitude (deeper water depth)
    spotterAll(i).location.z_msl = spotterAll(i).location.z_msl - interp_tidal_elevation;

end



%%
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ---------------------------- MAKE FIGURES ----------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%% Plot coverage

% % %
% % figure
% %     %
% %     scatter(bathyCSUMB.CyPt_PtPn.utm_east(1:100:end), ...
% %             bathyCSUMB.CyPt_PtPn.utm_north(1:100:end), ...
% %             60, bathyCSUMB.CyPt_PtPn.z_msl(1:100:end), 'filled')
% %     %
% %     hold on
% % 
% %     %
% %     lspotters = (mooringLocs.instrument == "spotter") | ...
% %                 (mooringLocs.instrument == "smartspotter");
% %     plot(mooringLocs.easting(lspotters), ...
% %          mooringLocs.northing(lspotters), '.k', 'MarkerSize', 32)
% % 
% % % %     %
% % % %     plot(easting_aux(lintrim), northing_aux(lintrim), '.r')
% % 
% %     %
% %     axis equal


%% Plot bathy around each Spotter

% % %
% % for i = 1:length(list_files)
% % 
% %     %
% %     figure
% %         %
% %         scatter(spotterAll(i).bathy.easting, ...
% %                 spotterAll(i).bathy.northing, ...
% %                 120, spotterAll(i).bathy.z_msl, 'filled')
% %         %
% %         hold on
% % 
% %         % Plot Spotter location (a subset so the plot is less crowded) 
% %         indsplt_aux = find(spotterAll(i).ltrimedges);
% %         %
% %         plot(spotterAll(i).location.easting(indsplt_aux(1:20:end)), ...
% %              spotterAll(i).location.northing(indsplt_aux(1:20:end)), '.k')
% %     
% %         % Plot reference location
% %         lmatch_mooringtable = strcmp(mooringLocs.mooringID, [spotterAll(i).dataID(1:3) 's']);
% %         plot(mooringLocs.easting(lmatch_mooringtable), ...
% %              mooringLocs.northing(lmatch_mooringtable), '.r', 'MarkerSize', 32)
% %     
% %         %
% %         hcb = colorbar;
% %             hcb.Label.String = 'z mean sea level [m]';
% %             hcb.Label.Interpreter = 'Latex';
% %             hcb.Label.FontSize = 28;
% % 
% %         %
% %         set(gca, 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% %         set(gca, 'DataAspectRatio', [1, 1, 1])
% % 
% %         %
% %         xlabel(['Easting (' num2str(round(diff(xlim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
% %         ylabel(['Northing (' num2str(round(diff(ylim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
% %         %
% %         title([spotterAll(i).dataID(1:3) ' SN ' spotterAll(i).dataID(9:end)], ...
% %               'Interpreter', 'Latex', 'FontSize', 20)
% % 
% % end

%% Same as above, but with gridded unsmoothed bathymetry

% % %
% % for i = 1:length(list_files)
% % 
% %     %
% %     figure
% %         %
% %         pcolor(spotterAll(i).bathy.gridded.easting, ...
% %                spotterAll(i).bathy.gridded.northing, ...
% %                spotterAll(i).bathy.gridded.z_msl)
% %         shading flat
% %         %
% %         hold on
% % 
% %         % Plot Spotter location (a subset so the plot is less crowded) 
% %         indsplt_aux = find(spotterAll(i).ltrimedges);
% %         %
% %         plot(spotterAll(i).location.easting(indsplt_aux(1:20:end)), ...
% %              spotterAll(i).location.northing(indsplt_aux(1:20:end)), '.k')
% %     
% %         % Plot reference location
% %         lmatch_mooringtable = strcmp(mooringLocs.mooringID, [spotterAll(i).dataID(1:3) 's']);
% %         plot(mooringLocs.easting(lmatch_mooringtable), ...
% %              mooringLocs.northing(lmatch_mooringtable), '.r', 'MarkerSize', 32)
% %     
% %         %
% %         hcb = colorbar;
% %             hcb.Label.String = 'z mean sea level [m]';
% %             hcb.Label.Interpreter = 'Latex';
% %             hcb.Label.FontSize = 28;
% % 
% %         %
% %         caxis([min(spotterAll(i).location.z_msl), ...
% %                max(spotterAll(i).location.z_msl)])
% % 
% %         %
% %         set(gca, 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% %         set(gca, 'DataAspectRatio', [1, 1, 1])
% % 
% %         %
% %         xlabel(['Easting (' num2str(round(diff(xlim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
% %         ylabel(['Northing (' num2str(round(diff(ylim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
% %         %
% %         title([spotterAll(i).dataID(1:3) ' SN ' spotterAll(i).dataID(9:end)], ...
% %               'Interpreter', 'Latex', 'FontSize', 20)
% % 
% % end

%% Same as above, but with gridded smoothed bathymetry

%
for i = 1:length(list_files)

    %
    hfig_aux = figure;
        %
        pcolor(spotterAll(i).bathy.gridded.easting, ...
               spotterAll(i).bathy.gridded.northing, ...
               spotterAll(i).bathy.gridded.z_msl_smoothed)
        shading flat
        %
        hold on

        % Plot Spotter location (a subset so the plot is less crowded) 
        indsplt_aux = find(spotterAll(i).ltrimedges);
        %
        plot(spotterAll(i).location.easting(indsplt_aux(1:20:end)), ...
             spotterAll(i).location.northing(indsplt_aux(1:20:end)), '.k')
    
        % Plot reference location
        lmatch_mooringtable = strcmp(mooringLocs.mooringID, [spotterAll(i).dataID(1:3) 's']);
        plot(mooringLocs.easting(lmatch_mooringtable), ...
             mooringLocs.northing(lmatch_mooringtable), '.r', 'MarkerSize', 32)
    
        %
        hcb = colorbar;
            hcb.Label.String = 'z mean sea level [m]';
            hcb.Label.Interpreter = 'Latex';
            hcb.Label.FontSize = 28;

        %
        caxis([min(spotterAll(i).location.z_msl), ...
               max(spotterAll(i).location.z_msl)])

        %
        set(gca, 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
        set(gca, 'DataAspectRatio', [1, 1, 1])

        %
        xlabel(['Easting (' num2str(round(diff(xlim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
        ylabel(['Northing (' num2str(round(diff(ylim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
        %
        title([spotterAll(i).dataID(1:3) ' SN ' spotterAll(i).dataID(9:end)], ...
              'Interpreter', 'Latex', 'FontSize', 20)

    %
    exportgraphics(hfig_aux, fullfile(dir_data, [list_files{i}(1:end-4) '.png']), 'Resolution', 300)

end


%% Plot timeseries of bathymetry for each Spotter
% (NOT ADDING TIDAL ELEVATION THOUGH!!!)

%
for i = 1:length(spotterAll)

    %
    figure
        %
        plot(spotterAll(i).location.time, ...
             spotterAll(i).location.z_msl, '.-')

        %
        set(gca, 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')

% %         %
% %         xlabel(['Easting (' num2str(round(diff(xlim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
% %         ylabel(['Northing (' num2str(round(diff(ylim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
        %
        title([spotterAll(i).dataID(1:3) ' SN ' spotterAll(i).dataID(9:end)], ...
              'Interpreter', 'Latex', 'FontSize', 20)
end


%%
% -----------------------------------------------
% ----- GET SMART MOORING DATA AND ONLY USE -----
% -------- TIDAL ELEVATION AND PRESSURE ---------
% ---------- BECAUSE WE DON'T HAVE YET ----------
% ---- GREAT BATHYMETRY AROUND THESE MOORINGS ---
% -----------------------------------------------
%
% NOTE: Median pressure is slightly (<= 10 cm) higher
% than mean pressure (I checked the median to make
% sure the gaps were not a big problem for defining
% the "true mean" depth of the Spotters.


%% Since I don't have the bathymetry around Smart moorings,
% I will just use the pressure data for now

%
dir_SmartMoorings_pressure = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Smart_Level1/gridded/';
dir_SmartMoorings_buoy = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';

%
list_SmartMoorings = {'E01_spot1851', 'E02_spot1859', ...
                      'E05_spot1853', 'E07_spot1855', 'E07_spot1857', ...
                      'E08_spot1852', 'E09_spot1850', 'E09_spot1856', ...
                      'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};

% % % Significant wave height and trimming for Smart Moorings
% % %
% % for i = 1:length(list_SmartMoorings)
% % % %     %
% % % %     data_aux = load(fullfile(dir_SmartMoorings_pressure, ...
% % % %                     ['smart_mooring_' list_SmartMoorings{i}(1:3) 'sp_' list_SmartMoorings{i}(9:12) '_L1_gridded.mat']));
% % 
% %     %
% %     data_aux = load(fullfile(dir_SmartMoorings_buoy, ...
% %                     [list_SmartMoorings{i} '.mat']));
% %     data_aux = data_aux.s;
% %     %
% %     data_aux.bulkparameters.time.TimeZone = 'America/Los_Angeles';
% % 
% %     %
% %     lmatch_aux = strcmp(dplySpotters.SN, list_SmartMoorings{i}(9:12));
% % 
% %     %
% %     figure
% %         %
% %         plot(data_aux.bulkparameters.time, ...
% %              data_aux.bulkparameters.("Significant Wave Height"), '.-')
% %         %
% %         hold on
% % 
% %             
% %         %
% %         time_trim_1 = dplySpotters.time_begin_trim(lmatch_aux);
% %         time_trim_2 = dplySpotters.time_end_trim(lmatch_aux);
% %         %
% %         time_trim_1 = datetime(datenum(time_trim_1, "yyyy/mm/dd HH:MM:SS"), 'ConvertFrom', 'datenum');
% %         time_trim_2 = datetime(datenum(time_trim_2, "yyyy/mm/dd HH:MM:SS"), 'ConvertFrom', 'datenum');
% %         %
% %         time_trim_1.TimeZone = 'America/Los_Angeles';
% %         time_trim_2.TimeZone = 'America/Los_Angeles';
% % 
% %         %
% %         overlayline('v', [time_trim_1, time_trim_2], '--k')
% % 
% %         %
% %         set(gca, 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% %         %
% %         title([list_SmartMoorings{i}(1:3) ' SN ' list_SmartMoorings{i}(9:12)], ...
% %               'Interpreter', 'Latex', 'FontSize', 20)
% %         %
% %         set(gcf, 'units', 'normalized')
% %         set(gcf, 'Position', [0.3082, 0.4785, 0.4492, 0.2944])
% % 
% % end

%%

%
for i = 1:length(list_SmartMoorings)


    % Buoy data
    %
    data_BUOY = load(fullfile(dir_SmartMoorings_buoy, ...
                    [list_SmartMoorings{i} '.mat']));
    data_BUOY = data_BUOY.s;
    %
    data_BUOY.bulkparameters.time.TimeZone = 'America/Los_Angeles';
    data_BUOY.location.time.TimeZone = 'America/Los_Angeles';


    % Pressure data
    data_PRES = load(fullfile(dir_SmartMoorings_pressure, ...
                     ['smart_mooring_' list_SmartMoorings{i}(1:3) ...
                      'sp_' list_SmartMoorings{i}(9:12) '_L1_gridded.mat']));
    data_PRES = data_PRES.spotsmart;

    
% %     % Check mean pressure and median pressure and if gaps are a problem
% %     figure
% %         plot(data_PRES.dtime, data_PRES.pressure)
% %         %
% %         grid on
% %         %
% %         overlayline('h', mean(data_PRES.pressure, 'omitnan'), '-k')
% %         overlayline('h', median(data_PRES.pressure, 'omitnan'), '-r')

   
    % Compute easting/northing coordinates
    [data_BUOY.location.easting, ...
     data_BUOY.location.northing] = ll2utm(data_BUOY.location.("latitude (decimal degrees)"), ...
                                           data_BUOY.location.("longitude (decimal degrees)"));

    %
    lmatch_aux = strcmp(dplySpotters.SN, list_SmartMoorings{i}(9:12));
    %
    time_trim_1 = dplySpotters.time_begin_trim(lmatch_aux);
    time_trim_2 = dplySpotters.time_end_trim(lmatch_aux);
    %
    time_trim_1 = datetime(datenum(time_trim_1, "yyyy/mm/dd HH:MM:SS"), 'ConvertFrom', 'datenum');
    time_trim_2 = datetime(datenum(time_trim_2, "yyyy/mm/dd HH:MM:SS"), 'ConvertFrom', 'datenum');
    %
    time_trim_1.TimeZone = 'America/Los_Angeles';
    time_trim_2.TimeZone = 'America/Los_Angeles';
    %
    if strcmp(list_SmartMoorings{i}, 'E09_spot1850')
        time_trim_2 = datetime(2022, 06, 24, 06, 00, 00, 'TimeZone', 'America/Los_Angeles');
    end
    %
    ltrimedges = (data_BUOY.location.time >= time_trim_1) & ...
                 (data_BUOY.location.time <= time_trim_2);

    %
    data_BUOY.location.z_msl = NaN(size(data_BUOY.location, 1), 1);
    
    %
    tides_interp_aux = interp1(tidal_elevation.dtime, ...
                               tidal_elevation.zMSL, ...
                               data_BUOY.location.time(ltrimedges));

    % Save depth to Smart moorings (without taking the bathymetry into
    % account)
    data_BUOY.location.z_msl(ltrimedges) = -mean(data_PRES.pressure, 'omitnan') - tides_interp_aux;
    % RIGOROUSLY, SHOULD INCLUDE THE CONVERSION FROM DBAR TO M
    
    %
    spotterSmartAll(i).dataID = list_SmartMoorings{i};
    spotterSmartAll(i).location = data_BUOY.location;
    spotterSmartAll(i).timetrimedges = [time_trim_1, time_trim_2];
    spotterSmartAll(i).ltrimedges = ltrimedges;
    
end




%%
% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------ SAVE OUTPUT THAT IS REQUIRED FOR DIRECTIONAL SPECTRA ------
% ------------------------------------------------------------------
% ------------------------------------------------------------------

%% Now save output for all spotters in a good format

%
dir_output = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';

%
spotter_location = struct('dataID', [], 'location', [], ...
                          'timetrimedges', [], 'ltrimedges', [], ...
                          'bathy', []);

%
for i = 1:length(spotterAll)
    %
    spotter_location(i).dataID = spotterAll(i).dataID;
    spotter_location(i).location = spotterAll(i).location;
    spotter_location(i).timetrimedges = spotterAll(i).timetrimedges;
    spotter_location(i).ltrimedges = spotterAll(i).ltrimedges;
    spotter_location(i).bathy = spotterAll(i).bathy;

end


%
for i = 1:length(spotterSmartAll)
    %
    ind_aux = length(spotterAll) + i;
    %
    spotter_location(ind_aux).dataID = spotterSmartAll(i).dataID;
    spotter_location(ind_aux).location = spotterSmartAll(i).location;
    spotter_location(ind_aux).timetrimedges = spotterSmartAll(i).timetrimedges;
    spotter_location(ind_aux).ltrimedges = spotterSmartAll(i).ltrimedges;

end


%
disp(' '), disp(' ')
disp(['Saving Spotter location and depth at ' dir_output])
save('-v7.3', fullfile(dir_output, 'Spotter_all_location_depth.mat'), 'spotter_location')


