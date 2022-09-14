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
dir_data = '/Volumes/LaCie/RAW/Spotters/SDcards/';

%
list_files = {'B01_spot1150.mat', 'B01_spot1158.mat', 'B03_spot1152.mat', ...
              'B05_spot1153.mat', 'X01_spot1151.mat', 'X03_spot1157.mat', ...
              'X04_spot1155.mat'};


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
    spotterAll(i).llocationtrim = (spotterAll(i).location.time >= time_trim_1) & ...
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
    mean_easting = mean(spotterAll(i).location.easting(spotterAll(i).llocationtrim));
    mean_northing = mean(spotterAll(i).location.northing(spotterAll(i).llocationtrim));

    %
    minmax_easting = [min(spotterAll(i).location.easting(spotterAll(i).llocationtrim)), max(spotterAll(i).location.easting(spotterAll(i).llocationtrim))];
    minmax_northing = [min(spotterAll(i).location.northing(spotterAll(i).llocationtrim)), max(spotterAll(i).location.northing(spotterAll(i).llocationtrim))];
    
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

%%
% -------------------------------------------------------
% -------------------------------------------------------
% -------------------------------------------------------

%% Plot coverage

%
figure
    %
    scatter(bathyCSUMB.CyPt_PtPn.utm_east(1:100:end), ...
            bathyCSUMB.CyPt_PtPn.utm_north(1:100:end), ...
            60, bathyCSUMB.CyPt_PtPn.z_msl(1:100:end), 'filled')
    %
    hold on

    %
    lspotters = (mooringLocs.instrument == "spotter") | ...
                (mooringLocs.instrument == "smartspotter");
    plot(mooringLocs.easting(lspotters), ...
         mooringLocs.northing(lspotters), '.k', 'MarkerSize', 32)

% %     %
% %     plot(easting_aux(lintrim), northing_aux(lintrim), '.r')

    %
    axis equal


%% Plot bathy around each Spotter

%
for i = 1:length(list_files)

    %
    figure
        %
        scatter(spotterAll(i).bathy.easting, ...
                spotterAll(i).bathy.northing, ...
                120, spotterAll(i).bathy.z_msl, 'filled')
        %
        hold on

        % Plot Spotter location (a subset so the plot is less crowded) 
        indsplt_aux = find(spotterAll(i).llocationtrim);
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
        set(gca, 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
        set(gca, 'DataAspectRatio', [1, 1, 1])

        %
        xlabel(['Easting (' num2str(round(diff(xlim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
        ylabel(['Northing (' num2str(round(diff(ylim))) ' m)'], 'Interpreter', 'Latex', 'FontSize', 20)
        %
        title([spotterAll(i).dataID(1:3) ' SN ' spotterAll(i).dataID(9:end)], ...
              'Interpreter', 'Latex', 'FontSize', 20)

end

%% Same as above, but with gridded unsmoothed bathymetry

%
for i = 1:length(list_files)

    %
    figure
        %
        pcolor(spotterAll(i).bathy.gridded.easting, ...
               spotterAll(i).bathy.gridded.northing, ...
               spotterAll(i).bathy.gridded.z_msl)
        shading flat
        %
        hold on

        % Plot Spotter location (a subset so the plot is less crowded) 
        indsplt_aux = find(spotterAll(i).llocationtrim);
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

end

%% Same as above, but with gridded smoothed bathymetry

%
for i = 1:length(list_files)

    %
    figure
        %
        pcolor(spotterAll(i).bathy.gridded.easting, ...
               spotterAll(i).bathy.gridded.northing, ...
               spotterAll(i).bathy.gridded.z_msl_smoothed)
        shading flat
        %
        hold on

        % Plot Spotter location (a subset so the plot is less crowded) 
        indsplt_aux = find(spotterAll(i).llocationtrim);
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


