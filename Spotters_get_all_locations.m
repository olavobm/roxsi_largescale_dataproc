%% Script to get reference location of all Spotters (including
% Smart moorings). Here, I "reference location" means what will be
% used as the location of pressure sensors on Smart Moorings. This
% is NOT the average location because the Spotters does not spend
% time homogenously along its watch circle. Instead a least-squares
% fit is used to calculate the "best circle" and the coordinate at
% its center is taken.
%
% Since I need to know when the Spotter was in the water to
% calculate the watch circle, in this script I will also get
% the deployment times semi-manually.
%
% This script uses data from location.csv files, which are generated
% by the Sofar's python script that parses the data for a Spotter
% (note this DOES NOT (!) include the pressure data from Smart Moorings).
%
% These longitude/latitudes are given at about once a minute.

clear
close all



%%

%
dataInfo.standardspotters.parentdir = ['/Volumes/GoogleDrive/Shared drives/ROXSI' ...
                                       '/LargeScale_Data_2022/RAW/Spotters/SDcards/'];
dataInfo.standardspotters.list_spotters = {'B01_spot1150', 'B01_spot1158', 'B03_spot1152', ...
                                           'B05_spot1153', ...
                                           'X01_spot1151', 'X03_spot1157', 'X04_spot1155'};

%
dataInfo.smartmoorings.parentdir = ['/Volumes/GoogleDrive/Shared drives/ROXSI' ...
                                    '/LargeScale_Data_2022/RAW/Spotters_Smart/SDcards/'];
dataInfo.smartmoorings.list_spotters = {'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
                                           'E07_spot1855', 'E07_spot1857', ...
                                           'E08_spot1852', ...
                                           'E09_spot1850', 'E09_spot1856', ...
                                           'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};


%% Deployment times for each Spotter and Smart Mooring
% (mannually obtained from the plot in the code block
% that comes later. These are rough, but may be useful
% for refining the deployment/recovery times)
%
% Times in local/PDT.

%
deplySpot.B01_spot1150 = [datenum(2022, 07, 06, 13, 00, 0), datenum(2022, 07, 20, 02, 20, 0)];
deplySpot.B01_spot1158 = [datenum(2022, 06, 15, 21, 50, 0), datenum(2022, 07, 02, 22, 15, 0)];
%
deplySpot.B03_spot1152 = [datenum(2022, 06, 16, 01, 55, 0), datenum(2022, 07, 19, 12, 07, 0)];
deplySpot.B05_spot1153 = [datenum(2022, 06, 15, 10, 00, 0), datenum(2022, 07, 19, 22, 00, 0)];
%
deplySpot.X01_spot1151 = [datenum(2022, 06, 15, 19, 20, 0), datenum(2022, 07, 19, 22, 00, 0)];
deplySpot.X03_spot1157 = [datenum(2022, 06, 15, 18, 00, 0), datenum(2022, 07, 20, 06, 00, 0)];
deplySpot.X04_spot1155 = [datenum(2022, 06, 15, 13, 00, 0), datenum(2022, 07, 19, 08, 40, 0)];

%
deplySpot.E01_spot1851 = [datenum(2022, 06, 17, 12, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
deplySpot.E02_spot1859 = [datenum(2022, 06, 17, 14, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
deplySpot.E05_spot1853 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E07_spot1855 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 05, 05, 00, 0)];
deplySpot.E07_spot1857 = [datenum(2022, 07, 05, 13, 40, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E08_spot1852 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E09_spot1850 = [datenum(2022, 06, 17, 09, 00, 0), datenum(2022, 06, 24, 05, 00, 0)];
deplySpot.E09_spot1856 = [datenum(2022, 06, 24, 13, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E10_spot1848 = [datenum(2022, 06, 17, 14, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
deplySpot.E11_spot1860 = [datenum(2022, 06, 17, 15, 20, 0), datenum(2022, 07, 21, 03, 00, 0)];    % recovered a day later after the rest of the array
deplySpot.E13_spot1849 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];


%%

%
list_fields_type = fieldnames(dataInfo);

% Loop over types of Spotters
for i1 = 1:length(list_fields_type)

    %
    Nspotters = length(dataInfo.(list_fields_type{i1}).list_spotters);

    % Loop over Spotters
    for i2 = 1:Nspotters

        % ----------------------------------------------------
        %
        locationfile_aux = fullfile(dataInfo.(list_fields_type{i1}).parentdir, ...
                                    dataInfo.(list_fields_type{i1}).list_spotters{i2}, ...
                                    'parsed', 'location.csv');

        %
        fid = fopen(locationfile_aux);
    
        %
        alldata = textscan(fid, '%d%d%d %d%d%d %d %f %f', ...
                                'Delimiter', ',', 'HeaderLines', 1);
    
    
        % Parse the time and location data    
        spotdata.dtime = datenum(double(alldata{1}), double(alldata{2}), double(alldata{3}), ...
                                 double(alldata{4}), double(alldata{5}), double(alldata{6} + (alldata{7}./1000)));
        %
        spotdata.latitude = alldata{8};
        spotdata.longitude = alldata{9};
    
        % Convert time from UTC to PDT
        spotdata.dtime = spotdata.dtime - (7/24);
    
        % ----------------------------------------------------
        % First eliminate unreasonable points that were tests
        % in the lab and/or by Sofar
        
        lin_timeoverall = (spotdata.dtime >= datenum(2022, 06, 14)) & ...
                          (spotdata.dtime <= datenum(2022, 07, 30));
        %
        spotdata.dtime = spotdata.dtime(lin_timeoverall);
        spotdata.latitude = spotdata.latitude(lin_timeoverall);
        spotdata.longitude = spotdata.longitude(lin_timeoverall);
    
        % ----------------------------------------------------
        % Compute horizontal speed
        arc_len_aux = distance(spotdata.latitude(1:end-1), spotdata.longitude(1:end-1), ...
                               spotdata.latitude(2:end), spotdata.longitude(2:end));
        dist_aux = deg2km(arc_len_aux);
    % %     %
    % %     horiz_speed_aux = dist_aux./(24*3600*diff(spotdata.dtime));
    
    % %     %
    % %     [ARCLEN, AZ] = distance(LAT1,LON1,LAT2,LON2);
    
        % ----------------------------------------------------
        % Simple figure to check the deployment start and end times
        if ~exist("deplySpot", 'var')
            figure
                %
                set(gcf, 'units', 'normalized')
                set(gcf, 'Position', [0.3320, 0.2965, 0.4328, 0.3014])
        
                %
                haxs_1 = subplot(1, 2, 1);
                    plot(haxs_1, spotdata.longitude, spotdata.latitude, '.k')
                %
                haxs_2 = subplot(1, 2, 2);
                    plot(haxs_2, spotdata.dtime, spotdata.longitude, '.-')
        
                %
                set([haxs_1, haxs_2], 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
        
                %
                title(haxs_1, dataInfo.(list_fields_type{i1}).list_spotters{i2}(1:3), 'FontSize', 20)
                title(haxs_2, dataInfo.(list_fields_type{i1}).list_spotters{i2}(1:3), 'FontSize', 20)
            %
            disp('Check the data!!')
            keyboard
        end

        % ----------------------------------------------------
        % With start/end deployment defined, compute watch
        % circle and plot the result
        if exist("deplySpot", 'var')
            %
            lindeply_aux = (spotdata.dtime >= deplySpot.(dataInfo.(list_fields_type{i1}).list_spotters{i2})(1)) & ...
                           (spotdata.dtime <= deplySpot.(dataInfo.(list_fields_type{i1}).list_spotters{i2})(2));
            % 
            [lat0, lon0, r0] = Spotter_fitwatchcircle(spotdata.latitude(lindeply_aux), ...
                                                      spotdata.longitude(lindeply_aux));

            % Compute circle to plot
            [easting0, northing0, UTMZone] = lltoUTM(lat0, lon0);
            %
            circumplt.easting = easting0 + r0.*cos(linspace(0, 2*pi, 50));
            circumplt.northing = northing0 + r0.*sin(linspace(0, 2*pi, 50));
            %
            [circumplt.latitude, circumplt.longitude] = UTMtoll(circumplt.northing, circumplt.easting, UTMZone);

            %
            newFigDims([12.17, 7.3])
                %
                plot(spotdata.longitude(lindeply_aux), spotdata.latitude(lindeply_aux), '.k')
                hold on
  
                %
                plot(mean(spotdata.longitude(lindeply_aux)), ...
                     mean(spotdata.latitude(lindeply_aux)), ...
                     '.', 'Color', [0.301, 0.745, 0.933], 'MarkerSize', 62)

                %
                plot(circumplt.longitude, ...
                     circumplt.latitude, 'r', 'LineWidth', 2)
                %
                plot(lon0, lat0, 'pr', 'MarkerSize', 42, 'MarkerFaceColor', 'r')
            
                %
                hleg = legend('spotter location', 'average location', ...
                              ['circle fit (r = ' num2str(r0, '%.0f') ' m)'], ...
                              'circle center', 'Location', 'EastOutside');
                    hleg.FontSize = 24;
            
                %
                axis equal
            
                %
                grid on
                set(gca, 'FontSize', 16, 'DataAspectRatio', [1, 1, 1])
                %
                xlim(lon0 + (r0 * 1.8e-5.*[-1, 1]))
                ylim(lat0 + (r0 * 1.8e-5.*[-1, 1]))
            
                %
                xlabel('Longitude', 'Interpreter', 'Latex', 'FontSize', 22)
                ylabel('Latitude', 'Interpreter', 'Latex', 'FontSize', 22)
                %
                title([dataInfo.(list_fields_type{i1}).list_spotters{i2}(1:3) ' - SN ' ...
                       dataInfo.(list_fields_type{i1}).list_spotters{i2}(9:end) ...
                      ': estimate of Spotter center'], 'Interpreter', 'Latex', 'FontSize', 22)

            % ----------------------------------------------------
            % Save fit reference location / watchcircle radius
            spottersDeployment.(dataInfo.(list_fields_type{i1}).list_spotters{i2}).rough_time_lims = deplySpot.(dataInfo.(list_fields_type{i1}).list_spotters{i2});
            %
            spottersDeployment.(dataInfo.(list_fields_type{i1}).list_spotters{i2}).latitude = lat0;
            spottersDeployment.(dataInfo.(list_fields_type{i1}).list_spotters{i2}).longitude = lon0;
            %
            spottersDeployment.(dataInfo.(list_fields_type{i1}).list_spotters{i2}).radius_watchcircle = r0;


        end

    end
end


%% Print to the screen 

%
fieldbla = {'longitude'};

% % %
% % [spottersDeployment.X01_spot1151.(fieldbla{1}); ...
% %  spottersDeployment.X03_spot1157.(fieldbla{1}); ...
% %  spottersDeployment.X04_spot1155.(fieldbla{1}); ...
% %  spottersDeployment.B01_spot1150.(fieldbla{1}); ...
% %  spottersDeployment.B03_spot1152.(fieldbla{1}); ...
% %  spottersDeployment.B05_spot1153.(fieldbla{1}); ...
% %  spottersDeployment.E01_spot1851.(fieldbla{1}); ...
% %  spottersDeployment.E02_spot1859.(fieldbla{1}); ...
% %  spottersDeployment.E05_spot1853.(fieldbla{1}); ...
% %  spottersDeployment.E07_spot1855.(fieldbla{1}); ...
% %  spottersDeployment.E08_spot1852.(fieldbla{1}); ...
% %  spottersDeployment.E09_spot1850.(fieldbla{1}); ...
% %  spottersDeployment.E10_spot1848.(fieldbla{1}); ...
% %  spottersDeployment.E11_spot1860.(fieldbla{1}); ...
% %  spottersDeployment.E13_spot1849.(fieldbla{1})]

%
[sprintf('%.6f', spottersDeployment.X01_spot1151.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.X03_spot1157.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.X04_spot1155.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.B01_spot1150.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.B03_spot1152.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.B05_spot1153.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E01_spot1851.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E02_spot1859.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E05_spot1853.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E07_spot1855.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E08_spot1852.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E09_spot1850.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E10_spot1848.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E11_spot1860.(fieldbla{1})); ...
 sprintf('%.6f', spottersDeployment.E13_spot1849.(fieldbla{1}))]

