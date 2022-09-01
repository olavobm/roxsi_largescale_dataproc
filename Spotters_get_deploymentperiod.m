%% Get fine-tuned deployment start and end times
% for all the Spotters (standard and smart moorings)
%
% In other words, based on deployment and recovery days, plot different
% variables to get a fairly precise time when Spotters were deployed
% and recovery. By the way, all Spotters were NOT turned on and off
% in a consistent manner, so that's not a great approach. In particular,
% some do not have a lot of data before/after deployment to make it obvious
% when instrument went in and out of the water.
%
% This script creates a table with begin and end times for each Spotter.
% There are deployment and recovery times (just before Spotter going in
% the water, and after Spotter is recovered), and trim times (with Spotters
% in the water, removing a few hours at the start of Smart mooring until
% after divers checked/moved Smart moorings, and after SST sensor has
% "equilibrated" for standard Spotters (or maybe a boat effect on SST
% that takes several minutes to disperse???). SST can take several minutes
% after Spotter is in the water to be "good"/trustworthy.

% % clear
% % close all


%%
% -------------------------------------------------------
% -------------------------------------------------------
% -------------------------------------------------------


%%

%
dataInfo.standardspotters.parentdir = ['/Volumes/GoogleDrive/Shared drives/ROXSI' ...
                                       '/LargeScale_Data_2022/RAW/Spotters/SDcards/'];
dataInfo.standardspotters.list_spotters = {'B01_spot1158', 'B01_spot1150', 'B03_spot1152', ...
                                           'B05_spot1153', ...
                                           'X01_spot1151', 'X03_spot1157', 'X04_spot1155'};

%
% % dataInfo.smartmoorings.parentdir = ['/Volumes/GoogleDrive/Shared drives/ROXSI' ...
% %                                     '/LargeScale_Data_2022/RAW/Spotters_Smart/SDcards/'];
dataInfo.smartmoorings.parentdir = ['/Volumes/LaCie/RAW/Spotters_Smart/SDcards/'];

%
dataInfo.smartmoorings.list_spotters = {'E01_spot1851', 'E02_spot1859', 'E05_spot1853', ...
                                        'E07_spot1855', 'E07_spot1857', ...
                                        'E08_spot1852', ...
                                        'E09_spot1850', 'E09_spot1856', ...
                                        'E10_spot1848', 'E11_spot1860', 'E13_spot1849'};


%% DAYS (with only a reference time) of when each Spotter
% was deployed and recovered
%
% Time in local/PDT.

% ------------------------------------------------------------
%
deplySpot.B01_spot1158 = [datenum(2022, 06, 15, 21, 50, 0), datenum(2022, 07, 02, 22, 15, 0)];
deplySpot.B01_spot1150 = [datenum(2022, 07, 06, 13, 00, 0), datenum(2022, 07, 20, 02, 20, 0)];
%
deplySpot.B03_spot1152 = [datenum(2022, 06, 16, 01, 55, 0), datenum(2022, 07, 19, 12, 07, 0)];
deplySpot.B05_spot1153 = [datenum(2022, 06, 15, 10, 00, 0), datenum(2022, 07, 19, 22, 00, 0)];
%
deplySpot.X01_spot1151 = [datenum(2022, 06, 15, 19, 20, 0), datenum(2022, 07, 19, 22, 00, 0)];
deplySpot.X03_spot1157 = [datenum(2022, 06, 15, 18, 00, 0), datenum(2022, 07, 20, 06, 00, 0)];
deplySpot.X04_spot1155 = [datenum(2022, 06, 15, 13, 00, 0), datenum(2022, 07, 19, 08, 40, 0)];

% ------------------------------------------------------------
%
deplySpot.E01_spot1851 = [datenum(2022, 06, 17, 12, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
deplySpot.E02_spot1859 = [datenum(2022, 06, 17, 14, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
deplySpot.E05_spot1853 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E07_spot1855 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 05, 08, 00, 0)];
deplySpot.E07_spot1857 = [datenum(2022, 07, 05, 13, 40, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E08_spot1852 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E09_spot1850 = [datenum(2022, 06, 17, 09, 00, 0), datenum(2022, 06, 24, 05, 00, 0)];
deplySpot.E09_spot1856 = [datenum(2022, 06, 24, 13, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
%
deplySpot.E10_spot1848 = [datenum(2022, 06, 17, 14, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];
deplySpot.E11_spot1860 = [datenum(2022, 06, 17, 15, 20, 0), datenum(2022, 07, 21, 06, 00, 0)];    % recovered a day later after the rest of the array
deplySpot.E13_spot1849 = [datenum(2022, 06, 17, 15, 00, 0), datenum(2022, 07, 20, 04, 00, 0)];



%%
% -------------------------------------------------------
% ------- FOR STANDARD SPOTTERS, LOOK AT COORDINATES, ---
% ------------------ SST, DISPLACEMENT ------------------
% -------------------------------------------------------

%
Nspotters = length(dataInfo.standardspotters.list_spotters);


% Loop over Spotters in either Standard or Smart Mooring
tic
for i1 = [1:5, 7]% 6%3:Nspotters

    % ----------------------------------------------------
    %
    locationfile_aux = fullfile(dataInfo.standardspotters.parentdir, ...
                                dataInfo.standardspotters.list_spotters{i1}, ...
                                'parsed', 'location.csv');

    %
    fid = fopen(locationfile_aux);

    %
    alldata = textscan(fid, '%d%d%d %d%d%d %d %f %f', ...
                            'Delimiter', ',', 'HeaderLines', 1);


    % Parse the time and location data -- MAKE SURE THAT EACH
    % ARE INVIDIVIDUALLY CONVERTED TO DOUBLE!!!  
    spotdata.dtime = datenum(double(alldata{1}), double(alldata{2}), double(alldata{3}), ...
                             double(alldata{4}), double(alldata{5}), double(alldata{6}) + (double(alldata{7})./1000));
    %
    spotdata.latitude = alldata{8};
    spotdata.longitude = alldata{9};

    % Convert time from UTC to PDT
    spotdata.dtime = spotdata.dtime - (7/24);

    %
    fclose(fid);


    % ----------------------------------------------------
    % Load different variable
    % 
    % Displacement can also be done for Smart mooring

    %
    displacementfile_aux = fullfile(dataInfo.standardspotters.parentdir, ...
                                    dataInfo.standardspotters.list_spotters{i1}, ...
                                    'parsed', 'displacement.csv');
    %
    fid = fopen(displacementfile_aux);

    %
    alldata = textscan(fid, '%d%d%d %d%d%d %d %f%f%f', ...
                            'Delimiter', ',', 'HeaderLines', 1);


    % Parse the time and location data -- MAKE SURE THAT EACH
    % ARE INVIDIVIDUALLY CONVERTED TO DOUBLE!!!
    spotdata.displacement.dtime = datenum(double(alldata{1}), double(alldata{2}), double(alldata{3}), ...
                                          double(alldata{4}), double(alldata{5}), double(alldata{6}) + (double(alldata{7})./1000));

    %
    spotdata.displacement.x = alldata{8};
    spotdata.displacement.y = alldata{9};
    spotdata.displacement.z = alldata{10};
    
    % Convert time from UTC to PDT
    spotdata.displacement.dtime = spotdata.displacement.dtime - (7/24);
    %
    fclose(fid);


    % ------------
    % Now SST
    %
    sstfile_aux = fullfile(dataInfo.standardspotters.parentdir, ...
                           dataInfo.standardspotters.list_spotters{i1}, ...
                           'parsed', 'sst.csv');
    %
    fid = fopen(sstfile_aux);

    %
    alldata = textscan(fid, '%d%d%d %d%d%d %d %f', ...
                            'Delimiter', ',', 'HeaderLines', 1);


    % Parse the time and location data -- MAKE SURE THAT EACH
    % ARE INVIDIVIDUALLY CONVERTED TO DOUBLE!!!    
    spotdata.sst.dtime = datenum(double(alldata{1}), double(alldata{2}), double(alldata{3}), ...
                                 double(alldata{4}), double(alldata{5}), double(alldata{6}) + (double(alldata{7})./1000));
    %
    spotdata.sst.SST = alldata{8};
    
    % Convert time from UTC to PDT
    spotdata.sst.dtime = spotdata.sst.dtime - (7/24);

    %
    fclose(fid);


    % ----------------------------------------------------
    % Plot variables to identify deployment/recovery

    %
    time_pltlims_deployment = [(deplySpot.(dataInfo.standardspotters.list_spotters{i1})(1) - 1), ...
                               deplySpot.(dataInfo.standardspotters.list_spotters{i1})(1)];
    %
    time_pltlims_recovery = [deplySpot.(dataInfo.standardspotters.list_spotters{i1})(2), ...
                             (deplySpot.(dataInfo.standardspotters.list_spotters{i1})(2) + (5/24))];

    %
    figure
        %
        set(gcf, 'units', 'normalized')
        set(gcf, 'Position', [0.25, 0.25, 0.5, 0.4])
        %
        haxs_1 = axes('Position', [0.1, 0.6, 0.35, 0.33]);
        haxs_2 = axes('Position', [0.55, 0.6, 0.35, 0.35]);
        haxs_3 = axes('Position', [0.1, 0.1, 0.35, 0.35]);
        haxs_4 = axes('Position', [0.55, 0.1, 0.35, 0.35]);

        %
        plot(haxs_1, datetime(spotdata.displacement.dtime, 'ConvertFrom', 'datenum'), spotdata.displacement.z, '.-')
        plot(haxs_2, datetime(spotdata.displacement.dtime, 'ConvertFrom', 'datenum'), spotdata.displacement.z, '.-')
        %
        plot(haxs_3, datetime(spotdata.sst.dtime, 'ConvertFrom', 'datenum'), spotdata.sst.SST, '.-')
        plot(haxs_4, datetime(spotdata.sst.dtime, 'ConvertFrom', 'datenum'), spotdata.sst.SST, '.-')

    %
    set([haxs_1, haxs_2, haxs_3, haxs_4], 'FontSize', 16, 'Box', 'on', ...
                                          'XGrid', 'on', 'YGrid', 'on')
    set([haxs_1, haxs_3], 'XLim', datetime(time_pltlims_deployment, 'ConvertFrom', 'datenum'))
    set([haxs_2, haxs_4], 'XLim', datetime(time_pltlims_recovery, 'ConvertFrom', 'datenum'))

    %
    xlabel(haxs_3, 'Time [local/PDT]', 'Interpreter', 'Latex', 'FontSize', 20)
    xlabel(haxs_4, 'Time [local/PDT]', 'Interpreter', 'Latex', 'FontSize', 20)
    %
    ylabel(haxs_1, '[m]', 'Interpreter', 'Latex', 'FontSize', 20)
    ylabel(haxs_3, '[$^\circ$C]', 'Interpreter', 'Latex', 'FontSize', 20)

    %
    title(haxs_1, ['Spotter: ' dataInfo.standardspotters.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.standardspotters.list_spotters{i1}(9:end) ' at deployment, z displacement (m)'], ...
                   'Interpreter', 'Latex', 'FontSize', 20)
    title(haxs_2, ['Spotter: ' dataInfo.standardspotters.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.standardspotters.list_spotters{i1}(9:end) ' at recovery, z displacement (m)'], ...
                   'Interpreter', 'Latex', 'FontSize', 18)
    %
    title(haxs_3, ['Spotter: ' dataInfo.standardspotters.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.standardspotters.list_spotters{i1}(9:end) ' at deployment, SST ($^\circ$C)'], ...
                   'Interpreter', 'Latex', 'FontSize', 20)
    title(haxs_4, ['Spotter: ' dataInfo.standardspotters.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.standardspotters.list_spotters{i1}(9:end) ' at recovery, SST ($^\circ$C)'], ...
                   'Interpreter', 'Latex', 'FontSize', 18)


    %
    linkaxes([haxs_1, haxs_3], 'x')
    linkaxes([haxs_2, haxs_4], 'x')

    %
    toc
    %
    
end


%%
return

%%
% ------------------------------------------------------
% ------ SIMILAR AS ABOVE, BUT FOR SMART MOORINGS ------
% ------------------------------------------------------

%%

%
Nspotters = length(dataInfo.smartmoorings.list_spotters);

% Loop over Spotters in either Standard or Smart Mooring
tic
for i1 = 1:Nspotters

    % ----------------------------------------------------
    % Load lon/lat
    locationfile_aux = fullfile(dataInfo.smartmoorings.parentdir, ...
                                dataInfo.smartmoorings.list_spotters{i1}, ...
                                'parsed', 'location.csv');

    %
    fid = fopen(locationfile_aux);

    %
    alldata = textscan(fid, '%d%d%d %d%d%d %d %f %f', ...
                            'Delimiter', ',', 'HeaderLines', 1);


    % Parse the time and location data -- MAKE SURE THAT EACH
    % ARE INVIDIVIDUALLY CONVERTED TO DOUBLE!!!
    spotdata.dtime = datenum(double(alldata{1}), double(alldata{2}), double(alldata{3}), ...
                             double(alldata{4}), double(alldata{5}), double(alldata{6}) + (double(alldata{7})./1000));
    %
    spotdata.latitude = alldata{8};
    spotdata.longitude = alldata{9};

    % Convert time from UTC to PDT
    spotdata.dtime = spotdata.dtime - (7/24);

    % ----------------------------------------------------
    % Load displacement from Spotter/Smart Mooring

    %
    displacementfile_aux = fullfile(dataInfo.smartmoorings.parentdir, ...
                                    dataInfo.smartmoorings.list_spotters{i1}, ...
                                    'parsed', 'displacement.csv');
    %
    fid = fopen(displacementfile_aux);

    %
    alldata = textscan(fid, '%d%d%d %d%d%d %d %f%f%f', ...
                            'Delimiter', ',', 'HeaderLines', 1);


    % Parse the time and location data -- MAKE SURE THAT EACH
    % ARE INVIDIVIDUALLY CONVERTED TO DOUBLE!!!
    spotdata.displacement.dtime = datenum(double(alldata{1}), double(alldata{2}), double(alldata{3}), ...
                                          double(alldata{4}), double(alldata{5}), double(alldata{6}) + (double(alldata{7})./1000));

    %
    spotdata.displacement.x = alldata{8};
    spotdata.displacement.y = alldata{9};
    spotdata.displacement.z = alldata{10};
    
    % Convert time from UTC to PDT
    spotdata.displacement.dtime = spotdata.displacement.dtime - (7/24);
    %
    fclose(fid);


    % ----------------------------------------------------
    % Load pressure from Smart Mooring
    %
    spotter_pres = Spotter_read_SMD_allSDcard(fullfile(dataInfo.smartmoorings.parentdir, ...
                                                       dataInfo.smartmoorings.list_spotters{i1}));

    % Put relevant variables in the data structure
    % and only feep valid data
%     lkeep_aux = (spotter_pres.allfiles.unixEpoch~=0) & (spotter_pres.allfiles.link==1);
    lkeep_aux = (spotter_pres.allfiles.link==1);
    %
    spotdata.pdata.dtime = spotter_pres.allfiles.dtime(lkeep_aux) - (7/24);
    spotdata.pdata.pressure = spotter_pres.allfiles.pressure(lkeep_aux);

    
    % ----------------------------------------------------
    % Make figure around deployment and recovery time

    %
    time_pltlims_deployment = [(deplySpot.(dataInfo.smartmoorings.list_spotters{i1})(1) - 1), ...
                                deplySpot.(dataInfo.smartmoorings.list_spotters{i1})(1)];
    %
    time_pltlims_recovery = [deplySpot.(dataInfo.smartmoorings.list_spotters{i1})(2), ...
                            (deplySpot.(dataInfo.smartmoorings.list_spotters{i1})(2) + (5/24))];


    %
    figure
        %
        set(gcf, 'units', 'normalized')
        set(gcf, 'Position', [0.25, 0.25, 0.5, 0.4])
        %
        haxs_1 = axes('Position', [0.1, 0.6, 0.35, 0.325]);
        haxs_2 = axes('Position', [0.55, 0.6, 0.35, 0.325]);
        haxs_3 = axes('Position', [0.1, 0.1, 0.35, 0.325]);
        haxs_4 = axes('Position', [0.55, 0.1, 0.35, 0.325]);

        %
        plot(haxs_1, datetime(spotdata.displacement.dtime, 'ConvertFrom', 'datenum'), spotdata.displacement.z, '.-')
        plot(haxs_2, datetime(spotdata.displacement.dtime, 'ConvertFrom', 'datenum'), spotdata.displacement.x, '.-')
        %
        plot(haxs_3, datetime(spotdata.pdata.dtime, 'ConvertFrom', 'datenum'), spotdata.pdata.pressure, '.-')
        plot(haxs_4, datetime(spotdata.pdata.dtime, 'ConvertFrom', 'datenum'), spotdata.pdata.pressure, '.-')

    %
    set([haxs_1, haxs_2, haxs_3, haxs_4], 'FontSize', 16, 'Box', 'on', ...
                                          'XGrid', 'on', 'YGrid', 'on')
    set([haxs_1, haxs_3], 'XLim', datetime(time_pltlims_deployment, 'ConvertFrom', 'datenum'))
    set([haxs_2, haxs_4], 'XLim', datetime(time_pltlims_recovery, 'ConvertFrom', 'datenum'))

    %
    xlabel(haxs_3, 'Time [local/PDT]', 'Interpreter', 'Latex', 'FontSize', 20)
    xlabel(haxs_4, 'Time [local/PDT]', 'Interpreter', 'Latex', 'FontSize', 20)
    %
    ylabel(haxs_1, '[m]', 'Interpreter', 'Latex', 'FontSize', 20)
    ylabel(haxs_3, 'Pressure', 'Interpreter', 'Latex', 'FontSize', 20)

    %
    title(haxs_1, ['Smart mooring: ' dataInfo.smartmoorings.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.smartmoorings.list_spotters{i1}(9:end) ' z displacement around deployment'], ...
                   'Interpreter', 'Latex', 'FontSize', 20)
    title(haxs_2, ['Smart mooring: ' dataInfo.smartmoorings.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.smartmoorings.list_spotters{i1}(9:end) ' z displacement around recovery'], ...
                   'Interpreter', 'Latex', 'FontSize', 20)
    %
    title(haxs_3, ['Smart mooring: ' dataInfo.smartmoorings.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.smartmoorings.list_spotters{i1}(9:end) ' pressure around deployment'], ...
                   'Interpreter', 'Latex', 'FontSize', 20)
    title(haxs_4, ['Smart mooring: ' dataInfo.smartmoorings.list_spotters{i1}(1:3) ' ' ...
                   '- SN ' dataInfo.smartmoorings.list_spotters{i1}(9:end) ' pressure around recovery'], ...
                   'Interpreter', 'Latex', 'FontSize', 20)

    %
    linkaxes([haxs_1, haxs_3], 'x')
    linkaxes([haxs_2, haxs_4], 'x')

    %
    toc
    %
    pause(1)
end

%%
return

%%
% ------------------------------------------------------------------
% ------------------------------------------------------------------
% --------- CREATE DEPLOYMENT/RECOVERY TABLE FOR SPOTTERS ----------
% ------------------------------------------------------------------
% ------------------------------------------------------------------


%% Now create table with deployment recovery info
%
% All times are local.


%
mooringID = ["B01s"; "B01s"; ...
             "B03s"; "B05s"; "X01s"; "X03s"; "X04s"; ...
             "E01sp"; "E02sp"; "E05sp"; ...
             "E07sp"; "E07sp"; ...
             "E08sp"; ...
             "E09sp"; "E09sp"; ...
             "E10sp"; "E11sp"; "E13sp"];
%
SN = ["1158"; "1150"; ...
                  "1152"; "1153"; "1151"; "1157"; "1155"; ...
                  "1851"; "1859"; "1853"; ...
                  "1855"; "1857"; ...
                  "1852"; ...
                  "1850"; "1856"; ...
                  "1848"; "1860"; "1849"];

% String name as in the mooring table to be consistent
spottertype = [repmat("spotter", 7, 1); ...
               repmat("smartspotter", 11, 1)];



%% Define being/end deployment and trim times


% For the standard spotters neither displacement nor SST give
% a clear signal of when it was pulled out of the water.
% Lat/long can be used (and maybe that's what I should have done)
% but a few Spotters or Smart moorings do not have measurements
% after recovery.
%
% Message on WhatsApp from Charlotte said that they were coming back with
% all 6 standard spotters at 09:53am. 6 standard Spotters were recovered
% between 07:30 and 09:53 (where 07:30 is about the time when the NPS
% boat arrived at the dock to offload instruments before departing to
% get the Spotters). I'll set trimming time at 07:30 on 07/20.

%
% - E01 1851: weird problem at the beginning, with a 4.5-hour gap in
%             pressure data. In the SMD file, the character RBRD, turned
%             into RBRU during the gap. This is the only smart mooring
%             with this issue.


% ------------------------------------------
%
time_begin_deployment = ["2022/06/15 08:41:00"; ...    % Standard spotter / B01 1158 -- issues needed to be recovered
                         "2022/07/06 08:59:30"; ...    % B01 1150
                         "2022/06/15 08:47:00"; ...    % B03 1152
                         "2022/06/15 08:55:00"; ...    % B05 1153
                         "2022/06/15 09:00:00"; ...    % X01 1151
                         "2022/06/15 09:12:00"; ...    % X03 1157
                         "2022/06/15 09:17:00"; ...    % X04 1155
                         "2022/06/17 05:54:00"; ...    % Smart mooring / E01 1851
                         "2022/06/17 06:01:30"; ...    % E02 1859
                         "2022/06/17 06:13:00"; ...    % E05 1853 -- this one stopped working at once.
                         "2022/06/17 08:12:30"; ...    % E07 1855 -- spotter that went bad and was replaced.
                         "2022/07/05 10:20:00"; ...     % E07 1857 -- this was a replacement Spotter.
                         "2022/06/17 08:26:00"; ...    % E08 1852
                         "2022/06/17 06:23:00"; ...    % E09 1850
                         "2022/06/24 07:53:15"; ...    % E09 1856
                         "2022/06/17 06:36:00"; ...    % E10 1848
                         "2022/06/17 08:32:30"; ...    % E11 1860
                         "2022/06/17 08:05:00"];    % E13 1849

time_end_deployment = ["2022/07/06 09:00:00"; ...    % Standard spotter / B01 1158 -- issues needed to be recovered. Time based on the deployment of the other one because 1158 stopped recording (solar power/chargin failed)
                       "2022/07/20 10:00:00"; ...    % B01 1150
                       "2022/07/20 10:00:00"; ...    % B03 1152
                       "2022/07/20 10:00:00"; ...    % B05 1153
                       "2022/07/20 10:00:00"; ...    % X01 1151
                       "2022/07/20 10:00:00"; ...    % X03 1157
                       "2022/07/20 10:00:00"; ...    % X04 1155
                       "2022/07/20 06:12:00"; ...    % Smart mooring / E01 1851
                       "2022/07/20 06:15:00"; ...     % E02 1859
                       "2022/07/17 07:10:00"; ...    % E05 1853 -- this one stopped working at once. deployment uncertain, unless there is location data
                       "2022/07/05 10:35:00"; ...     % E07 1855 -- spotter that went bad and was replaced.
                       "2022/07/20 06:30:00"; ...     % E07 1857 -- this was a replacement Spotter.
                       "2022/07/20 06:55:00"; ...     % E08 1852
                       "2022/06/24 07:52:00"; ...    % E09 1850 
                       "2022/07/20 07:00:00"; ...    % E09 1856
                       "2022/07/20 07:10:00"; ...    % E10 1848
                       "2022/07/21 09:10:00"; ...    % E11 1860
                       "2022/07/20 06:40:00"];    % E13 1849

% ------------------------------------------
%
time_begin_trim = ["2022/06/15 09:00:00"; ...    % Standard spotter / B01 1158 -- issues needed to be recovered
                   "2022/07/06 09:10:00"; ...    % B01 1150
                   "2022/06/15 09:05:00"; ...    % B03 1152
                   "2022/06/15 09:10:00"; ...    % B05 1153
                   "2022/06/15 09:20:00"; ...    % X01 1151
                   "2022/06/15 09:30:00"; ...    % X03 1157
                   "2022/06/15 09:35:00"; ...    % X04 1155
                   "2022/06/17 10:39:30"; ...    % Smart mooring / E01 1851
                   "2022/06/17 07:25:00"; ...    % E02 1859    06:11:30 anchor on the bottom  % not clear if divers moved it, maybe around 07:10
                   "2022/06/17 07:35:00"; ...    % E05 1853 -- this one stopped working at once.  06:20:30 anchor on the bottom
                   "2022/06/17 08:40:00"; ...    % E07 1855 -- spotter that went bad and was replaced. 08:23:30 anchor on the bottom
                   "2022/07/05 10:45:00"; ...     % E07 1857 -- 10:26:00 pressure sensor on the bottom, but was putting this in place of the other one by divers (i.e. this was a replacement Spotter, and moorings were swaped at the anchor point by divers).
                   "2022/06/17 08:40:00"; ...    % E08 1852    08:29:00 anchor on the bottom    % doesn't look like divers moved this one
                   "2022/06/17 08:08:00"; ...    % E09 1850    06:33:15 anchor on the bottom
                   "2022/06/24 08:56:00"; ...    % E09 1856    07:59:00 anchor on the bottom
                   "2022/06/17 08:20:00"; ...    % E10 1848    06:39:30 anchor on the bottom
                   "2022/06/17 08:50:00"; ...    % E11 1860 08:35:30 anchor on the bottom
                   "2022/06/17 08:32:00"];    % E13 1849,   08:09:30 anchor on the bottom

%
time_end_trim = ["2022/07/02 01:22:00"; ...    % Standard spotter / B01 1158 -- issues needed to be recovered. This end trim time may be refined after looking at data and issues that happened towards the end.
                 "2022/07/20 07:30:00"; ...    % B01 1150
                 "2022/07/20 07:30:00"; ...    % B03 1152
                 "2022/07/20 07:30:00"; ...    % B05 1153
                 "2022/07/20 07:30:00"; ...    % X01 1151
                 "2022/07/20 07:30:00"; ...    % X03 1157
                 "2022/07/20 07:30:00"; ...    % X04 1155
                 "2022/07/20 06:00:00"; ...    % Smart mooring / E01 1851
                 "2022/07/20 06:05:00"; ...     % E02 1859
                 "2022/07/17 10:13:00"; ...    % E05 1853 -- this one stopped getting pressure at once at ~10:20:30. Trimming a few minutes before the data is visibly bad. Spotter buoy itself is getting measurements all the way too the end, so trimming for that can be done a couple hours before time_end_deployment.
                 "2022/07/05 10:00:00"; ...     % E07 1855 -- spotter that went bad and was replaced.
                 "2022/07/20 06:12:00"; ...     % E07 1857 
                 "2022/07/20 06:42:00"; ...     % E08 1852
                 "2022/06/24 07:40:00"; ...    % E09 1850 
                 "2022/07/20 06:52:00"; ...    % E09 1856
                 "2022/07/20 06:56:00"; ...    % E10 1848
                 "2022/07/21 08:55:00"; ...    % E11 1860
                 "2022/07/20 06:27:00"];    % E13 1849


%% Now create table with deployment information


%
dployInfo_Spotters = table(mooringID, SN, spottertype, ...
                           time_begin_deployment, time_end_deployment, ...
                           time_begin_trim, time_end_trim);


%% Save the table

% %
% save(fullfile(repo_dirpath(), ...
%               'deploymentInfo_Spotters_ROXSI2022.mat'), ...
%               'dployInfo_Spotters'))

