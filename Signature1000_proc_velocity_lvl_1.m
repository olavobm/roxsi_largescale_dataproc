%% Script that process ROXSI's large-scale array Signature1000
% data from RAW to level 1.
% 
%
% Signature1000 filenames have different formats
% (which is a bit problematic), but function Signature_orderfiles.m
% takes care of that..
%
% All Signatures (ADCPs) were programmed to start
% sampling at 2022/06/21 18:00:00 (local). Signature
% SN 103045, at B10, was only deployed later due to
% weather conditions.

clear
close all


%%
% --------------------------------------
% -------------- PREAMBLE --------------
% --------------------------------------

%%

roxsi_add_libraries()

%%

%
% % dirparent_data = data_dirpath();
dirparent_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
%
dir_rawdata_parent = fullfile(dirparent_data, 'RAW', 'Signature1000');


%%

%
% % dir_output_parent = data_dirpath();
% % dir_output_parent = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/';
% % dir_output_parent = '/home/omarques/Documents/MATLAB/ROXSIproc_output/';
dir_output_parent = pwd;
%
% dir_output_data_L1 = fullfile(dir_output_parent, 'Level1_Data', 'Signature1000_Level1');
% dir_output_figs_L1 = fullfile(dir_output_parent, 'Level1_Data', 'Signature1000_Level1', 'qc_plots');
dir_output_data_L1 = dir_output_parent;

% Logical switches to save or not save data and figures
lsave_file = true;
lsave_fig = true;


%% Define time limits and time step for processing velocity data. Since
% there are a lot of Signature1000, it seems appropriate to break
% it apart in smaller chunks of time. Roughly speaking,
% files with 1/4 of a day (6 hours) will be around 300 MB.
%
% Note that the pressure (and heading/tilt/pitch/temperature) are
% still loaded for the full timeseries (for general QC) and only
% the processing of velocity will be within these time limits.
%
% Note that in the definition I will use, each file will have data
% starting and including the first time limit of each segment
% (given by time_lims_proc) and ending just before (not including)
% the next time limit.
%
% Note these time limits are "real times" -- i.e. they should be
% compared to clock-drift-corrected timestamps.

% -----------------------------------
% % %
% % time_beginend_proc = [datetime(2022, 06, 29, 00, 00, 00), ...
% %                       datetime(2022, 06, 29, 06, 00, 00)];
% % time_beginend_proc.TimeZone = 'America/Los_Angeles';
% % 
% % %
% % dtstep_proc = hours(6);
% % 
% % %
% % time_lims_proc = time_beginend_proc(1) : dtstep_proc : time_beginend_proc(2);

% -----------------------------------
%
% % time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 06, 00, 00); ...
% %                   datetime(2022, 07, 03, 18, 00, 00), datetime(2022, 07, 04, 00, 00, 00)];
%
% % time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 06, 00, 00)];    % 6 hour chunk
%
time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 07, 08, 00, 00, 00)];   % around 10 days
%
% time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 00, 05, 00)];   % 5 minutes
% time_lims_proc = [datetime(2022, 06, 29, 00, 00, 00), datetime(2022, 06, 29, 00, 01, 00)];   % 1 minute
time_lims_proc.TimeZone = 'America/Los_Angeles';

%
Ndatasegments = size(time_lims_proc, 1);


%% Load ADCP deployment information

% Just to be clear: file and variable have
% the same name (though not a requirement)
%
% dir_coderepo = repo_dirpath();
dir_coderepo = '/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc';
load(fullfile(dir_coderepo, 'deploymentInfo_ROXSI2022.mat'), 'deploymentInfo_ROXSI2022')


%% List of Signatures that will be processed

% % All Signatures
% list_Signature = {'A01_103043', ...
%                   'B10_103045', ...
%                   'B13_103046', ...
%                   'B15_103056', ...
%                   'B17_101923', ...
%                   'C01_102128', ...
%                   'X05_100231', ...
%                   'X11_101941'};

% Just a test
% % list_Signature = {'A01_103043'};
list_Signature = {'B10_103045'};

%
Nsignatures = length(list_Signature);


%% The list of ADCPs processed with either 4 or 5 beams

%
% % list_4beams = {'B10', 'B13', 'B15', 'B17', 'X11'};    % B10 has 5th beam, but on HR mode (up to 4 meters above the transducer). Others did the echosounder.
% % list_5beams = {'A01', 'C01', 'X05'};
%
list_4beams = {'B10'};    % B10 has 5th beam, but on HR mode (up to 4 meters above the transducer). Others did the echosounder.
list_5beams = {'A01', 'B13', 'B15', 'B17', 'C01', 'X05', 'X11'};

% Check all of the list of Signatures is in either of the
% lists for processing with 4 or 5 beams
for i = 1:length(list_Signature)

    %
    linanylist = any(contains(list_4beams, list_Signature{i}(1:3))) || ...
                 any(contains(list_5beams, list_Signature{i}(1:3)));

    %
    if ~linanylist

        error(['Signature ' list_Signature{i} ' is not present in ' ...
               'any of the lists specifying velocity processing.'])

    end

end

    

%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
% -------------------------------------
% --------- LOAD ATM PRESSURE ---------
% -------------------------------------


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


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------------- DO DATA PROCESSING ----------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------


%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_Signature_procL1_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_output_data_L1, log_file_name))

%
totalRunTime = tic;

% % 
% % %% Fields to be deleters

%%

list_Data_fields = {'Burst_Time', ...
                    'Burst_VelBeam1', 'Burst_VelBeam2', 'Burst_VelBeam3', 'Burst_VelBeam4', ...
                    'Burst_AmpBeam1', 'Burst_AmpBeam2', 'Burst_AmpBeam3', 'Burst_AmpBeam4', ...
                    'Burst_Heading', 'Burst_Pitch', 'Burst_Roll', ...
                    'Burst_Temperature', 'Burst_Pressure'};


%% Display on the screen:

%
disp(' '), disp(' ')
disp('------------------------------ Processing data from Signature1000s ------------------------------')
disp('List of Signature1000s being processed:')
%
for i = 1:Nsignatures
    disp([num2str(i) ' - ' list_Signature{i}])
end



% Loop over Signature1000's in the list
for i1 = 1:Nsignatures


    %%
    %
    disp(' '), disp(' ')
    disp(['----- Start processing raw Signature1000 data: ' list_Signature{i1} ' -----'])


    %% Get the Signature1000 files in the correct order

    %
    dir_data_aux = fullfile(dir_rawdata_parent, list_Signature{i1}, 'converted');
    %
    list_dir_aux = Signature_orderfiles(dir_data_aux, list_Signature{i1}(1:3));

    %
    Nfiles = length(list_dir_aux);


    %% First add basic metadata to structure variable

    %
    sig1000.SN = convertCharsToStrings(list_Signature{i1}(5:end));
    sig1000.mooringID = convertCharsToStrings(list_Signature{i1}(1:3));

    % Latitude/longitude
    info_mooringtable = ROXSI_mooringlocation(sig1000.mooringID, "ADCP");
    %
    sig1000.latitude = info_mooringtable.latitude;
    sig1000.longitude = info_mooringtable.longitude;


    %% Get trimming times for the deployment period of this Signature

    %
    ind_row_match = find(strcmp(deploymentInfo_ROXSI2022.SN, list_Signature{i1}(5:end)));
    
    %
    time_1 = deploymentInfo_ROXSI2022.time_begin_trim(ind_row_match);
    time_2 = deploymentInfo_ROXSI2022.time_end_trim(ind_row_match);

    %
    time_1 = datetime(datenum(time_1, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    time_2 = datetime(datenum(time_2, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');

% %     %
% %     lin_deployment = (senAQDP_aux.time >= time_1) & (senAQDP_aux.time <= time_2);



    %% Define trimming edges as time processing edges

    %
    time_lims_proc = [time_1, time_2];
    %
    Ndatasegments = size(time_lims_proc, 1);


    %%
    % ---------------------------------------------------------------
    % -------- NOW GET HEADING AND TILT AND DO FIRST QC PLOT --------
    % ---------------------------------------------------------------

    %% Load Configuration from first file

    % Just the Config structure variable
    dataread_aux = load(fullfile(dir_data_aux, list_dir_aux(1).name), 'Config');

    %
    sig1000.Config = dataread_aux.Config;
    

    %% First get the pressure timeseries (and heading, pitch, and
    % roll only to ???? trimming??? clock drift???

% %     % Pre-allocate in cell arrays
% %     prealloc_aux = cell(1, Nfiles);
% %     %
% %     sig1000.timedatenum = prealloc_aux;
% %     sig1000.pressure = prealloc_aux;
% %     sig1000.temperature = prealloc_aux;
% %     %
% %     sig1000.heading = prealloc_aux;
% %     sig1000.pitch = prealloc_aux;
% %     sig1000.roll = prealloc_aux;

% %     %
% %     disp('--- Load all scalar variables ---')
% % 
% %     tic
% %     % Loop over files
% %     for i2 = 1:Nfiles

% %         %% Load data in file i2'th
% %         %
% %         dataread_aux = load(fullfile(dir_data_aux, list_dir_aux(i2).name));
% % 
% %         %% Delete variables...
% % 
% %         %
% %         dataread_aux = rmfield(dataread_aux, {'Units', 'Descriptions'});
% %         %
% %         list_fields_in_Data = fieldnames(dataread_aux.Data);
% %         %
% %         for i3 = 1:length(list_fields_in_Data)
% %             %
% %             if ~any(strcmp(list_fields_in_Data{i3}, list_Data_fields))
% %                 dataread_aux.Data = rmfield(dataread_aux.Data, list_fields_in_Data{i3});
% %             end
% %         end
% % 
% % 
% %         %% Put time, pressure (which is already in dbar),
% %         % heading, pitch, and roll in output data structure
% % 
% %         %
% %         sig1000.timedatenum{i2} = dataread_aux.Data.Burst_Time;
% %         sig1000.pressure{i2} = dataread_aux.Data.Burst_Pressure;
% %         sig1000.temperature{i2} = dataread_aux.Data.Burst_Temperature;
% %         %
% %         sig1000.heading{i2} = dataread_aux.Data.Burst_Heading;
% %         sig1000.pitch{i2} = dataread_aux.Data.Burst_Pitch;
% %         sig1000.roll{i2} = dataread_aux.Data.Burst_Roll;
% % 
% %         %%
% %         clear dataread_aux
% % 
% %     end
% %     disp('--- Done with getting all scalar variables ---')
% %     toc
% % 
% %     % Concatenate cell array into a long column vector
% %     sig1000.timedatenum = cat(1, sig1000.timedatenum{:});
% %     sig1000.pressure = cat(1, sig1000.pressure{:});
% %     sig1000.temperature = cat(1, sig1000.temperature{:});
% %     %
% %     sig1000.heading = cat(1, sig1000.heading{:});
% %     sig1000.pitch = cat(1, sig1000.pitch{:});
% %     sig1000.roll = cat(1, sig1000.roll{:});
% % 
% % 
% %     disp('--- Done getting timeseries of scalar variables ---')


    %% Convert time to date time

    %
% %     sig1000.dtime = datetime(sig1000.timedatenum, 'ConvertFrom', 'datenum');
% %     sig1000.dtime.TimeZone = 'America/Los_Angeles';
% %     disp('Done converting datenum to date time')


    %% As opposed to Aquadopps, Signatures don't have
    % a customizable pressure offset (i.e. the measured pressure is
    % the water pressure (+- atmospheric variability????).
    %
    % Apply correction due to atmospheric pressure variability.


    %% Calculate bottom depth assuming the pressure is hydrostatic

    %
    %sig1000.bottomdepthfrompres = 1e4*sig1000.pressure ./ (1030*9.8);


    %% Construct height above the bottom of bin centers
    % (assuming the ADCP is at the bottom)

    % In meters (based on the Solidworks drawing)
    sig1000.transducerHAB = (31.88)/100;

    % In meters
    sig1000.binsize = sig1000.Config.Burst_CellSize;
    
    % Height of the first cell center relative to transducer
    % (based on the Principles of Operation manual by Nortek, page 12)
    cellcenter_first_bin = sig1000.Config.Burst_BlankingDistance + ...
                           sig1000.Config.Burst_CellSize;

% %  I expected there was a factor of 1/2
% %     cellcenter_first_bin = sig1000.Config.Burst_BlankingDistance + ...
% %                                          (sig1000.Config.Burst_CellSize/2);

    %
    sig1000.cellcenter = cellcenter_first_bin + ...
                        (0:1:(double(sig1000.Config.Burst_NCells) - 1)) .* sig1000.binsize;

    %
    sig1000.zhab = sig1000.transducerHAB + sig1000.cellcenter;


    %% Find all bins that are entirely below the maximum bottom depth

% %     %
% %     lin_verticalrange = ((sig1000.zhab + (sig1000.binsize/2)) < ...
% %                          max(sig1000.bottomdepthfrompres));
    % Or dummy for code developing purposes
    lin_verticalrange = true(1, length(sig1000.zhab));

    %
    sig1000.zhab = sig1000.zhab(lin_verticalrange);
% % % 
% % % 
% % %     %%
% % %     % ---------------------------------------------------------------
% % %     % ------------- DO QC ON PRESSURE, HEADING, AND TILT ------------
% % %     % ---------------------------------------------------------------
% % % 
% % %     %%
% % % 
% % %     disp('--- Making first QC plot with timeseries of scalars ---')
% % % 
% % %     %
% % %     fig_L1_QC_tilt = figure;
% % %     set(fig_L1_QC_tilt, 'units', 'normalized')
% % %     set(fig_L1_QC_tilt, 'Position', [0.2, 0.2, 0.4, 0.6])
% % %         %
% % %         haxs_1 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.73, 0.8, 0.17]);
% % %         haxs_2 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.52, 0.8, 0.17]);
% % %         haxs_3 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.31, 0.8, 0.17]);
% % %         haxs_4 = axes(fig_L1_QC_tilt, 'Position', [0.1, 0.10, 0.8, 0.17]);
% % %         %
% % %         haxs_all = [haxs_1, haxs_2, haxs_3, haxs_4];
% % %         hold(haxs_all, 'on')
% % %         %
% % %         plot(haxs_1, sig1000.dtime, sig1000.pressure, '-k')
% % %         plot(haxs_2, sig1000.dtime, sig1000.heading, '-k')
% % %         plot(haxs_3, sig1000.dtime, sig1000.pitch, '-k')
% % %         plot(haxs_4, sig1000.dtime, sig1000.roll, '-k')
% % % 
% % %     %
% % %     set(haxs_all, 'FontSize', 16, 'Box', 'on', ...
% % %                   'XGrid', 'on', 'YGrid', 'on')
% % %     %
% % %     set(haxs_all, 'XLim', sig1000.dtime([1, end]) + [-hours(12); hours(12)])
% % %     %
% % %     lin_deployment = (sig1000.dtime >= time_1) & (sig1000.dtime <= time_2);
% % %     %
% % %     pres_sub_aux = sig1000.pressure(lin_deployment);
% % %     ylim(haxs_1, [min(pres_sub_aux), max(pres_sub_aux)])
% % % % %     ylim(haxs_2, [0, 360])
% % %     ylim(haxs_3, [-6, 6])
% % %     ylim(haxs_4, [-6, 6])
% % % 
% % % 
% % %     %
% % %     ylabel(haxs_1, '[dbar]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     ylabel(haxs_2, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     ylabel(haxs_3, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     %
% % %     title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
% % %                    char(sig1000.SN) ': pressure, heading, pitch, and roll'], ...
% % %                   'Interpreter', 'Latex', 'FontSize', 16)
% % %     %
% % %     linkaxes([haxs_1, haxs_2, haxs_3, haxs_4], 'x')
% % % 
% % % 
% % %     %
% % %     for i2 = 1:length(haxs_all)
% % %         ylims_aux = ylim(haxs_all(i2));
% % %         %
% % %         plot(haxs_all(i2), [time_1, time_1], ylims_aux, '--r')
% % %         plot(haxs_all(i2), [time_2, time_2], ylims_aux, '--r')
% % %         %
% % %         ylim(haxs_all(i2), ylims_aux)
% % %     end
% % % 
% % % 
% % %     %
% % %     exportgraphics(fig_L1_QC_tilt, fullfile(pwd, ['sig_pres_tilt_' list_Signature{i1} '.png']), 'Resolution', 300)
% % % 
% % %     %
% % %     pause(10)
% % % 
% % %     %
% % %     close(fig_L1_QC_tilt)


    %%

    %
    sig1000.coordsystem = 'ENU';

    % In clockwise degrees from the true north
    sig1000.magdec = 12.86;

    %%
    % ---------------------------------------------------------------
    % ----------------------- NOW GET VELOCITY ----------------------
    % ---------------------------------------------------------------

    %
    for i2 = 1:Ndatasegments

% %         % For B10 (deployed later than planned), only start
% %         % getting the files when there is data in the deployment
% %         if strcmp(list_Signature{i1}(1:3), 'B10')
% %             % If the latest time stamp if before the beginning of the B10
% %             % deployment then, skip and go to next file
% %             if max(dataread_aux.Burst_Time) < datenum(2022, 06, 24, 07, 38, 00)    % based on trimming at 40 min mark
% %                 %
% %                 warning(['Skipping file ' list_dir_aux(i2).name ' because ' ...
% %                          'the data is from before the deployment.'])
% %                 %
% %                 continue
% %             end 
% %         end


        % If the time segment is a time that is fully outside the
        % deployment trimming edges, than skip this time segment
        if (time_1 >= time_lims_proc(i2, 2)) || (time_2 < time_lims_proc(i2, 1))
            continue
        end 
       
        %
        time_lims_aux = datenum(time_lims_proc(i2, :));
        %
        [listfiles_perseg, struct_timelims] = Signature1000_filesintimelims(list_Signature{i1}(5:end), time_lims_aux);
        Nfilesperseg = length(listfiles_perseg);

        %
        disp(['--- Raw data that will processed is split into ' num2str(Nfilesperseg) ' data files. Loading data from ---'])
        listfiles_perseg

        %
        tic
        % ------------------------------------------------

        % Pre-allocate in cell arrays
        prealloc_aux = cell(Nfilesperseg, 1);
% % %         %
% % %         list_beam_vars = {'vel1', 'vel2', 'vel3', 'vel4', ...
% % %                           'amp1', 'amp2', 'amp3', 'amp4', ...
% % %                           'corr1', 'corr2', 'corr3', 'corr4'};

        % ------------------------------------------------
        sig1000.timedatenum = prealloc_aux;
        sig1000.pressure = prealloc_aux;
        sig1000.temperature = prealloc_aux;
        %
        sig1000.heading = prealloc_aux;
        sig1000.pitch = prealloc_aux;
        sig1000.roll = prealloc_aux;
        % ------------------------------------------------

        %
        sig1000.timednum_fourbeams = prealloc_aux;
        %
        list_beam_vars = {'vel1', 'vel2', 'vel3', 'vel4', ...
                          'amp1', 'amp2', 'amp3', 'amp4'};
        % Loop over variables
        for i3 = 1:length(list_beam_vars)
            sig1000.(list_beam_vars{i3}) = prealloc_aux;
        end
        %
        if any(contains(list_5beams, list_Signature{i1}(1:3)))
            sig1000.timednum_beam5 = prealloc_aux;
            sig1000.vel5_raw = prealloc_aux;
            sig1000.amp5_raw = prealloc_aux;
            sig1000.corr5_raw = prealloc_aux;
        end
        

        % Load over all files with data in the i2'th segment
        for i3 = 1:Nfilesperseg

            %
            dataread_aux = load(fullfile(dir_rawdata_parent, ...
                                         list_Signature{i1}, ...
                                         'converted', ...
                                         listfiles_perseg(i3)));

            %
            lin_proclims_beam4time_aux = (dataread_aux.Data.Burst_Time >= datenum(time_lims_proc(i2, 1))) & ...
                                         (dataread_aux.Data.Burst_Time <  datenum(time_lims_proc(i2, 2)));

            % ------------------------------------------------
            %  Get all scalars in the same way as velocity
            % (INSTEAD OF THE CODE HIGHER UP)
            %
            sig1000.timedatenum{i3} = dataread_aux.Data.Burst_Time(lin_proclims_beam4time_aux);
            sig1000.pressure{i3} = dataread_aux.Data.Burst_Pressure(lin_proclims_beam4time_aux);
            sig1000.temperature{i3} = dataread_aux.Data.Burst_Temperature(lin_proclims_beam4time_aux);
            %
            sig1000.heading{i3} = dataread_aux.Data.Burst_Heading(lin_proclims_beam4time_aux);
            sig1000.pitch{i3} = dataread_aux.Data.Burst_Pitch(lin_proclims_beam4time_aux);
            sig1000.roll{i3} = dataread_aux.Data.Burst_Roll(lin_proclims_beam4time_aux);
            % ------------------------------------------------


            % Dummy/for code development
            lin_verticalrange = true(1, size(dataread_aux.Data.Burst_VelBeam1, 2));

            % Get data from the 4 beams
            sig1000.vel1{i3} = dataread_aux.Data.Burst_VelBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sig1000.vel2{i3} = dataread_aux.Data.Burst_VelBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sig1000.vel3{i3} = dataread_aux.Data.Burst_VelBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sig1000.vel4{i3} = dataread_aux.Data.Burst_VelBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            %
            sig1000.amp1{i3} = dataread_aux.Data.Burst_AmpBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sig1000.amp2{i3} = dataread_aux.Data.Burst_AmpBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sig1000.amp3{i3} = dataread_aux.Data.Burst_AmpBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sig1000.amp4{i3} = dataread_aux.Data.Burst_AmpBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            %
            if any(strcmp(list_beam_vars, 'corr1'))
                sig1000.corr1{i3} = dataread_aux.Data.Burst_CorBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
                sig1000.corr2{i3} = dataread_aux.Data.Burst_CorBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
                sig1000.corr3{i3} = dataread_aux.Data.Burst_CorBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
                sig1000.corr4{i3} = dataread_aux.Data.Burst_CorBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            end

            %
            sig1000.timednum_fourbeams{i3} = dataread_aux.Data.Burst_Time(lin_proclims_beam4time_aux);

            % Get the 5th beam
            if any(contains(list_5beams, list_Signature{i1}(1:3)))
                %
                lin_proclims_beam5time_aux = (dataread_aux.Data.IBurst_Time >= datenum(time_lims_proc(i2, 1))) & ...
                                             (dataread_aux.Data.IBurst_Time <  datenum(time_lims_proc(i2, 2)));
                %
                sig1000.vel5_raw{i3} = dataread_aux.Data.IBurst_VelBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
                sig1000.amp5_raw{i3} = dataread_aux.Data.IBurst_AmpBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
                sig1000.corr5_raw{i3} = dataread_aux.Data.IBurst_CorBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
    
                %
                sig1000.timednum_beam5{i3} = dataread_aux.Data.IBurst_Time(lin_proclims_beam5time_aux);
            end
            
            %
            disp(['--- Done loading data from file ' num2str(i3) ' out of ' num2str(Nfilesperseg) ' ---'])
        end

        % ------------------------------------------------
        % Concatenate cell array into a long column vector
        sig1000.timedatenum = cat(1, sig1000.timedatenum{:});
        sig1000.pressure = cat(1, sig1000.pressure{:});
        sig1000.temperature = cat(1, sig1000.temperature{:});
        %
        sig1000.heading = cat(1, sig1000.heading{:});
        sig1000.pitch = cat(1, sig1000.pitch{:});
        sig1000.roll = cat(1, sig1000.roll{:});
        % ------------------------------------------------

        % Concatenate cell arrays into matrices (loop over variables)
        for i3 = 1:length(list_beam_vars)
            sig1000.(list_beam_vars{i3}) = cat(1, sig1000.(list_beam_vars{i3}){:});
        end
        sig1000.timednum_fourbeams = cat(1, sig1000.timednum_fourbeams{:});

        % Concatenate 5th beam data
        if any(contains(list_5beams, list_Signature{i1}(1:3)))
            %
            sig1000.timednum_beam5 = cat(1, sig1000.timednum_beam5{:});
            %
            sig1000.vel5_raw = cat(1, sig1000.vel5_raw{:});
            sig1000.amp5_raw = cat(1, sig1000.amp5_raw{:});
            sig1000.corr5_raw = cat(1, sig1000.corr5_raw{:});
        end

        %
        disp(['--- Done with loading all of the data for the current time segment in: ---'])
        toc


        %% Quick pcolor plot to check the data

% %         %
% % %         indsplt_aux = 1:(8*3600);
% % %         indsplt_aux = 1:(8*60*1);    % a subset
% %         indsplt_aux = 1:length(sig1000.timednum_fourbeams);   % no subset
% % 
% %         %
% %         time_plt_aux = datetime(sig1000.timednum_fourbeams(indsplt_aux), 'ConvertFrom', 'datenum');
% % 
% %         %
% %         hfig_quickview = figure;
% % % %             %
% % % %             haxs = makeSubPlots(0.1, 0.1, 0.05, ...
% % % %                                 0.1, 0.15, 0.04, 2, 4);
% %             %
% %             haxs_1 = axes('Position', [0.1000    0.7425    0.3750    0.1575]);
% %             haxs_2 = axes('Position', [0.5250    0.7425    0.3750    0.1575]);
% %             haxs_3 = axes('Position', [0.1000    0.5450    0.3750    0.1575]);
% %             haxs_4 = axes('Position', [0.5250    0.5450    0.3750    0.1575]);
% %             haxs_5 = axes('Position', [0.1000    0.3475    0.3750    0.1575]);
% %             haxs_6 = axes('Position', [0.5250    0.3475    0.3750    0.1575]);
% %             haxs_7 = axes('Position', [0.1000    0.1500    0.3750    0.1575]);
% %             haxs_8 = axes('Position', [0.5250    0.1500    0.3750    0.1575]);
% %             %
% %             haxs = [haxs_1, haxs_2, haxs_3, haxs_4, haxs_5, haxs_6, haxs_7, haxs_8];
% %             for i3 = 1:length(haxs)
% %                 hold(haxs(i3), 'on')
% %             end
% % 
% %                 %
% %                 pcolor(haxs(1), time_plt_aux, sig1000.zhab, sig1000.vel1(indsplt_aux, :).')
% %                 pcolor(haxs(3), time_plt_aux, sig1000.zhab, sig1000.vel2(indsplt_aux, :).')
% %                 pcolor(haxs(5), time_plt_aux, sig1000.zhab, sig1000.vel3(indsplt_aux, :).')
% %                 pcolor(haxs(7), time_plt_aux, sig1000.zhab, sig1000.vel4(indsplt_aux, :).')
% %                 %
% %                 pcolor(haxs(2), time_plt_aux, sig1000.zhab, sig1000.amp1(indsplt_aux, :).')
% %                 pcolor(haxs(4), time_plt_aux, sig1000.zhab, sig1000.amp2(indsplt_aux, :).')
% %                 pcolor(haxs(6), time_plt_aux, sig1000.zhab, sig1000.amp3(indsplt_aux, :).')
% %                 pcolor(haxs(8), time_plt_aux, sig1000.zhab, sig1000.amp4(indsplt_aux, :).')
% % 
% %                 for i3 = 1:length(haxs)
% %                     plot(haxs(i3), time_plt_aux, sig1000.pressure(indsplt_aux), '-k')
% %                 end
% % 
% %         %
% %         for i3 = 1:length(haxs)
% %             shading(haxs(i3), 'flat')
% %         end
% %         
% % 
% % % %         %
% % % %         callCbrewer([], haxs(1:2:7))
% % 
% %         %
% % % %         set(haxs, 'FontSize', 16, 'Box', 'on', ...
% % % %                   'YLim', [0, (sig1000.zhab(end) + 1)])
% %         set(haxs, 'Box', 'on', 'YLim', [0, (sig1000.zhab(end) + 1)])
% %         set(haxs(1:2:7), 'CLim', 0.1.*[-1, 1])
% %         %
% %         ampclrlims = get(haxs(2), 'CLim');
% %         set(haxs(2:2:8), 'CLim', [40, (ampclrlims(2)+10)])
% %         %
% %         set(haxs, 'Color', 0.6.*[1, 1, 1])
% % 
% %         %
% %         title(haxs(1), ['Sig1000 - ' list_Signature{i1}(1:3) ' SN ' list_Signature{i1}(5:end) ' v1, v2, v3, and v4'])
% %         title(haxs(2), ['amp1, amp2, amp3, and amp4'])
% % 
% %         %
% %         set(hfig_quickview, 'units', 'normalized')
% %         set(hfig_quickview, 'Position', [0.4, 0.63, 0.38, 0.3])
% % 
% %         %
% %         linkaxes(haxs, 'xy')
% % 
% %         %
% %         time_lims_plt = get(haxs(1), 'XLim');
% %         time_xticks = linspace(time_lims_plt(1), time_lims_plt(2), 3);
% %         %
% %         set(haxs, 'XTick', time_xticks)
% %         set(haxs(1:end-2), 'XTickLabel', [])
% % 
% %         %
% % % %         exportgraphics(hfig_quickview, fullfile(pwd, ['sig_quickview_' list_Signature{i1} '_' num2str(i2, '%.3d') '.png']), 'Resolution', 300)


        %%
% %         %
% %         pause(10)
% %         %
% %         close(hfig_quickview)

% %         % DEBUGGING / QUICK PLOTS -- remove variables from previous loop
% %         sig1000 = rmfield(sig1000, list_beam_vars);

        % #################################
        % EVERYTHING BELOW WOULD GO IN HERE IF I WANTED TO DO ALL OF THE
        % PROCESSING FOR SEPARATE TIME CHUNKS!
        % #################################


    end


    %% Interpolate 5th beam to the same timestamps as the other 4

    %
    if any(contains(list_5beams, list_Signature{i1}(1:3)))
        tic
        %
        disp(['--- 5 beams will be used for velocity transformation. ' ...
              'First will interpolate 5th beam to time stamps of ' ...
              'other 4 beams ---'])

        %
        sig1000.vel5 = NaN(size(sig1000.vel1));
        sig1000.amp5 = sig1000.vel5;

        % Loop over bins of the 5th beam
        for i2 = 1:size(sig1000.vel5, 2)
            
            %
            sig1000.vel5(:, i2) = interp1(sig1000.timednum_beam5, ...
                                          sig1000.vel5_raw(:, i2), ...
                                          sig1000.timednum_fourbeams);
            %
            sig1000.amp5(:, i2) = interp1(sig1000.timednum_beam5, ...
                                          sig1000.amp5_raw(:, i2), ...
                                          sig1000.timednum_fourbeams);
        end
        disp('Done with 5th beam interpolation.') 
        toc
    end
    

    %% Transpose matrices so that row dimension is along bins

% %     %
% %     list_transpose = {'timedatenum', 'pressure', 'temperature', ...
% %                       'heading', 'pitch', 'roll', ...
% %                       'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
% %                       'amp1', 'amp2', 'amp3', 'amp4', 'amp5', ...
% %                       'corr1', 'corr2', 'corr3', 'corr4', 'corr5'};

    %
    list_transpose = {'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
                      'amp1', 'amp2', 'amp3', 'amp4', 'amp5', ...
                      'corr1', 'corr2', 'corr3', 'corr4', 'corr5'};

    %
    for i2 = 1:length(list_transpose)
        if isfield(sig1000, list_transpose{i2})
            sig1000.(list_transpose{i2}) = sig1000.(list_transpose{i2}).';
        end
    end


    %% Compute magnetic-ENU 3 components of velocity
    % (using library ADCPtools)
    %
    % Compute ENU velocities in different ways depending
    % whether 4 or 5 beams should be used (check the preamble
    % of this script. The preamble also makes sure that
    % all Signatures being processed are in the processing lists)
    %
    % Very long timeseries (e.g. a couple of weeks of 8 Hz
    % Signature1000 data (need to be broken into smaller chunks
    % otherwise janus5beam2earth will use all memory and crash Matlab).

    %
    npts_rot_TH = 1000000;
    Npts_alldata = length(sig1000.timedatenum);
    %
    if Npts_alldata<=npts_rot_TH
        indbreak_rot = [1; Npts_alldata];
    else
        %
        inds_edges_breakrot = 1 : npts_rot_TH : Npts_alldata;
        if (inds_edges_breakrot(end)~=Npts_alldata)
            inds_edges_breakrot = [inds_edges_breakrot, Npts_alldata];
        end
        %
        indbreak_rot = [inds_edges_breakrot(1:end-1); ...
                        (inds_edges_breakrot(2:end) - 1)];
        % The last one shouldn't have 1 subtracted. Add again
        indbreak_rot(2, end) = indbreak_rot(2, end) + 1;
    end


    %
    disp('--- Converting along-beam velocities to magnetic ENU ---')


    %
    disp(['--- Coordinate transformation will be computed ' ...
          'for ' num2str(size(indbreak_rot, 2)) ' separate chunks ' ...
          'of the timeseries (to avoid crashing Matlab) ---'])

    %
    sig1000.Ue = NaN(size(sig1000.vel1));
    sig1000.Vn = NaN(size(sig1000.vel1));
    sig1000.Wup = NaN(size(sig1000.vel1));
    sig1000.Wbeam5 = NaN(size(sig1000.vel1));


    tic
    % If the i1'th Signature is in the 5-beam list
    if any(contains(list_5beams, list_Signature{i1}(1:3)))
        % Compute velocity using 5 beams
        %
        disp('Using 5 beams')

        %
        for i2 = 1:size(indbreak_rot, 2)

            %
            ind_sub_aux = indbreak_rot(1, i2) : indbreak_rot(2, i2);

            %
            [sig1000.Ue(:, ind_sub_aux), ...
             sig1000.Vn(:, ind_sub_aux), ...
             sig1000.Wup(:, ind_sub_aux), ...
             sig1000.Wbeam5(:, ind_sub_aux)] = ...
                    janus5beam2earth((sig1000.heading(ind_sub_aux).' - 90), ...
                                     sig1000.roll(ind_sub_aux).', -sig1000.pitch(ind_sub_aux).', ...
                                     25, ...
                                     -sig1000.vel1(:, ind_sub_aux), -sig1000.vel3(:, ind_sub_aux), ...
                                     -sig1000.vel4(:, ind_sub_aux), -sig1000.vel2(:, ind_sub_aux), ...
                                     -sig1000.vel5(:, ind_sub_aux));
    
    
            % make sure that column dimension is time dimension for all variables
        end

    % 
    else
        % Or use 4 beams instead
        %
        disp('Using 4 beams')

        %
        for i2 = 1:size(indbreak_rot, 2)

            %
            ind_sub_aux = indbreak_rot(1, i2) : indbreak_rot(2, i2);
            %
            [sig1000.Ue(:, ind_sub_aux), ...
             sig1000.Vn(:, ind_sub_aux), ...
             sig1000.Wup(:, ind_sub_aux)] = ...
                        janus2earth(sig1000.heading(ind_sub_aux).' - 90, ...
                                    sig1000.roll(ind_sub_aux).', -sig1000.pitch(ind_sub_aux).', ...
                                    25, ...
                                    -sig1000.vel1(:, ind_sub_aux), -sig1000.vel3(:, ind_sub_aux), ...
                                    -sig1000.vel4(:, ind_sub_aux), -sig1000.vel2(:, ind_sub_aux));
        end
    end

    %
    disp('--- Done with coordinate transformation to magnetic ENU. It took: ---')
    toc


    %% Compute bottom depth from pressure

    %
    sig1000.bottomdepthfrompres = 1e4*sig1000.pressure ./ (1030*9.8);



    %% NaN data at/above the ocean surface

    %
    disp('----- Removing data close/above the ocean surface -----')

    %
    zhab_halfstep = sig1000.zhab + (sig1000.binsize/2);
    Nbins = length(sig1000.zhab);
    %
    ind_abovesurface = 1:(size(sig1000.Ue, 1) * size(sig1000.Ue, 2));
    ind_abovesurface = reshape(ind_abovesurface, size(sig1000.Ue, 1), size(sig1000.Ue, 2));
    %
    for i2 = 1:length(sig1000.timedatenum)
        %
        ind_above_aux = find((zhab_halfstep < sig1000.bottomdepthfrompres(i2)), 1, 'last');
        %
        ind_abovesurface(1:ind_above_aux, i2) = NaN;
    end
    %
    ind_abovesurface = ind_abovesurface(:);
    ind_abovesurface = ind_abovesurface(~isnan(ind_abovesurface));
    %
    sig1000.Ue(ind_abovesurface) = NaN;
    sig1000.Vn(ind_abovesurface) = NaN;
    sig1000.Wup(ind_abovesurface) = NaN;
    %
    if any(contains(list_5beams, list_Signature{i1}(1:3)))
        sig1000.Wbeam5(ind_abovesurface) = NaN;
    end
    
    

    %% Rotate horizontal velocity from magnetic to true north

    %
    disp(['----- Rotating horizontal velocity to ' ...
          sig1000.coordsystem ', with mag. declination of ' ...
          num2str(sig1000.magdec, '%.2f') ' degrees -----'])

    %
    rotMatrix = [cosd(sig1000.magdec), sind(sig1000.magdec); ...
                 -sind(sig1000.magdec), cosd(sig1000.magdec)];

% %     % Check the rotation (i.e. velocity aligned with magnetic north
% %     % should have a small zonal component and large meridional
% %     % component in a geographical north coordinate system)
% %     rotMatrix * [0; 1]

    %
    for i2 = 1:size(sig1000.Ue, 1)
        %
        uv_aux = [sig1000.Ue(i2, :); sig1000.Vn(i2, :)];
        %
        uv_rot_aux = rotMatrix * uv_aux;

        %
        sig1000.Ue(i2, :) = uv_aux(1, :);
        sig1000.Vn(i2, :) = uv_aux(2, :);
    end
    

% %     keyboard
    
    %% Filter out velocity where amplitude is below a threshold value

% %     %
% %     disp('----- Removing velocity where backscatter is smaller than threshold value -----')
% % 
% %     %
% %     aquadoppL1.backscatterTH = 30;
% % 
% %     %
% %     l_belowTH = (aquadoppL1.a1 <= aquadoppL1.backscatterTH) | ...
% %                 (aquadoppL1.a2 <= aquadoppL1.backscatterTH) | ...
% %                 (aquadoppL1.a3 <= aquadoppL1.backscatterTH);
% %     %
% %     aquadoppL1.Ue(l_belowTH) = NaN;
% %     aquadoppL1.Vn(l_belowTH) = NaN;
% %     aquadoppL1.Wup(l_belowTH) = NaN;



    %%
    % ----------------------------------------------------
    % Turn all vectors into column vectors so that Matlab
    % can quickly displace the structure variable in the 
    % command window (Matlab displays first elements of
    % row vectors, and it takes longer)
%     list_fields_aux = fieldnames(aquadoppL1);
%     %
%     for i2 = 1:length(list_fields_aux)
%         %
%         if isvector(aquadoppL1.(list_fields_aux{i2})) && ...
%            ~isstruct(aquadoppL1.(list_fields_aux{i2})) && ...
%            ~ischar(aquadoppL1.(list_fields_aux{i2}))
%             %
%             aquadoppL1.(list_fields_aux{i2}) = aquadoppL1.(list_fields_aux{i2})(:);
%         end
%     end

    %% Compute low-frequency velocity?


    %% Plot to check processed velocity data

% % %     % ---------------------------------
% % %     %
% % %     hfig_procvel = figure;
% % % 
% % %         % Create axes of subplots on regular grid
% % %         mty = 0.1;
% % %         mby = 0.1;
% % %         miy = 0.05;
% % %         ny = 4;
% % %         heighsubplot = (1 - mty - mby - (ny-1)*miy) / ny;
% % %         %
% % %         pthisy = (ny+1) - [1:ny];
% % %         %
% % %         ypos = mby + ((pthisy-1) * (heighsubplot+miy));
% % %         %
% % %         axshndls = gobjects(1, (ny*2));
% % %         
% % %         % Now create the axes
% % %         for i3 = 1:ny
% % %             axshndls((2*i3 - 1)) = axes('Position', [0.1, ypos(i3), 0.35, heighsubplot]);
% % %             axshndls(2*i3) = axes('Position', [0.55, ypos(i3), 0.35, heighsubplot]);
% % %         end
% % %         
% % % 
% % %     % ---------------------------------
% % % 
% % %     %
% % %     indsplt_aux = 1:length(sig1000.timednum_fourbeams);   % no subset
% % %     %
% % %     time_plt_aux = datetime(sig1000.timednum_fourbeams(indsplt_aux), 'ConvertFrom', 'datenum');
% % % 
% % %     
% % %     % ---------------------------------
% % % 
% % %         %
% % %         for i3 = 1:length(axshndls)
% % %             hold(axshndls(i3), 'on')
% % %         end
% % % 
% % %         %
% % %         pcolor(axshndls(1), time_plt_aux, sig1000.zhab, sig1000.Ue(:, indsplt_aux))
% % %         pcolor(axshndls(3), time_plt_aux, sig1000.zhab, sig1000.Vn(:, indsplt_aux))
% % %         pcolor(axshndls(5), time_plt_aux, sig1000.zhab, sig1000.Wup(:, indsplt_aux))
% % %         pcolor(axshndls(7), time_plt_aux, sig1000.zhab, sig1000.Wbeam5(:, indsplt_aux))
% % %         %
% % %         pcolor(axshndls(2), time_plt_aux, sig1000.zhab, sig1000.amp1(:, indsplt_aux))
% % %         pcolor(axshndls(4), time_plt_aux, sig1000.zhab, sig1000.amp2(:, indsplt_aux))
% % %         pcolor(axshndls(6), time_plt_aux, sig1000.zhab, sig1000.amp3(:, indsplt_aux))
% % %         pcolor(axshndls(8), time_plt_aux, sig1000.zhab, sig1000.amp4(:, indsplt_aux))
% % % 
% % %         for i3 = 1:length(axshndls)
% % %             plot(axshndls(i3), time_plt_aux, sig1000.bottomdepthfrompres(indsplt_aux), '-k')
% % %         end
% % % 
% % %         %
% % %         for i3 = 1:length(axshndls)
% % %             shading(axshndls(i3), 'flat')
% % %         end
% % %         
% % % 
% % % % %         %
% % % % %         callCbrewer([], haxs(1:2:7))
% % % 
% % %         %
% % % % %         set(haxs, 'FontSize', 16, 'Box', 'on', ...
% % % % %                   'YLim', [0, (sig1000.zhab(end) + 1)])
% % %         set(axshndls, 'Box', 'on', 'YLim', [0, (sig1000.zhab(end) + 1)])
% % %         set(axshndls(1:2:7), 'CLim', 0.1.*[-1, 1])
% % %         %
% % %         ampclrlims = get(axshndls(2), 'CLim');
% % %         set(axshndls(2:2:8), 'CLim', [40, (ampclrlims(2)+10)])
% % %         %
% % %         set(axshndls, 'Color', 0.6.*[1, 1, 1])
% % % 
% % %         %
% % %         title(axshndls(1), ['Sig1000 - ' list_Signature{i1}(1:3) ' SN ' list_Signature{i1}(5:end) ' v1, v2, v3, and v4'])
% % %         title(axshndls(2), ['amp1, amp2, amp3, and amp4'])
% % % 
% % %         %
% % %         set(hfig_quickview, 'units', 'normalized')
% % %         set(hfig_quickview, 'Position', [0.4, 0.63, 0.38, 0.3])
% % % 
% % %         %
% % %         linkaxes(axshndls, 'xy')
% % % 
% % %         %
% % %         time_lims_plt = get(axshndls(1), 'XLim');
% % %         time_xticks = linspace(time_lims_plt(1), time_lims_plt(2), 3);
% % %         %
% % %         set(axshndls, 'XTick', time_xticks)
% % %         set(axshndls(1:end-2), 'XTickLabel', [])


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ------------------------- FINAL STEPS --------------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------

    %% Add README

% %     %
% %     time_dataproc = datetime('now', 'TimeZone', 'Local');
% %     time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');
% %     % Add a general README
% %     sig1000.README = ['Level 1 Signature1000 data from ROXSI 2022. The data is from Signature ' ...
% %                          ' with serial number SN and deployed at mooring site mooringID. ' ...
% %                          'Data processed by script ' mfilename() '.m on ' time_dataproc_char ' (TimeZone ' time_dataproc.TimeZone '). ' ...
% %                          'Horizontal velocity components are relative to ??????, where the magnetic ' ...
% %                          'decliation was taken from www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml. Data ' ...
% %                          'in the Level 1 structure has been trimmed for the deployment ' ...
% %                          'period, as defined in the table deploymentInfo_ROXSI2022.mat. ' ...
% %                          'Pressure is in dbar, where atmospheric pressure has been removed.'];


    %%


    % ----------------------------------------------------
    % Save level 1 data and QC figures

    %
    if lsave_file
        %
        disp('----- Saving level 1 data -----')
        tic
        %
        str_filename = ['roxsi_signature_L1_' char(sig1000.mooringID) '_' char(sig1000.SN) '_velocity'];
        %
        save(fullfile(dir_output_data_L1, [str_filename '.mat']), 'sig1000', '-v7.3')

        %
        disp('----- Done saving data. It took: -----')
        toc
        
    end

    %
    if lsave_fig
% %         %
% %         disp('----- Saving level 1 QC plot figures -----')
% % 
% %         %
% %         str_filename = ['roxsi_aquadopp_L1_' list_Signature{i1}(1:3) '_' char(sig1000.SN) '_QC_1'];
% %         % Save figure as *.png
% %         exportgraphics(fig_L1_QC_tilt, fullfile(dir_output_figs_L1, [str_filename '.png']), 'Resolution', 300)
% % 
% %         %
% %         str_filename = ['roxsi_aquadopp_L1_' list_Signature{i1}(1:3) '_' char(sig1000.SN) '_QC_2'];
% %         % Save figure as *.png
% %         exportgraphics(fig_L1_QC, fullfile(dir_output_figs_L1, [str_filename '.png']), 'Resolution', 300)
    
    end


    %%

% %     % Clear some variables to avoid issues in the next loop iteration
% %     close(fig_L1_QC_tilt), close(fig_L1_QC)
% %     clear aquadoppL1 header_aux beamAQDP_aux senAQDP_aux
    %
    clear sig1000


    %%
    
    %
    disp(['----- DONE with Signature data proc: ' list_Signature{i1} ' -----'])
    toc(totalRunTime)


end

%
disp('###################### Done with data processing for all Signature1000 ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');

