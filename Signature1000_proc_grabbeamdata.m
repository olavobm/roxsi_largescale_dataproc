%% Script that grabs the beam data from Signature1000: along beam
% velocity, backscatter and correlation variables. This is done
% separately from the primary L1 processing (Signature1000_proc_lvl_1.m)
% because otherwise the required memory for processing is too much.
%
% These variables will be saved in a mat file, but not as fields
% of a structure, so that specific variables can be loaded separately.


clear
close all


%%
% --------------------------------------
% -------- SET DIRECTORY PATHS ---------
% --------------------------------------


%% Raw data directory

dirparent_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/';
dir_data_raw = fullfile(dirparent_data, 'RAW', 'Signature1000');


%% Output directory

% % dir_output_L1 = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Signature_Level1/';
dir_output_L1 = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/';


%%
% -------------------------------------------
% --- DEFINE VARIABLES FOR DATA PROCESSING --
% -------------------------------------------


%% Load ADCP deployment information

% dir_coderepo = repo_dirpath();
dir_coderepo = '/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc';
load(fullfile(dir_coderepo, 'deploymentInfo_ROXSI2022.mat'), 'deploymentInfo_ROXSI2022')


%% List of Signatures that will be processed

% % % All Signatures
% % list_Signature = {'A01_103043', ...
% %                   'B10_103045', ...
% %                   'B13_103046', ...
% %                   'B15_103056', ...
% %                   'B17_101923', ...
% %                   'C01_102128', ...
% %                   'X11_101941', ...
% %                   'X05_100231'};
% % 
% % % Just a test
list_Signature = {'A01_103043'};
% % list_Signature = {'B10_103045'};
% % list_Signature = {'B13_103046'};

%
Nsignatures = length(list_Signature);


%% The list of ADCPs processed with either 4 or 5 beams

%
list_4beams = {'B10'};    % B10 has 5th beam, but on HR mode (up to 4 meters above the transducer). Others did the echosounder.
list_5beams = {'A01', 'B13', 'B15', 'B17', 'C01', 'X05', 'X11'};

% Check all of the Signatures in the list are in either
% list for Signatures with 4 or 5 beams
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


%% List of all variables/field names in the data structure
% as data is first stored in data structure

% % %
% % list_rawdata = {'timedatenum', 'pressure', 'temperature', ...
% %                 'heading', 'pitch', 'roll', ...
% %                 'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
% %                 'amp1', 'amp2', 'amp3', 'amp4', 'amp5', ...
% %                 'cor1', 'cor2', 'cor3', 'cor4', 'cor5', ...
% %                 'timedatenum5'};

% Remove variables to free up memory
list_rawdata = {'timedatenum', 'pressure', 'temperature', ...
                'heading', 'pitch', 'roll', ...
                'vel1', 'vel2', 'vel3', 'vel4', 'vel5', ...
                'timedatenum5'};



%%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% -------------------------- DO DATA PROCESSING --------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% Initialize a log file with what is printed to the
% command window and timer for running the whole script

%
log_file_name = ['log_Signature_procL1_at_' datestr(datetime('now', 'TimeZone', 'Local'), 'yyyymmdd_HHMMSS') '.txt'];
%
diary(fullfile(dir_output_L1, log_file_name))

%
totalRunTime = tic;


%% Display on the screen:

%
disp(' '), disp(' ')
disp('------------------------------ Processing Signature1000 data from RAW to L1 ------------------------------')
disp('List of Signature1000s being processed:')
%
for i = 1:Nsignatures
    disp([num2str(i) ' - ' list_Signature{i}])
end


%% Loop over Signatures in the list and process dat

% Loop over Signature1000's in the list
for i1 = 1:Nsignatures


    %%
    % ------------------------------------------
    % ----- PREAMBLE OF THE DATA PROCESSING ----
    % ------------------------------------------

    %% Display progress messahe on the screen
    %
    disp(' '), disp(' ')
    disp(['--- Start processing Signature1000 data from: ' list_Signature{i1} ' ---'])


    %% Get the Signature1000 files in the correct order

    %
    dir_data_aux = fullfile(dir_data_raw, list_Signature{i1}, 'converted');
    %
    list_dir_aux = Signature_orderfiles(dir_data_aux, list_Signature{i1}(1:3));

    %
    Nfiles = length(list_dir_aux);


    %% First add basic metadata to structure variable

    %
    sigL1.SN = convertCharsToStrings(list_Signature{i1}(5:end));
    sigL1.mooringID = convertCharsToStrings(list_Signature{i1}(1:3));

    % Latitude/longitude
    info_mooringtable = ROXSI_mooringlocation(sigL1.mooringID, "ADCP");
    %
    sigL1.latitude = info_mooringtable.latitude;
    sigL1.longitude = info_mooringtable.longitude;

    %
    sigL1.site = info_mooringtable.roxsiarray;

    %
    [sigL1.X, sigL1.Y] = ROXSI_lltoxy(sigL1.latitude, sigL1.longitude, sigL1.site);

    % In clockwise degrees from the true north
    sigL1.magdec = 12.86;


    %%

    %
    if any(contains(list_5beams, list_Signature{i1}(1:3)))
        sigL1.l5beams = true;
    else
        sigL1.l5beams = false;
    end

    
    %% Get trimming times for the deployment period of this Signature

    %
    ind_row_match = find(strcmp(deploymentInfo_ROXSI2022.SN, list_Signature{i1}(5:end)));
    
    %
    time_1 = deploymentInfo_ROXSI2022.time_begin_trim(ind_row_match);
    time_2 = deploymentInfo_ROXSI2022.time_end_trim(ind_row_match);

    %
    time_1 = datetime(datenum(time_1, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
    time_2 = datetime(datenum(time_2, 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');

    %
    time_lims_proc = [time_1, time_2];
    Ndatasegments = size(time_lims_proc, 1);
    

    %%
    % ------------------------------------------
    % ------ FIRST LOAD FILE WITH SCALARS ------
    % ------------------------------------------


    %%

    %
    sigL1scalars = load(fullfile(dir_output_L1, ['roxsi_signature_L1_' char(sigL1.mooringID) '_' char(sigL1.SN) '_scalars.mat']));
    field_aux = fieldnames(sigL1scalars);
    sigL1scalars = sigL1scalars.(field_aux{1});

    %% Copy some variables

    %
    sigL1metadata.SN = sigL1scalars.SN;
    sigL1metadata.binsize = sigL1scalars.binsize;
    sigL1metadata.zhab = sigL1scalars.zhab;
    

    %% Find the last bin that is not outside of the ocean
    % at all times (this must match Signature1000_proc_lvl_1.m)


    % Find all bins that are below the maximum bottom depth
    lin_verticalrange = ((sigL1metadata.zhab + (sigL1metadata.binsize)) < ...
                         max(sigL1scalars.bottomdepthfrompres));

    keyboard
    %%
    % ------------------------------------------
    % ------ LOAD CONFIGURATION INFO FROM ------
    % ------------- ONE DATA FILE --------------
    % ------------------------------------------

% %     %%
% % 
% %     % Just the Config structure variable
% %     dataread_aux = load(fullfile(dir_data_aux, list_dir_aux(1).name), 'Config');
% % 
% %     %
% %     sigL1.Config = dataread_aux.Config;
% % 
% % 
% %     %% Get height above the bottom (zhab) of the ADCP bins
% % 
% %     % In meters (based on the Solidworks drawing)
% %     sigL1.transducerHAB = (31.88)/100;
% % 
% %     % In meters
% %     sigL1.binsize = dataread_aux.Config.Burst_CellSize;
% %     
% %     % Height of the first cell center relative to transducer
% %     % (based on the Principles of Operation manual by Nortek, page 12)
% %     cellcenter_first_bin = dataread_aux.Config.Burst_BlankingDistance + ...
% %                            dataread_aux.Config.Burst_CellSize;
% % 
% %     %
% %     sigL1.cellcenter = cellcenter_first_bin + ...
% %                         (0:1:(double(dataread_aux.Config.Burst_NCells) - 1)) .* sigL1.binsize;
% % 
% %     %
% %     sigL1.zhab = sigL1.transducerHAB + sigL1.cellcenter;


    %%
    % ------------------------------------------
    % ------------- LOAD THE DATA --------------
    % ------------------------------------------


    %%

    %
    for i2 = 1:Ndatasegments

        %%

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

% %         %
% %         disp(['--- Raw data that will processed is split into ' num2str(Nfilesperseg) ' data files. Loading data from ---'])
% %         listfiles_perseg

        %% Pre-allocate space for variables that will be read

        % Pre-allocate in cell arrays
        prealloc_aux = cell(Nfilesperseg, 1);

        %
        sigL1.timedatenum = prealloc_aux;
        sigL1.samplingrateHz = NaN;    % only added at the end
        %
        sigL1.dtime = [];    % variable not filled in the loop. I just
                             % want to have its position here
% %         sigL1.pressure = prealloc_aux;

        %
        sigL1.vel1 = prealloc_aux;
        sigL1.vel2 = prealloc_aux;
        sigL1.vel3 = prealloc_aux;
        sigL1.vel4 = prealloc_aux;
        %
        sigL1.amp1 = prealloc_aux;
        sigL1.amp2 = prealloc_aux;
        sigL1.amp3 = prealloc_aux;
        sigL1.amp4 = prealloc_aux;
        %
        sigL1.cor1 = prealloc_aux;
        sigL1.cor2 = prealloc_aux;
        sigL1.cor3 = prealloc_aux;
        sigL1.cor4 = prealloc_aux;

        %
        if sigL1.l5beams
            %
            sigL1.timedatenum5 = prealloc_aux;
            sigL1.dtime5 = prealloc_aux;
            %
            sigL1.vel5 = prealloc_aux;
            sigL1.amp5 = prealloc_aux;
            sigL1.cor5 = prealloc_aux;
        end
        

        %% Finally load the raw data

        %
        tic
        % Load over all files with data in the i2'th segment
        for i3 = 1:Nfilesperseg

            %
            dataread_aux = load(fullfile(dir_data_aux, ...
                                         listfiles_perseg(i3)));

            %
            lin_proclims_beam4time_aux = (dataread_aux.Data.Burst_Time >= datenum(time_lims_proc(i2, 1))) & ...
                                         (dataread_aux.Data.Burst_Time <  datenum(time_lims_proc(i2, 2)));

            % -------------------------------
            % First get all the scalars
            %
            sigL1.timedatenum{i3} = dataread_aux.Data.Burst_Time(lin_proclims_beam4time_aux);
% %             sigL1.pressure{i3} = dataread_aux.Data.Burst_Pressure(lin_proclims_beam4time_aux);
            

            % -------------------------------
            % Then get velocity data
            
            % Dummy/for code development
            lin_verticalrange = true(1, size(dataread_aux.Data.Burst_VelBeam1, 2));

            % Get data from the 4 beams
            sigL1.vel1{i3} = dataread_aux.Data.Burst_VelBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel2{i3} = dataread_aux.Data.Burst_VelBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel3{i3} = dataread_aux.Data.Burst_VelBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.vel4{i3} = dataread_aux.Data.Burst_VelBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            %
            sigL1.amp1{i3} = dataread_aux.Data.Burst_AmpBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.amp2{i3} = dataread_aux.Data.Burst_AmpBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.amp3{i3} = dataread_aux.Data.Burst_AmpBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.amp4{i3} = dataread_aux.Data.Burst_AmpBeam4(lin_proclims_beam4time_aux, lin_verticalrange);
            %
            sigL1.cor1{i3} = dataread_aux.Data.Burst_CorBeam1(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.cor2{i3} = dataread_aux.Data.Burst_CorBeam2(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.cor3{i3} = dataread_aux.Data.Burst_CorBeam3(lin_proclims_beam4time_aux, lin_verticalrange);
            sigL1.cor4{i3} = dataread_aux.Data.Burst_CorBeam4(lin_proclims_beam4time_aux, lin_verticalrange);

            %
            if sigL1.l5beams
                %
                lin_proclims_beam5time_aux = (dataread_aux.Data.IBurst_Time >= datenum(time_lims_proc(i2, 1))) & ...
                                             (dataread_aux.Data.IBurst_Time <  datenum(time_lims_proc(i2, 2)));
                %
                sigL1.vel5{i3} = dataread_aux.Data.IBurst_VelBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
                sigL1.amp5{i3} = dataread_aux.Data.IBurst_AmpBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
                sigL1.cor5{i3} = dataread_aux.Data.IBurst_CorBeam5(lin_proclims_beam5time_aux, lin_verticalrange);
    
                %
                sigL1.timedatenum5{i3} = dataread_aux.Data.IBurst_Time(lin_proclims_beam5time_aux);
            end
            
            %
            disp(['--- Done loading data from ' num2str(i3) ' out of ' num2str(Nfilesperseg) ' files ---'])

            % Clear to free up some memory
            clear dataread_aux
        end


        %% Concatenate data in cell arrays into matrices

        %
        for i3 = 1:length(list_rawdata)
            %
            if isfield(sigL1, list_rawdata{i3})
                %
                sigL1.(list_rawdata{i3}) = cat(1, sigL1.(list_rawdata{i3}){:});

                % Turn integers into single precision variables
                % so they can be interpolated (correlations are
                % given as integers)
                if isinteger(sigL1.(list_rawdata{i3}))    
                    sigL1.(list_rawdata{i3}) = single(sigL1.(list_rawdata{i3}));
                end
            end
        end

        %
        disp('--- Done with loading all of the data for the current time segment in: ---')
        toc

    end

    keyboard

    %%
    % ------------------------------------------
    % ---------- BASIC DATA PROCESSING ---------
    % ------------------------------------------

    %% Clock drift correction

    %
    disp('----- Correcting for clock drift -----')

    % Now correct for clock drift
    time_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, ...
                                            list_Signature{i1}(5:end), ...
                                            sigL1.timedatenum);

    % Now correct for clock drift
    if sigL1.l5beams
        time5_aux = ROXSI_rescaletime_instrument(deploymentInfo_ROXSI2022, ...
                                                 list_Signature{i1}(5:end), ...
                                                 sigL1.timedatenum5);
    end


    %% Convert time to date time

    %
    sigL1.dtime = datetime(time_aux, 'ConvertFrom', 'datenum', ...
                                     'TimeZone', 'America/Los_Angeles');

    %
    if sigL1.l5beams
        %
        sigL1.dtime5 = datetime(time5_aux, 'ConvertFrom', 'datenum', ...
                                           'TimeZone', 'America/Los_Angeles');
    end




    %% Check clock/gaps

    % Check for something that should never happen at this point
    if any(isnan(sig1000.timedatenum))
        warning(['###### Signature ' list_Signature{i1} ' has ' ...
                 'invalid (NaN) timestamps ######'])
    end

    % Now do a plot for cheching diff(time)
    disp('--- QC plot checking diff time ---')

    %
    inds_time = 1:length(sig1000.dtime);
    inds_difftime = (inds_time(1:end-1) + inds_time(2:end))./2;

    %
    fig_L1_QC_clock = figure;
    set(fig_L1_QC_clock, 'units', 'normalized')
    set(fig_L1_QC_clock, 'Position', [0.2, 0.2, 0.4, 0.6])
        %
        haxs_1 = axes(fig_L1_QC_clock, 'Position', [0.2, 0.575, 0.6 , 0.325]);
        haxs_2 = axes(fig_L1_QC_clock, 'Position', [0.2, 0.150, 0.6, 0.325]);
        %
        haxs_all = [haxs_1, haxs_2];
        hold(haxs_all, 'on')
        %
        plot(haxs_1, inds_time, sig1000.dtime, '-k')
        plot(haxs_2, inds_difftime, seconds(diff(sig1000.dtime)), '-k')

    %
    set(haxs_all, 'FontSize', 12, 'Box', 'on', ...
                  'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_all, 'XLim', [0, (inds_time(end) + 1)])
    %
    ylim(haxs_1, sig1000.dtime([1, end]))

    %
    xlabel(haxs_2, 'Indices', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    ylabel(haxs_1, 'Time', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, 'seconds', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, ['ROXSI 2022: Signature  ' char(sig1000.mooringID) ' - SN ' ...
                   char(sig1000.SN) ': time and diff(time) (in seconds)'], ...
                  'Interpreter', 'Latex', 'FontSize', 16)
    %
    linkaxes(haxs_all, 'x')


    % Plot horizontal lines for trimming edges
    xlims_aux = xlim(haxs_all(1));
    %
    plot(haxs_all(1), xlims_aux, [time_1, time_1], '--r')
    plot(haxs_all(1), xlims_aux, [time_2, time_2], '--r')
    %
    xlim(haxs_all(1), xlims_aux)
    
    
    %% Interpolate 5th beam to the same timestamps as the other 4

    %
    if sigL1.l5beams

        tic
        %
        disp(['--- 5 beams will be used for velocity transformation. ' ...
              'Interpolating 5th along-beam velocity (and backscatter ' ...
              'and correlation) to time stamps of other 4 beams ---'])

        %
        vel5_aux = NaN(size(sigL1.vel1));
% %         amp5_aux = vel5_aux;
% %         cor5_aux = vel5_aux;

        % Loop over bins of the 5th beam
        for i2 = 1:size(sigL1.vel5, 2)
            
            %
            vel5_aux(:, i2) = interp1(sigL1.timedatenum5, ...
                                      sigL1.vel5(:, i2), ...
                                      sigL1.timedatenum);
% %             %
% %             amp5_aux(:, i2) = interp1(sigL1.timedatenum5, ...
% %                                       sigL1.amp5(:, i2), ...
% %                                       sigL1.timedatenum);
% %             %
% %             cor5_aux(:, i2) = interp1(sigL1.timedatenum5, ...
% %                                       single(sigL1.cor5(:, i2)), ...
% %                                       sigL1.timedatenum);
        end

        % Replace
        sigL1.vel5 = vel5_aux;
% %         sigL1.amp5 = amp5_aux;
% %         sigL1.cor5 = cor5_aux;
        
        %
        disp('Done with 5th beam interpolation.') 
        toc
    end


    %% Interpolate variables to gridded time (after
    % making sure there are no major issues above)
    %
    % At this point, time is always along the row dimension

    tic
    %
    disp('--- Gridding variables to time grid ---')

    % As time edges of the grid, round to the whole minute (done in a more
    % complicated way so that it doesn't rely on Matlab versions newer than
    % at least 2021b)
    dtime_edge_1 = datetime(sigL1.dtime(1).Year, sigL1.dtime(1).Month, sigL1.dtime(1).Day, ...
                            sigL1.dtime(1).Hour, (sigL1.dtime(1).Minute + 1), 00);
    dtime_edge_2 = datetime(sigL1.dtime(end).Year, sigL1.dtime(end).Month, sigL1.dtime(end).Day, ...
                            sigL1.dtime(end).Hour, sigL1.dtime(end).Minute, 00);
    %
    df_sampling = double(sigL1.Config.Burst_SamplingRate);    % in Hertz

    %
    dtime_grid = dtime_edge_1 : seconds(1/df_sampling) : dtime_edge_2;
    dtime_grid.TimeZone = sigL1.dtime.TimeZone;
    dtime_grid = dtime_grid(:);

    %
    Nlengthtimeseries = length(sigL1.dtime);
    %
    list_time_vars = {'dtime', 'dtime5', 'timedatenum', 'timedatenum5'};
    list_all_fields = fieldnames(sigL1);

    %
    for i2 = 1:length(list_all_fields)

        % Only interpolate fields that are NOT contained in list_time_vars
        % and are arrays with expected number of points along time
        if ~any(contains(list_time_vars, list_all_fields{i2})) && ...
           (size(sigL1.(list_all_fields{i2}), 1) == Nlengthtimeseries)
           
            %
            sigL1.(list_all_fields{i2}) = ...
                            interp1(sigL1.dtime, ...              
                                    sigL1.(list_all_fields{i2}), ...
                                    dtime_grid);
            % PS: matrices are interpolated column-wise

% %             % Turn to column vector
% %             sigL1.(list_variables_aux{i2}) = sigL1.(list_variables_aux{i2})(:);

        end
    end

    % Replace measured time stamps by time grid
    sigL1.dtime = dtime_grid;

    %
    disp('--- Done with time gridding ---')
    toc

% %     % Remove datenum time
% %     sigL1 = rmfield(sigL1, 'timedatenum');

    % Add sampling rate as a field
    sigL1.samplingrateHz = df_sampling;


end


%%

%
disp('###################### Done with processing from RAW to L1 for all Signature1000 ######################')

%
disp(' '), disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)

% Close the log file
diary('off');

