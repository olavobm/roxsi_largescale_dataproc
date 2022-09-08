%% Define variable with time limits to remove bad Smart Mooring data.
% There are a few different reasons of this bad data.
%
%
% Spotters that broke and and had intermmitently bad data are:
%       E07_spot1855, E09_spot1850, E13_spot1849
%
% Spotters that have clock going backwards but not coming back as expected
% (this whole thing happening several minutes after the mode change gap):
%       E05_spot1853, E08_spot1852, E09_spot1850
%
% Spotter with very weird wrong pressure data that lasts for 30s around
% 12:31pm on deployment day:
%       E10_spot1848
%
% E09_spot1850 is so bad after cable broke, that I do not keep any data
% after cable broke and instead I sent the time_end_trim to remove all
% of that.


clear
close all


%%

%
dir_rawdata_parent = '/Volumes/LaCie/RAW/Spotters_Smart/SDcards';


%%

%
dplySpotters = load('/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/code_proc/deploymentInfo_Spotters_ROXSI2022');
dplySpotters = dplySpotters.dployInfo_Spotters;



%%
% ------------------------------------------------
% ---- FOR FAULTY SMART MOORINGS, IDENTIFY BAD ---
% ------ DATA TIMES WHEN CONTINUOUS DATA ONLY ----
% ---------- HAPPENS IN A SHORT SEGMENT ----------
% ------------------------------------------------


%% Read E07 1855 data

%
raw_readdata = Spotter_readmulti_SMD(fullfile(dir_rawdata_parent, 'E07_spot1855'));

% In datenum, convert from UTC to local time (PDT)
raw_readdata.allfiles.dtime = raw_readdata.allfiles.dtime - (7/24);

%
l_gooddata = (raw_readdata.allfiles.link==1);
%
faultySmartMoorings.E07_spot1855.dtime = raw_readdata.allfiles.dtime(l_gooddata);
faultySmartMoorings.E07_spot1855.pressure = raw_readdata.allfiles.pressure(l_gooddata);



%% Read E13 1849 data

%
raw_readdata = Spotter_readmulti_SMD(fullfile(dir_rawdata_parent, 'E13_spot1849'));

% In datenum, convert from UTC to local time (PDT)
raw_readdata.allfiles.dtime = raw_readdata.allfiles.dtime - (7/24);

%
l_gooddata = (raw_readdata.allfiles.link==1);
%
faultySmartMoorings.E13_spot1849.dtime = raw_readdata.allfiles.dtime(l_gooddata);
faultySmartMoorings.E13_spot1849.pressure = raw_readdata.allfiles.pressure(l_gooddata);


%%

%
list_fields = fieldnames(faultySmartMoorings);

% *****
% AT LEAST IN 2022, THE SMART MOORINGS
% HAD A BUG/FEATURE WHERE ALL PRESSURE VALUES
% ARE RECORDED WITH AN EXTRA ZERO AS THE LAST
% DIGIT. THAT IS, THE PRESSURE IN PASCAL IS
% THE VALUE STORED IN THE SD CARD DIVIDED BY 10.
% *****
%
for i = 1:length(list_fields)
    %
    faultySmartMoorings.(list_fields{i}).pressure = faultySmartMoorings.(list_fields{i}).pressure./10;
end


%% Trim according to deployment table


%
for i = 1:length(list_fields)


    % Find the Smart Mooring in the deployment table
    lmatch = strcmp(dplySpotters.SN, list_fields{i}(end-3:end));

    %
    time_lims_aux = [datenum(dplySpotters(lmatch, :).time_begin_trim, "yyyy/mm/dd HH:MM:SS"), ...
                     datenum(dplySpotters(lmatch, :).time_end_trim, "yyyy/mm/dd HH:MM:SS")];

    %
    lintime_lims_aux = (faultySmartMoorings.(list_fields{i}).dtime >= time_lims_aux(1)) & ...
                       (faultySmartMoorings.(list_fields{i}).dtime <= time_lims_aux(2));

    %
    faultySmartMoorings.(list_fields{i}).dtime = faultySmartMoorings.(list_fields{i}).dtime(lintime_lims_aux);
    faultySmartMoorings.(list_fields{i}).pressure = faultySmartMoorings.(list_fields{i}).pressure(lintime_lims_aux);
    
        
end
        


%% Define a gap threshold (i.e. define what a "continuous" timeseries is)

% In seconds
gap_length_TH = 5;

% Continuous segment requirement/threshold (in minutes)
minlength_segment_TH = 120;


%% Remove clearly wrong spikes

%
for i = 1:length(list_fields)

    %
    if strcmp(list_fields{i}(end-3:end), '1849')    % E13

        %
        lims_good_data = [1.9e5, 2.5e5];

    elseif strcmp(list_fields{i}(end-3:end), '1855')    % E07

        %
        lims_good_data = [1.7e5, 2.3e5];
    end

    %
    lingoodlims_aux = (faultySmartMoorings.(list_fields{i}).pressure >= lims_good_data(1)) & ...
                      (faultySmartMoorings.(list_fields{i}).pressure <= lims_good_data(2));
    %
    faultySmartMoorings.(list_fields{i}).dtime = faultySmartMoorings.(list_fields{i}).dtime(lingoodlims_aux);
    faultySmartMoorings.(list_fields{i}).pressure = faultySmartMoorings.(list_fields{i}).pressure(lingoodlims_aux);



end

    

%% Loop over Smart moorings and find times bounding
% segments that do not have long enough continuous data

%
for i = 1:length(list_fields)
    

    %
    ind_gaps_aux = find(24*3600*diff(faultySmartMoorings.(list_fields{i}).dtime) > gap_length_TH);

    % Find bounding times of continuous segments (notice that if there
    % is an isolated point, the numbers in the first and second column
    % will be the same in a row -- that happens for most rows in the
    % faulty smart moorings)
    segments_timedges = [faultySmartMoorings.(list_fields{i}).dtime(ind_gaps_aux(1:end-1)+1), ...
                         faultySmartMoorings.(list_fields{i}).dtime(ind_gaps_aux(2:end))];

    % Get which of the continuous segments are long enough
    llong_enough_aux = (24*60*(segments_timedges(:, 2) - segments_timedges(:, 1))) >= minlength_segment_TH;

    % Select the time edges of good/long continuous data
    time_trim_data = segments_timedges(llong_enough_aux, :);
    % Edit the edges of the time_trim_data
    time_trim_data(1, 1) = faultySmartMoorings.(list_fields{i}).dtime(1);
    %
    l_keepdata = false(1, length(faultySmartMoorings.(list_fields{i}).dtime));

    %
    for i2 = 1:size(time_trim_data, 1)
        %
        lin_segment = (faultySmartMoorings.(list_fields{i}).dtime >= time_trim_data(i2, 1)) & ...
                      (faultySmartMoorings.(list_fields{i}).dtime <= time_trim_data(i2, 2));
        %
        l_keepdata(lin_segment) = true;
    end

    %
    faultySmartMoorings.(list_fields{i}).dtime = faultySmartMoorings.(list_fields{i}).dtime(l_keepdata);
    faultySmartMoorings.(list_fields{i}).pressure = faultySmartMoorings.(list_fields{i}).pressure(l_keepdata);


    % Plot to check the result
    figure
        %
        set(gcf, 'Units', 'normalized')
        set(gcf, 'Position', [0.3273, 0.4257, 0.4445, 0.2576])
        %
        plot(datetime(faultySmartMoorings.(list_fields{i}).dtime, 'ConvertFrom', 'datenum'), ...
             faultySmartMoorings.(list_fields{i}).pressure, '.-k')
        hold on

        %
        set(gca, 'FontSize', 16, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
        ylabel(gca, 'Total pressure [Pa]', 'Interpreter', 'Latex', 'FontSize', 20)
        %
        title(['ROXSI 2022 - keeping only long segments of good data ' ...
               list_fields{i}(1:3) ' - SN ' list_fields{i}(end-3:end)], ...
               'Interpreter', 'Latex', 'FontSize', 22)

        % Plot beginning of time edges in blue
        for i2 = 1:size(time_trim_data, 1)
            %
            time_aux = datetime(time_trim_data(i2, 1), 'ConvertFrom', 'datenum');
            time_plt = [time_aux, time_aux];
            %
            plot(time_plt, ylim, '-b', 'LineWidth', 2)
        end
        % Plot end of time edges in red
        for i2 = 1:size(time_trim_data, 1)
            %
            time_aux = datetime(time_trim_data(i2, 2), 'ConvertFrom', 'datenum');
            time_plt = [time_aux, time_aux];
            %
            plot(time_plt, ylim, '-r', 'LineWidth', 2)
        end

    % ----------------------------------------
    % Store time_trim_data
    removeBadData.(list_fields{i}).time_trim_data = time_trim_data;

end

%%
% ------------------------------------------------
% ------- MANNUALLY ADD TRIM TIME EDGES FOR ------
% ------------- OTHER SMART MOORINGS -------------
% ------------------------------------------------


%% E02 1859 clock reversal that doesn't reset in a consistent way
% (after resetting Spotter mode on dashboard)

removeBadData.E02_spot1859.time_trim_data = [datenum(2022, 06, 17, 07, 25, 00), datenum(2022, 06, 17, 17, 02, 20); ...
                                             datenum(2022, 06, 17, 17, 08, 15), datenum(2022, 07, 20, 06, 05, 00)];


%% E05 1853 clock reversal that doesn't reset in a consistent way
% (after resetting Spotter mode on dashboard)

removeBadData.E05_spot1853.time_trim_data = [datenum(2022, 06, 17, 07, 35, 00), datenum(2022, 06, 17, 17, 15, 00); ...
                                             datenum(2022, 06, 17, 17, 20, 30), datenum(2022, 07, 17, 10, 13, 00)];


%% E07 1857 clock reversal that doesn't reset in a consistent way
% (after resetting Spotter mode on dashboard)

removeBadData.E07_spot1857.time_trim_data = [datenum(2022, 07, 05, 10, 45, 00), datenum(2022, 07, 05, 11, 20, 14); ...
                                             datenum(2022, 07, 05, 11, 24, 27), datenum(2022, 07, 20, 06, 12, 00)];



%% E08 1852 clock reversal that doesn't reset in a consistent way
% (after resetting Spotter mode on dashboard)

removeBadData.E08_spot1852.time_trim_data = [datenum(2022, 06, 17, 08, 40, 00), datenum(2022, 06, 17, 17, 25, 45); ...
                                             datenum(2022, 06, 17, 17, 36, 40), datenum(2022, 07, 20, 06, 42, 00)];



%% E09 1850 clock reversal that doesn't reset in a consistent way
% (after resetting Spotter mode on dashboard)

removeBadData.E09_spot1850.time_trim_data = [datenum(2022, 06, 17, 08, 08, 00), datenum(2022, 06, 17, 17, 22, 40); ...
                                             datenum(2022, 06, 17, 17, 41, 30), datenum(2022, 06, 21, 14, 00, 00)];


%% E10 1848 has weird short bad data
% (after resetting Spotter mode on dashboard)

removeBadData.E10_spot1848.time_trim_data = [datenum(2022, 06, 17, 08, 20, 00), datenum(2022, 06, 17, 12, 31, 15); ...
                                             datenum(2022, 06, 17, 12, 31, 48), datenum(2022, 07, 20, 06, 56, 00)];



%% Add gap associated with resetting mode for E13 1849 -- this one 
% doesn't seem to have bad clock reversal after reset gap, but
% it does have a seemingly double-reversal before the reset gap.
% Set time limits to deal with that

% It would be like this if there was not other gap
% % removeBadData.E13_spot1849.time_trim_data = [datenum(2022, 06, 17, 08, 32, 00), datenum(2022, 06, 17, 17, 04, 30); ...
% %                                              datenum(2022, 06, 17, 17, 07, 20), datenum(2022, 07, 20, 06, 27, 00)];

% It would be like this if there was not other gap -- extend gap to deal
% with clock issue
% % removeBadData.E13_spot1849.time_trim_data = [datenum(2022, 06, 17, 08, 32, 00), datenum(2022, 06, 17, 17, 04, 30); ...
% %                                              datenum(2022, 06, 17, 17, 18, 50), datenum(2022, 07, 20, 06, 27, 00)];


% Appending with what has already been created:
temporary_var = [removeBadData.E13_spot1849.time_trim_data(1, :); ...
                 NaN(1, 2); ...
                 removeBadData.E13_spot1849.time_trim_data(2:end, :)];

% Replace NaN with existing time
temporary_var(2, 2) = temporary_var(1, 2);
% Insert the new time edges
temporary_var(1, 2) = datenum(2022, 06, 17, 17, 04, 30);
temporary_var(2, 1) = datenum(2022, 06, 17, 17, 18, 50);

%
removeBadData.E13_spot1849.time_trim_data = temporary_var;


%% Add gap associated with resetting mode for E07 1855 -- more uncertain
% where clock reversal problem is, but maybe it could be a 2-step
% reset that takes longer than most cases


% % %
% % removeBadData.E07_spot1855.time_trim_data = [datenum(2022, 06, 17, 08, 40, 00), datenum(2022, 06, 17, 17, 12, 05); ...
% %                                              datenum(2022, 06, 17, 17, 20, 30), datenum(2022, 07, 05, 10, 00, 00)];


% Appending with what has already been created:
temporary_var = [removeBadData.E07_spot1855.time_trim_data(1, :); ...
                 NaN(1, 2); ...
                 removeBadData.E07_spot1855.time_trim_data(2:end, :)];

% Replace NaN with existing time
temporary_var(2, 2) = temporary_var(1, 2);
% Insert the new time edges
temporary_var(1, 2) = datenum(2022, 06, 17, 17, 12, 05);
temporary_var(2, 1) = datenum(2022, 06, 17, 17, 20, 30);

%
removeBadData.E07_spot1855.time_trim_data = temporary_var;



%%
% ------------------------------------------------
% - ADD SOME MORE INFORMATION AND SAVE STRUCTURE -
% ------------------------------------------------

% Sort fields in increasing number of mooring ID
removeBadData = orderfields(removeBadData);

%
removeBadData.timezone = 'PDT';

%
removeBadData.README = ['Time limits (in datenum) defining times where ' ...
                        'is at least OK. Data OUTSIDE these ranges ' ...
                        'should be removed. The are multiple reasons ' ...
                        'for the source of problems outside of these ranges'];

%% Save variable

%
save(fullfile(repo_dirpath(), 'smart_mooring_extratrim_gooddata.mat'), 'removeBadData');

