%% Merge L1 Smart Mooring data from the Smart Moorings (E07 and E09)
% that were replaced in the middle of the experiment.
%
% The gridded data is merged.
%
% The gridded data has the same time grid (at least for each mooring)
% so I'll use that in this code.


clear
close all

%%
% -------------------------------------------
% ------------- DEFINE VARIABLES ------------
% -------------------------------------------

%%

%
% dir_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Smart_Level1/gridded/';
dir_data = ['/home/omarques/Documents/obm_ROXSI/obm_DataLocal' ...
            '/Level1_Data/Spotter_Smart_Level1/gridded/'];

%
moordID_merge = {'E07', 'E09'};

%
name_string_beginning = 'roxsi_smartmooring_L1_';

%
% % dir_output = dir_data;
dir_output = ['/home/omarques/Documents/obm_ROXSI/obm_DataLocal' ...
              '/Level1_Data/Spotter_Smart_Level1/merged/'];


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ---------------- .................. ----------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------



for i1 = 1:length(moordID_merge)

    %% LOAD DATA 
    
    %
    list_dirSpotter = dir(fullfile(dir_data, [name_string_beginning, moordID_merge{i1} '*.mat']));

    %
    Nfiles = length(list_dirSpotter);

    %
    if Nfiles<2
        error(['There are NOT at least 2 Smart Mooring data files ' ...
               'at mooring ID ' moordID_merge{i1} ' to be merged.'])
    end
    
    
    %% SORT THE ORDER OF THE FILES

    %
    time_edges = NaT(Nfiles, 2);

    %
    for i2 = 1:Nfiles
        %
        dataAll(i2) = load(fullfile(list_dirSpotter(i2).folder, list_dirSpotter(i2).name));

        %
        if i2==1
            time_edges.TimeZone = dataAll(i2).spotsmartL1.dtime.TimeZone;
        end

        %
        ind_edges_aux = [find(~isnan(dataAll(i2).spotsmartL1.pressure), 1, 'first'), ...
                         find(~isnan(dataAll(i2).spotsmartL1.pressure), 1, 'last')];
        %
        time_edges(i2, :) = dataAll(i2).spotsmartL1.dtime(ind_edges_aux);

    end

    %
    [~, ind_sort] = sort(time_edges(:, 1)); 

    
    %% Start with metadata

    %
    spotsmartL1.mooringID = moordID_merge{i1};

    %
    spotsmartL1.SN = [];
    for i2 = 1:Nfiles
        %
        if i2~=Nfiles
            spotsmartL1.SN = [spotsmartL1.SN, dataAll(ind_sort(i2)).spotsmartL1.SN ' + '];
        else 
            spotsmartL1.SN = [spotsmartL1.SN, dataAll(ind_sort(i2)).spotsmartL1.SN];
        end
    end

    %
    % spotsmartL1.site = dataAll(1).spotterL1.site;
    spotsmartL1.site = "ChinaRock";
    spotsmartL1.latitude = dataAll(1).spotsmartL1.latitude;
    spotsmartL1.longitude = dataAll(1).spotsmartL1.longitude;
    % % spotsmartL1.X = dataAll(1).spotterL1.X;
    % % spotsmartL1.Y = dataAll(1).spotterL1.Y;

    %
    spotsmartL1.zhab = dataAll(1).spotsmartL1.zhab;
    spotsmartL1.dt = dataAll(1).spotsmartL1.dt;
    spotsmartL1.gapTH = dataAll(1).spotsmartL1.gapTH;


    %%

    % % list_parentfields = {'spectra', 'displacement', 'location', 'SST'};

    list_fields_cat = {'dtime'; 'pressure'};
    
    %% Loop over files and merge in data structure spotterL1

    %
    for i2 = 1:Nfiles
        %
        if i2==1
            %
            for i3 = 1:length(list_fields_cat)
                spotsmartL1.(list_fields_cat{i3}) = dataAll(ind_sort(i2)).spotsmartL1.(list_fields_cat{i3})(:);
            end

        %
        else

            %
            ladddata_aux = ~isnan(dataAll(ind_sort(i2)).spotsmartL1.pressure);

            %
            spotsmartL1.pressure(ladddata_aux) = dataAll(ind_sort(i2)).spotsmartL1.pressure(ladddata_aux);

        end

    end

    %% Save merged data
    
    save(fullfile(dir_output, ['roxsi_smartmooring_L1_' moordID_merge{i1} '.mat']), 'spotsmartL1', '-v7.3')


    %%
    
    clear spotsmartL1 dataAll
    
end


%% Copy Smart Mooring files that don't need to be merged
% from the gridded to the merged folder -- just so there
% is a data file for each mooring (for consistency).


%
copyfile(fullfile(dir_data, 'roxsi_smartmooring_L1_E01sp_1851_gridded.mat'), ...
         fullfile(dir_output, 'roxsi_smartmooring_L1_E01.mat'))
%
copyfile(fullfile(dir_data, 'roxsi_smartmooring_L1_E02sp_1859_gridded.mat'), ...
         fullfile(dir_output, 'roxsi_smartmooring_L1_E02.mat'))
%
copyfile(fullfile(dir_data, 'roxsi_smartmooring_L1_E05sp_1853_gridded.mat'), ...
         fullfile(dir_output, 'roxsi_smartmooring_L1_E05.mat'))
%
copyfile(fullfile(dir_data, 'roxsi_smartmooring_L1_E08sp_1852_gridded.mat'), ...
         fullfile(dir_output, 'roxsi_smartmooring_L1_E08.mat'))
%
copyfile(fullfile(dir_data, 'roxsi_smartmooring_L1_E10sp_1848_gridded.mat'), ...
         fullfile(dir_output, 'roxsi_smartmooring_L1_E10.mat'))
%
copyfile(fullfile(dir_data, 'roxsi_smartmooring_L1_E11sp_1860_gridded.mat'), ...
         fullfile(dir_output, 'roxsi_smartmooring_L1_E11.mat'))
%
copyfile(fullfile(dir_data, 'roxsi_smartmooring_L1_E13sp_1849_gridded.mat'), ...
         fullfile(dir_output, 'roxsi_smartmooring_L1_E13.mat'))





