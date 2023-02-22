%% Merge L1 Spotter data from the two spotters that were deployed
% at the B01 mooring site (the first had issues and was replaced)

clear
close all

%%
% -------------------------------------------
% ------------- DEFINE VARIABLES ------------
% -------------------------------------------

%%

%
% dir_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';
dir_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Level1/not_gridded/';


%
% moordID_merge = {'B01'};
% moordID_merge = {'E09'};
moordID_merge = {'B01', 'E07', 'E09'};

%
name_string_beginning = 'roxsi_spotter_L1_';

%
% % dir_output = dir_data;
dir_output = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Level1/not_gridded/';


%%
% -------------------------------------------
% ---------------- ..... ----------------
% -------------------------------------------

%%

list_parentfields = {'spectra', 'displacement', 'location', 'SST'};

%%


for i1 = 1:length(moordID_merge)
    
    %% Check if a merged file already exists
    
    %
    list_checkmerged = dir(fullfile(dir_data, [name_string_beginning, moordID_merge{i1} '.mat']));
    %
    if ~isempty(list_checkmerged)
        warning(['Merged Spotter data file already exists for mooring ' moordID_merge{i1} '. Skipping to next one.'])
        continue
    end
    
    %%
    
    %
    list_dirSpotter = dir(fullfile(dir_data, [name_string_beginning, moordID_merge{i1} '*.mat']));
    
    %
    Nfiles = length(list_dirSpotter);

    %
    if Nfiles<2
        error(['There are NOT at least 2 Spotter data files ' ...
               'at mooring ID ' moordID_merge{i1} ' to be merged.'])
    end


    %%

    %
    time_edges = NaT(Nfiles, 2);

    %
    for i2 = 1:Nfiles
       %
       dataAll(i2) = load(fullfile(list_dirSpotter(i2).folder, list_dirSpotter(i2).name));

       %
       if i2==1
           time_edges.TimeZone = dataAll(i2).spotterL1.displacement.dtime.TimeZone;
       end

       %
       time_edges(i2, :) = dataAll(i2).spotterL1.displacement.dtime([1, end]);

    end

    %
    [~, ind_sort] = sort(time_edges(:, 1)); 

    
    %%

    %
    spotterL1.mooringID = moordID_merge{i1};

    %
    spotterL1.SN = [];
    for i2 = 1:Nfiles
        %
        if i2~=Nfiles
            spotterL1.SN = [spotterL1.SN, dataAll(ind_sort(i2)).spotterL1.SN ' + '];
        else 
            spotterL1.SN = [spotterL1.SN, dataAll(ind_sort(i2)).spotterL1.SN];
        end
    end

    %
    spotterL1.site = dataAll(1).spotterL1.site;
    spotterL1.latitude = dataAll(1).spotterL1.latitude;
    spotterL1.longitude = dataAll(1).spotterL1.longitude;
    spotterL1.X = dataAll(1).spotterL1.X;
    spotterL1.Y = dataAll(1).spotterL1.Y;

    
    %% Loop over files and merge in data structure spotterL1

    %
    for i2 = 1:Nfiles
        %
        if i2==1
            %
            for i3 = 1:length(list_parentfields)
                %
                if isfield(dataAll(ind_sort(i2)).spotterL1, list_parentfields{i3})
                    spotterL1.(list_parentfields{i3}) = dataAll(ind_sort(i2)).spotterL1.(list_parentfields{i3});
                end
            end

        %
        else

            %
            for i3 = 1:length(list_parentfields)

                if ~isfield(dataAll(ind_sort(i2)).spotterL1, list_parentfields{i3})
                    continue
                end

                %
                list_fields_aux = fieldnames(dataAll(ind_sort(i2)).spotterL1.(list_parentfields{i3}));

                %
                for i4 = 1:length(list_fields_aux)
                    %
                    if size(dataAll(ind_sort(i2)).spotterL1.(list_parentfields{i3}).(list_fields_aux{i4}), 1) == ...
                       length(dataAll(ind_sort(i2)).spotterL1.(list_parentfields{i3}).dtime)
                        %
                        spotterL1.(list_parentfields{i3}).(list_fields_aux{i4}) = ...
                                [spotterL1.(list_parentfields{i3}).(list_fields_aux{i4}); ...
                                 dataAll(ind_sort(i2)).spotterL1.(list_parentfields{i3}).(list_fields_aux{i4})];
                    end
                end
            end

        end

    end



    %% Save merged data

    %
    disp(['--- Saving merged Spotter data for mooring ' moordID_merge{i1} ' ---'])
    
    %
    save(fullfile(dir_output, [name_string_beginning, moordID_merge{i1} '.mat']), 'spotterL1', '-v7.3')


    %%
    
    clear dataAll spotterL1    % otherwise, unwanted fields will be copied
    
end








