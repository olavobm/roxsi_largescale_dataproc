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
dir_data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/Spotter_Level1/';

%
% moordID_merge = {'B01'};
moordID_merge = {'E09'};

%
name_string_beginning = 'roxsi_spotter_L1_';

%
% % dir_output = dir_data;
dir_output = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Level1/';


%%
% -------------------------------------------
% ---------------- LOAD DATA ----------------
% -------------------------------------------


%%

%
list_dirSpotter = dir(fullfile(dir_data, [name_string_beginning, moordID_merge{1} '*.mat']));

%
Nfiles = length(list_dirSpotter);

%
if Nfiles<2
	error(['There are NOT at least 2 Spotter data files ' ...
           'at mooring ID ' moordID_merge{1} ' to be merged.'])
end


%%
% -------------------------------------------
% ----------------MERGE DATA ----------------
% -------------------------------------------

%%

%
time_edges = NaT(Nfiles, 2);

%
for i = 1:Nfiles
   %
   dataAll(i) = load(fullfile(list_dirSpotter(i).folder, list_dirSpotter(i).name));
    
   %
   if i==1
       time_edges.TimeZone = dataAll(i).spotterL1.displacement.dtime.TimeZone;
   end
   
   %
   time_edges(i, :) = dataAll(i).spotterL1.displacement.dtime([1, end]);
   
end


%
[~, ind_sort] = sort(time_edges(:, 1)); 


%%

%
spotterL1.mooringID = moordID_merge{1};

%
spotterL1.SN = [];
for i = 1:Nfiles
    %
    if i~=Nfiles
        spotterL1.SN = [spotterL1.SN, dataAll(ind_sort(i)).spotterL1.SN ' + '];
    else 
        spotterL1.SN = [spotterL1.SN, dataAll(ind_sort(i)).spotterL1.SN];
    end
end

%
spotterL1.site = dataAll(1).spotterL1.site;
spotterL1.latitude = dataAll(1).spotterL1.latitude;
spotterL1.longitude = dataAll(1).spotterL1.longitude;
spotterL1.X = dataAll(1).spotterL1.X;
spotterL1.Y = dataAll(1).spotterL1.Y;


%%

list_parentfields = {'spectra', 'displacement', 'location', 'SST'};



%% Loop over files and merge in data structure spotterL1

%
for i1 = 1:Nfiles
    %
    if i1==1
        %
        for i2 = 1:length(list_parentfields)
            spotterL1.(list_parentfields{i2}) = dataAll(ind_sort(i1)).spotterL1.(list_parentfields{i2});
        end
        
	%
    else
        
        %
        for i2 = 1:length(list_parentfields)
            %
            list_fields_aux = fieldnames(dataAll(ind_sort(i1)).spotterL1.(list_parentfields{i2}));
            
            %
            for i3 = 1:length(list_fields_aux)
                %
                if size(dataAll(ind_sort(i1)).spotterL1.(list_parentfields{i2}).(list_fields_aux{i3}), 1) == ...
                   length(dataAll(ind_sort(i1)).spotterL1.(list_parentfields{i2}).dtime)
                    %
                    spotterL1.(list_parentfields{i2}).(list_fields_aux{i3}) = ...
                            [spotterL1.(list_parentfields{i2}).(list_fields_aux{i3}); ...
                             dataAll(ind_sort(i1)).spotterL1.(list_parentfields{i2}).(list_fields_aux{i3})];
                end
            end
        end
        
    end

end



%% Save merged data

%
save(fullfile(dir_output, [name_string_beginning, moordID_merge{1} '.mat']), 'spotterL1', '-v7.3')





