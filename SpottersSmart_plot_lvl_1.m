%% Plot all level 1 data from smart moorings

clear
close all


%%

%
dir_data_L1 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/code_proc';


%%

% % %  Select some files
% % list_files = ["smart_mooring_E01sp_1851_L1.mat"; ...
% %               "smart_mooring_E02sp_1859_L1.mat"; ...
% %               "smart_mooring_E08sp_1852_L1.mat"; ...
% %               "smart_mooring_E10sp_1848_L1.mat"];

% Or get all L1 files in the folder
%
dir_all_L1 = dir(fullfile(dir_data_L1, '*_L1.mat'));
%
list_files = cell(1, length(dir_all_L1));

%
for i = 1:length(list_files)
    %
    list_files{i} = dir_all_L1(i).name;
end


%% Loop over file names and load data

for i = 1:length(list_files)

    %
    data_loaded = load(fullfile(dir_data_L1, list_files{i}));

    %
    smartMoorlvl1.([list_files{i}(15:17) '_' list_files{i}(21:24)]) = data_loaded.spotsmart;

end

%
list_fields = fieldnames(smartMoorlvl1);


%% Plot timeseries of pressure from all Smart moorings

%
figure
    hold on
    %
    for i = 1:length(list_fields)
        plot(smartMoorlvl1.(list_fields{i}).dtime, ...
             smartMoorlvl1.(list_fields{i}).pressure + 2.*(i-1), '.-')
    end

    %
    ylabel('Water pressure [dbar]', 'Interpreter', 'Latex', 'FontSize', 22)
    %
    grid on
    set(gca, 'FontSize', 16)
    %
    set(gcf, 'units', 'normalized')
    set(gcf, 'Position', [0.2, 0.2, 0.6, 0.4])

    

%% Plot clock difference

% % %
% % for i = 1:length(list_fields)
% %     %
% %     figure
% %         %
% %         hold on
% %         plot(datetime((datenum(smartMoorlvl1.(list_fields{i}).dtime(1:end-1)) + ...
% %                        datenum(smartMoorlvl1.(list_fields{i}).dtime(2:end)))./2, 'ConvertFrom', 'datenum'), ...
% %              24*3600*diff(datenum(smartMoorlvl1.(list_fields{i}).dtime)), '.-k')
% % 
% %         %
% %         plot(xlim, [0, 0], '-k')
% % 
% %     %
% %     set(gca, 'FontSize', 16, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
% %              'YLim', [-0.5, 4])
% %     %
% %     ylabel('[s]', 'Interpreter', 'Latex', 'FontSize', 22)
% %     title(['Smart mooring ' list_fields{i} ': diff(time)'], 'Interpreter', 'Latex', 'FontSize', 22)
% %     %
% %     set(gcf, 'units', 'normalized')
% %     set(gcf, 'Position', [0.2, 0.2, 0.5, 0.3])
% % end

    
