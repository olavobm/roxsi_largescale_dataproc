%% Simple script to dheck deployment locations
% in ROXSI2022_mooring_locations.csv

clear
close all


%% Load and parse variables

%
file_fid = fopen('./ROXSI2022_mooring_locations.csv');

loaded_vars = textscan(file_fid, '%s%s %f%f %f%f %s %s', ...
                                 'Headerlines', 1, 'Delimiter', ',');

%
mooringsLocs.mooringID = loaded_vars{1};
mooringsLocs.array = loaded_vars{2};

%
mooringsLocs.planned.latitude = loaded_vars{3};
mooringsLocs.planned.longitude = loaded_vars{4};

%
mooringsLocs.actual.latitude = loaded_vars{5};
mooringsLocs.actual.longitude = loaded_vars{6};

%
mooringsLocs.mooringType = loaded_vars{7};
mooringsLocs.GPSinfo = loaded_vars{8};


%% Replace flag values with NaNs

%
mooringsLocs.planned.latitude(mooringsLocs.planned.latitude==9999999999) = NaN;
mooringsLocs.planned.longitude(mooringsLocs.planned.latitude==9999999999) = NaN;
%
mooringsLocs.actual.latitude(mooringsLocs.planned.latitude==9999999999) = NaN;
mooringsLocs.actual.longitude(mooringsLocs.planned.latitude==9999999999) = NaN;


%% Select whether to look at planned or actual locations

pltwhich = 'planned';
% pltwhich = 'actual';


%% All large-scale moorings (no ISPAR for now)

%
figure
    plot(mooringsLocs.(pltwhich).longitude(1:end-1), ...
         mooringsLocs.(pltwhich).latitude(1:end-1), '.k', 'MarkerSize', 42)

    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'DataAspectRatio', [1, 1, 1], ...
             'XGrid', 'on', 'YGrid', 'on')
    %
    xlabel('Longitude', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Latitude', 'Interpreter', 'Latex', 'FontSize', 18)


%% Just Asilomar

%
figure
    plot(mooringsLocs.(pltwhich).longitude(strcmp(mooringsLocs.array, 'Asilomar')), ...
         mooringsLocs.(pltwhich).latitude(strcmp(mooringsLocs.array, 'Asilomar')), ...
         '.k', 'MarkerSize', 42)

    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'DataAspectRatio', [1, 1, 1], ...
             'XGrid', 'on', 'YGrid', 'on')
    %
    xlabel('Longitude', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Latitude', 'Interpreter', 'Latex', 'FontSize', 18)


%% Just China Rock (no ISPAR for now, because I
% haven't got the actual ISPAR location yet)

%
lplot_aux = strcmp(mooringsLocs.array, 'ChinaRock') & ...
            ~strcmp(mooringsLocs.mooringID, 'ISPAR');

% Highlight one array
lplotHL_aux = strncmp(mooringsLocs.mooringID, 'B', 1);

%
figure
    plot(mooringsLocs.(pltwhich).longitude(lplot_aux), ...
         mooringsLocs.(pltwhich).latitude(lplot_aux), ...
         '.k', 'MarkerSize', 42)
    hold on
    plot(mooringsLocs.(pltwhich).longitude(lplotHL_aux), ...
         mooringsLocs.(pltwhich).latitude(lplotHL_aux), ...
         '.r', 'MarkerSize', 42)

    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'DataAspectRatio', [1, 1, 1], ...
             'XGrid', 'on', 'YGrid', 'on')
    %
    xlabel('Longitude', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Latitude', 'Interpreter', 'Latex', 'FontSize', 18)


%% Check distances between smart moorings

%
lsmartmoorings = strncmp(mooringsLocs.mooringID, 'E', 1);
for i = 1:length(lsmartmoorings)
    %
    if lsmartmoorings(i)
        %
        if strcmp(mooringsLocs.mooringID{i}(4), 'a')
            lsmartmoorings(i) = false;
        end
    end
end

%
ind_smartmoorings = find(lsmartmoorings);

%
1000*deg2km(distance(mooringsLocs.(pltwhich).latitude(ind_smartmoorings(1:end-1)), ...
                     mooringsLocs.(pltwhich).longitude(ind_smartmoorings(1:end-1)), ...
                     mooringsLocs.(pltwhich).latitude(ind_smartmoorings(2:end)), ...
                     mooringsLocs.(pltwhich).longitude(ind_smartmoorings(2:end))))
