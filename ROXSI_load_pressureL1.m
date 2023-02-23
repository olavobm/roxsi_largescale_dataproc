function dataL1 = ROXSI_load_pressureL1(mooringID, dtimelims)
%% dataL1 = ROXSI_LOAD_PRESSUREL1(mooringID, dtimelims)
%
%   inputs
%       - mooringID:
%       - dtimelims:
%
%   outputs
%       - dataL1:
%
%
%
%
%
% Olavo Badaro Marques


%% Make sure mooringID is a cell array

%
if isstring(mooringID)
    mooringID = cellstr(mooringID);
%
elseif ischar(mooringID)
    mooringID = {mooringID};
end


%%

%
if exist('dtimelims', 'var')

    %
    if isempty(dtimelims.TimeZone)
        dtimelims.TimeZone = 'America/Los_Angeles';
    end    
    %
    ltimesub = true;
    
%
else
    ltimesub = false;
end


%% 

%
dir_parent = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/';
%
% fullfile(dir_parent)


%%

listInstr.SoloD.dirdata = fullfile(dir_parent, 'SoloD_Level1');
listInstr.Signature.dirdata = fullfile(dir_parent, 'Signature_Level1');
listInstr.Aquadopp.dirdata = fullfile(dir_parent, 'Aquadopp_Level1');

%
% listInstr.SmartMooring.dirdata = fullfile(dir_parent, 'Spotter_Smart_Level1', 'gridded');
% Uses alternative folder until I merge with the main data repository
listInstr.SmartMooring.dirdata = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Smart_Level1/merged/';


%%

%
listInstr.SoloD.list_moorings = ["A02", "A04", "A05", "A06", "A07", "A09", ...
                                 "B04", "B06", "B07", "B09", "B12", "B14", "B16", "B18", ...
                                 "C02", "C03", "C05", "C06", "C07", "C08", "C09", ...
                                 "D01", "D02", "X07", "X08", "X09", "X10", "X12", "X14"];

%
listInstr.Signature.list_moorings = ["A01", "B10", "B13", "B15", "B17", "C01", "X05", "X11"];

%
listInstr.Aquadopp.list_moorings = ["B08", "B11", "E03", "E04", "E06", "E12", ...
                      "F01", "F02", "F03", "F04", "F05", "X06", "X13", ];

%
listInstr.SmartMooring.list_moorings = ["E01", "E02", "E05", "E07", "E08", "E09", "E10", "E11", "E13"];


%%

list_mooringTypes = fieldnames(listInstr);


%%
% ----------------------------------------
% ----------------------------------------
% ----------------------------------------


%%

% % % dataL1.


%%

%
for i1 = 1:length(mooringID)
    
    %
    for i2 = 1:length(list_mooringTypes)
    
        %
        if any(strcmp(listInstr.(list_mooringTypes{i2}).list_moorings, mooringID{i1}))
            %
            ind_type = i2;
            %
            filedir = dir(fullfile(listInstr.(list_mooringTypes{i2}).dirdata, ...
                                   ['*' mooringID{i1} '*.mat']));
        end
    end
    
    
    %%
    
    %
    data_aux = load(fullfile(filedir.folder, filedir.name));
    %
    namefield_aux = fieldnames(data_aux);
    
    %
    if length(namefield_aux)~=1
        error('Expected a structure variable!!!') 
    end
    data_aux = data_aux.(namefield_aux{1});
    
    
    %%
    
    %
    dataL1.(mooringID{i1}).sensor = list_mooringTypes{ind_type};
    
    %
    if strcmp(list_mooringTypes{ind_type}, "SoloD")
        
        %
        time_aux = data_aux.time_dnum;
        pressure_aux = data_aux.Pwater;
        
        %
%         dataL1.SN
        dataL1.(mooringID{i1}).latitude = data_aux.latitude;
        dataL1.(mooringID{i1}).longitude = data_aux.longitude;
        dataL1.(mooringID{i1}).X = data_aux.X;
        dataL1.(mooringID{i1}).Y = data_aux.Y;
        dataL1.(mooringID{i1}).site = data_aux.site;

        %
        dataL1.(mooringID{i1}).zhab = 0.09;
        
        %
        time_aux = datetime(time_aux, 'ConvertFrom', 'datenum');
        time_aux.TimeZone = 'America/Los_Angeles';
        
        %
        dataL1.(mooringID{i1}).dtime = time_aux(:);
        dataL1.(mooringID{i1}).pressure = (1025*9.8)*pressure_aux(:)./1e4;    % Falk saves pressure in meters. Convert to dbar
        
        % --------------------------------------------------
        % --------------------------------------------------
        % Interpolate pressure in datetime
        
        % In seconds
        gapTH = seconds(4);
        dt_grid = seconds(0.5);
        
        % ------------------------------
        % Start and end time grid on the whole second
        dtime_start = datetime(year(dataL1.(mooringID{i1}).dtime(1)), month(dataL1.(mooringID{i1}).dtime(1)), day(dataL1.(mooringID{i1}).dtime(1)), ...
                               hour(dataL1.(mooringID{i1}).dtime(1)), minute(dataL1.(mooringID{i1}).dtime(1)), 1 + second(dataL1.(mooringID{i1}).dtime(1)));
        %
        dtime_end = datetime(year(dataL1.(mooringID{i1}).dtime(end)), month(dataL1.(mooringID{i1}).dtime(end)), day(dataL1.(mooringID{i1}).dtime(end)), ...
                             hour(dataL1.(mooringID{i1}).dtime(end)), minute(dataL1.(mooringID{i1}).dtime(end)), -1 + second(dataL1.(mooringID{i1}).dtime(end)));
        
        %
        dtime_start.TimeZone = time_aux.TimeZone;
        dtime_end.TimeZone = time_aux.TimeZone;
        
        %
        dtime_grid = dtime_start : dt_grid : dtime_end;
        dtime_grid = dtime_grid(:);
        
        % ------------------------------
        %
        pinterp_aux = interp1_skipgaps(dataL1.(mooringID{i1}).dtime, ...
                                       dataL1.(mooringID{i1}).pressure, ...
                                       gapTH, dtime_grid);
        %
        dataL1.(mooringID{i1}).dtime = dtime_grid(:);
        dataL1.(mooringID{i1}).pressure = pinterp_aux(:);
        

	%
    else
        
        %
        dataL1.(mooringID{i1}).dtime = data_aux.dtime(:);
        dataL1.(mooringID{i1}).pressure = data_aux.pressure(:);
        
    end

    
    %%
    
    if ltimesub
        %
        lintime_aux = (dataL1.(mooringID{i1}).dtime >= dtimelims(1)) & ...
                      (dataL1.(mooringID{i1}).dtime <= dtimelims(2));
        %
        dataL1.(mooringID{i1}).dtime = dataL1.(mooringID{i1}).dtime(lintime_aux);
        dataL1.(mooringID{i1}).pressure = dataL1.(mooringID{i1}).pressure(lintime_aux);
    end

end



