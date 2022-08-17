function [x, y] = ROXSI_lltoxy(latitude, longitude, strarray)
%% [x, y] = ROXSI_LLTOXY(latitude, longitude, strarray)
%
%   inputs
%       - latitude: a vector (or scalar) of latitudes.
%       - longitude: a vector of longitudes of the same
%                    dimension as latitude.
%       - strarray (optional):
%
%   outputs
%       - x:
%       - y:
%
%
% ROXSI_LLTOXY.m.
%
% IMPORTANT NOTE (!!!): The *.mat file which defines the xy
% ROXSI grids must be saved in the same directory as this
% function (and the name/format must be what this function expects).
% Also, the string in this function that specified the arrays
% MUST MATCH the strings in the ROXSI grids.
%
% Olavo Badaro Marques.



%% This is only used if input strarray is not given,
% such that this function needs to figure out which
% array corresponds to each (lat, lon)

%
boxlimits.Asilomar.latitude = 36 + [((36/60) + (58/3600)), ((38/60) + (05/3600))];
boxlimits.Asilomar.longitude = -121 - [1, (55/60)];

%
boxlimits.ChinaRock.latitude = 36 + [(34/60), ((36/60) + (57/3600))];
boxlimits.ChinaRock.longitude = -121 - [1, (55/60)];

%
list_arrays = fieldnames(boxlimits);


%% Based on the coordinate, identify which grid the coordinates belong to

%
if ~exist('strarray', 'var')
    
    %
    lin_Asilomar = (latitude >= boxlimits.Asilomar.latitude(1)) & ...
                   (latitude <= boxlimits.Asilomar.latitude(2)) & ...
                   (longitude >= boxlimits.Asilomar.longitude(1)) & ...
                   (longitude <= boxlimits.Asilomar.longitude(2));

    %
    lin_ChinaRock = (latitude >= boxlimits.ChinaRock.latitude(1)) & ...
                    (latitude <= boxlimits.ChinaRock.latitude(2)) & ...
                    (longitude >= boxlimits.ChinaRock.longitude(1)) & ...
                    (longitude <= boxlimits.ChinaRock.longitude(2));

    %
    strarray = strings(size(latitude));

    %
    strarray(lin_Asilomar) = "Asilomar";
    strarray(lin_ChinaRock) = "ChinaRock";

    %
    Ndiffarrays = length(unique(strarray));

    %
    if Ndiffarrays==1
        strarray = unique(strarray);
    end

%
else
    % Just in case, make sure this is a 1x1 string array
    %
% %     strarray

    %
    Ndiffarrays = 1;
end


%% Get grid information

%
fcnfullpathname = mfilename('fullpath');
indslash = strfind(fcnfullpathname, '/');
%
fcnfolder = fcnfullpathname(1:indslash(end));

%
load(fullfile(fcnfolder, 'ROXSI_xygrids.mat'), 'roxsigrid');

% % % This function expects the loaded variable to look like:
% % roxsigrid.Asilomar
% % roxsigrid.ChinaRock


%% Convert from (latitude, longitude) to UTM coordinates

%
[UTMEasting, UTMNorthing, ~] = lltoUTM(latitude, longitude);


%% Convert from UTM coordinates to (x, y) grid

%
if Ndiffarrays==1

    %
    [originEasting, ...
     originNorthing] = lltoUTM(roxsigrid.(strarray).latref, ...
                               roxsigrid.(strarray).lonref);

    %
    angle_rotation = 270 - roxsigrid.(strarray).angleref;

    %
    [x, y] = UTMtoxy(UTMEasting, UTMNorthing, originEasting, originNorthing, angle_rotation);

%
else

    %
    x = NaN(size(latitude));
    y = x;

    %
    for i = 1:Ndiffarrays

        %
        lin_array_aux = (latitude >= boxlimits.(list_arrays{i}).latitude(1)) & ...
                        (latitude <= boxlimits.(list_arrays{i}).latitude(2)) & ...
                        (longitude >= boxlimits.(list_arrays{i}).longitude(1)) & ...
                        (longitude <= boxlimits.(list_arrays{i}).longitude(2));

        %
        [originEasting_aux, ...
         originNorthing_aux] = lltoUTM(roxsigrid.(list_arrays{i}).latref, ...
                                       roxsigrid.(list_arrays{i}).lonref);

        %
        angle_rotation_aux = 270 - roxsigrid.(list_arrays{i}).angleref;

        %
        [x_aux, y_aux] = UTMtoxy(UTMEasting(lin_array_aux), UTMNorthing(lin_array_aux), ...
                                 originEasting_aux, originNorthing_aux, angle_rotation_aux);

        %
        x(lin_array_aux) = x_aux;
        y(lin_array_aux) = y_aux;

    end

end
