function [easting, northing] = ROXSI_xytoUTM(strarray, x, y)
%% [easting, northing] = ROXSI_XYTOUTM(strarray, x, y)
%
%   inputs
%       - strarray (optional):
%       - x: 
%       - y: 
%
%   outputs
%       - easting:
%       - northing:
%
%
% ROXSI_XYTOUTM.m converts local (x, y) coordinates in the strarray
% coordinate system to easting and northing UTM coordinates. (x, y)
% are in meters relative to an origin and rotated to have the x-axis
% increasing onshore.
%
% IMPORTANT NOTE (!!!): The *.mat file which defines the xy
% ROXSI grids must be saved in the same directory as this
% function (and the name/format must be what this function expects).
% Also, the string in this function that specified the arrays
% MUST MATCH the strings in the ROXSI grids.
%
% Based on code by B. Woodward 14 Nov, 2008.
%
% Olavo Badaro Marques.


%% Check strarray

%
if ~(strcmp(strarray, "Asilomar") || strcmp(strarray, "ChinaRock"))
    %
    error('Specified grid in input was not recognized!!!')
end


%% Define CONSTANTS

%
ZoneNumber = '10S';    % UTM zone for much of California (can use just 10 as well)


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


%% Rotate (x, y) so that they are along the (easting, northing) grid

%
rotCCW = 270 - roxsigrid.(strarray).angleref;
rotCCW = -rotCCW;

%
% Compute x coordinate
x_rotated = x*cosd(rotCCW) + y*sind(rotCCW);

% Compute y coordinate
y_rotated = y*cosd(rotCCW) - x*sind(rotCCW);


%% Convert Origin lat,lon to UTM

[OriginE, OriginN, UTMZone] = lltoUTM(roxsigrid.(strarray).latref, ...
                                      roxsigrid.(strarray).lonref);
% % UTMZone == '10S'

%% Shift data relative back to UTM N,E:

easting = x_rotated + OriginE;
northing = y_rotated + OriginN;



