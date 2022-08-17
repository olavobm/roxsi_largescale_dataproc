function [x, y] = UTMtoxy(easting, northing, originEasting, originNorthing, rotCCW)
%% [x, y] = UTMTOXY(easting, northing, originEasting, originNorthing, rotCCW)
%
%   inputs
%       - easting:
%       - northing:
%       - originEasting:
%       - originNorthing:
%       - rotCCW: angle (in degrees) of the counterclockwise rotation
%                 of the coordinate system (from UTM to xy).
%
%   outputs
%       - x:
%       - y:
%
%

%
% Olavo Badaro Marques.



%% Shift data relative to origin

%
x_aux = easting - originEasting;
y_aux = northing - originNorthing;


%% Rotate coordinates to the desired angle

% Compute x coordinate
x = x_aux*cosd(rotCCW) + y_aux*sind(rotCCW);

% Compute y coordinate
y = y_aux*cosd(rotCCW) - x_aux*sind(rotCCW);
