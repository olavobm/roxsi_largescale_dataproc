function [latitude, longitude] = ROXSI_xytoll(strarray, x, y)
%% [latitude, longitude] = ROXSI_XYTOLL(strarray, x, y)
%
%   inputs
%       - strarray (optional):
%       - x: 
%       - y: 
%
%   outputs
%       - latitude:
%       - longitude:
%
%
% ROXSI_XYTOLL.m converts local (x, y) coordinates in the strarray
% coordinate system to latitude and longitude. (x, y) are in meters
% relative to an origin and rotated to have the x-axis increasing
% onshore. ROXSI_XYTOLL.m is a high-level function and the conversion
% is computed by two other functions.
%
% Based on code by B. Woodward 14 Nov, 2008.
%
% Olavo Badaro Marques.
%
% See also:
%   ROXSI_xytoUTM.m
%   UTMtoll.m



%% Convert (x, y) to UTM

[easting, northing] = ROXSI_xytoUTM(strarray, x, y);


%% Convert data N,E to lat,lon

[latitude, longitude] = UTMtoll(easting, northing, '10S');

