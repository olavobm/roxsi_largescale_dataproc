function [Lat, Long] = UTMtoll(UTMEasting, UTMNorthing, ZoneNumber)
%% [Lat, Long] = UTMTOLL(UTMEasting, UTMNorthing, ZoneNumber)
%
%   inputs
%       - UTMEasting: UTM easting coordinate.
%       - UTMNorthing: UTM northing coordinate.
%       - ZoneNumber: character referencing a single UTM zone number.
%
%   outputs
%       - Lat: latitude.
%       - Long: longitude.
%
%
% I have just copied this code. There is this unused variable
% "NorthernHemisphere = 1". So I don't know if there something
% in this function restricts it to the Northern Hemisphere.
%


%% Define constanst

%
deg2rad = pi/180;
rad2deg = 180/pi;
ellipsoidName = 'WGS-84';
a = 6378137;
eccSquared = 0.00669438;
k0 = 0.9996;


%% Get the numeric part of the UTM zone and convert to double

%
ZoneNumber_double = str2double(ZoneNumber(1:2));


%% Do the conversion from UTM to geographical coordinate

%
e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared));

x = UTMEasting - 500000.0; %%//remove 500,000 meter offset for longitude
y = UTMNorthing;

NorthernHemisphere = 1;  %%//point is in northern hemisphere

LongOrigin = (ZoneNumber_double - 1)*6 - 180 + 3;  %%//+3 puts origin in middle of zone

eccPrimeSquared = (eccSquared)/(1-eccSquared);

M = y ./ k0;
mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64-5*eccSquared*eccSquared*eccSquared/256));

phi1Rad = mu	+ (3*e1/2-27*e1*e1*e1/32)*sin(2*mu) ...
			+ (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)...
			+(151*e1*e1*e1/96)*sin(6*mu);
phi1 = phi1Rad*180/pi;

N1 = a./sqrt(1-eccSquared.*sin(phi1Rad).*sin(phi1Rad));
T1 = tan(phi1Rad).*tan(phi1Rad);
C1 = eccPrimeSquared.*cos(phi1Rad).*cos(phi1Rad);
R1 = a.*(1-eccSquared)./((1-eccSquared.*sin(phi1Rad).*sin(phi1Rad)).^1.5);
D = x./(N1.*k0);

Lat = phi1Rad - (N1.*tan(phi1Rad)./R1).*(D.*D./2-(5+3.*T1+10.*C1-4.*C1.*C1-9.*eccPrimeSquared).*D.*D.*D.*D./24 ...
                    +(61+90.*T1+298.*C1+45.*T1.*T1-252.*eccPrimeSquared-3.*C1.*C1).*D.*D.*D.*D.*D.*D./720);
Lat = Lat * 180/pi;

Long = (D-(1+2.*T1+C1).*D.*D.*D./6+(5-2.*C1+28.*T1-3.*C1.*C1+8.*eccPrimeSquared+24.*T1.*T1) ...
                     .*D.*D.*D.*D.*D./120)./cos(phi1Rad);
Long = LongOrigin + Long * 180/pi;



