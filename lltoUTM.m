function [UTMEasting, UTMNorthing, UTMZone] = lltoUTM(latitude, longitude)
%% [UTMEasting, UTMNorthing, UTMZone] = LLTOUTM(latitude, longitude)
%
%   inputs
%       - latitude = latitude in fraction degrees
%       - longitude = longitude in fraction degrees (West of Grenwich are negative!).
%
%   outputs
%       - UTMEasting: easting coordinate.
%       - UTMNorthing: northing coordinate.
%       - UTMZone: character variable with with the UTM zone
%                  of ONLY THE FIRST coordinate in the inputs.
%
%
%
% LLTOUTM.m converts (latitude, longitude) data to UTM coordinates.
% Uses subroutine written by Chuck Gantz downloaded from
% the GPSy web site.
%
% Note: longitudes west of greenwhich must be negative.
%
%
% Rewritten/copied on 16, August 2022. O. Marques
%
% Consolidated  29,August 2006.  B. Woodward
% 14 April 1999. T. C. Lippmann
%
%
%
%
%
%	refellip = ref. ellipsoid identifier (see refellip_menu.m)
%	a = Equatorial Radius (optional)
%	eccSquared = eccentricity squared (optional)
%
%	if a and eccSquared are included, then refellip is ignored and a & eccSquared are
% 	used for the radius and eccentricity squared
%


%% %%%%%%%%  CONSTANTS  %%%%%%%%%

%
deg2rad = pi/180;
rad2deg = 180/pi;
ellipsoidName = 'WGS-84';
a = 6378137;
eccSquared = 0.00669438;
k0 = 0.9996;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

%
LatRad = latitude * deg2rad;
LongRad = longitude * deg2rad;

%
LongOrigin = zeros(size(longitude));
%
dy = find(longitude > -6 & longitude <= 0);
LongOrigin(dy) = -3; 
%
dy = find(longitude < 6 & longitude > 0);
LongOrigin(dy) = 3; 
%
dy = find(abs(longitude) >= 6);
LongOrigin(dy) = sign(longitude(dy)).*floor(abs(longitude(dy))/6)*6 + 3*sign(longitude(dy));

%
LongOriginRad = LongOrigin * deg2rad;


%% Compute the UTM Zone and Grid from the latitude and longitude

%
UTMZone = sprintf('%d%c', floor((longitude(1) + 180)/6) + 1, UTMLetterDesignator(latitude(1)));

%
eccPrimeSquared = (eccSquared)./(1-eccSquared);

%
N = a./sqrt(1-eccSquared.*sin(LatRad).*sin(LatRad));
T = tan(LatRad).*tan(LatRad);
C = eccPrimeSquared.*cos(LatRad).*cos(LatRad);
A = cos(LatRad).*(LongRad-LongOriginRad);

%
M = a.*((1 - eccSquared/4 - 3*eccSquared*eccSquared/64- 5*eccSquared*eccSquared*eccSquared/256).*LatRad ...
    - (3*eccSquared/8 + 3*eccSquared*eccSquared/32 + 45*eccSquared*eccSquared*eccSquared/1024).*sin(2*LatRad) ...
    + (15*eccSquared*eccSquared/256 + 45*eccSquared*eccSquared*eccSquared/1024).*sin(4*LatRad) ...
    - (35*eccSquared*eccSquared*eccSquared/3072).*sin(6*LatRad));

%
UTMEasting = (k0.*N.*(A+(1-T+C).*A.*A.*A/6 ...
		+ (5-18.*T+T.*T+72.*C-58.*eccPrimeSquared).*A.*A.*A.*A.*A/120) ...
		+ 500000.0);
%
UTMNorthing = (k0.*(M+N.*tan(LatRad).*(A.*A/2+(5-T+9.*C+4.*C.*C).*A.*A.*A.*A/24 ...
		+ (61-58.*T+T.*T+600.*C-330.*eccPrimeSquared).*A.*A.*A.*A.*A.*A/720)));

%
dy = find(latitude < 0);
%
UTMNorthing(dy) = UTMNorthing(dy) + 10000000.0;    % 10000000 meter offset for southern hemisphere

%
return;
  

%%
function [LetterDesignator] = UTMLetterDesignator(Lat)
%% function [letdes] = UTMLetterDesignator(Lat)
% This routine determines the correct UTM letter designator for the given latitude
% returns 'Z' if latitude is outside the UTM limits of 80N to 80S
%
% //Written by Chuck Gantz- chuck.gantz@globalstar.com

        if((80 >= Lat) & (Lat > 72)) LetterDesignator = 'X';
        else if((72 >= Lat) & (Lat > 64)) LetterDesignator = 'W'; 
        else if((64 >= Lat) & (Lat > 56)) LetterDesignator = 'V'; 
        else if((56 >= Lat) & (Lat > 48)) LetterDesignator = 'U'; 
        else if((48 >= Lat) & (Lat > 40)) LetterDesignator = 'T'; 
        else if((40 >= Lat) & (Lat > 32)) LetterDesignator = 'S'; 
        else if((32 >= Lat) & (Lat > 24)) LetterDesignator = 'R'; 
        else if((24 >= Lat) & (Lat > 16)) LetterDesignator = 'Q'; 
        else if((16 >= Lat) & (Lat > 8)) LetterDesignator = 'P'; 
        else if(( 8 >= Lat) & (Lat > 0)) LetterDesignator = 'N'; 
        else if(( 0 >= Lat) & (Lat > -8)) LetterDesignator = 'M'; 
        else if((-8>= Lat) & (Lat > -16)) LetterDesignator = 'L'; 
        else if((-16 >= Lat) & (Lat > -24)) LetterDesignator = 'K'; 
        else if((-24 >= Lat) & (Lat > -32)) LetterDesignator = 'J'; 
        else if((-32 >= Lat) & (Lat > -40)) LetterDesignator = 'H'; 
        else if((-40 >= Lat) & (Lat > -48)) LetterDesignator = 'G'; 
        else if((-48 >= Lat) & (Lat > -56)) LetterDesignator = 'F'; 
        else if((-56 >= Lat) & (Lat > -64)) LetterDesignator = 'E'; 
        else if((-64 >= Lat) & (Lat > -72)) LetterDesignator = 'D'; 
        else if((-72 >= Lat) & (Lat > -80)) LetterDesignator = 'C'; 
        else LetterDesignator = 'Z'; %%This is here as an error flag 
				     %%to show that the Latitude is 
				     %%outside the UTM limits
	end; end; end; end; end; end; end; end; end; end;
	end; end; end; end; end; end; end; end; end; end;
return;
