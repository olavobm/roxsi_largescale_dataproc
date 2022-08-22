function [latcenter, loncenter, radius] = Spotter_fitwatchcircle(latitude, longitude)
%% [latcenter, loncenter, radius] = SPOTTER_FITWATCHCIRCLE(latitude, longitude)
%
%   inputs
%       - latitude: vector of latitudes.
%       - longitude: longitudes (same dimensions as latitude).
%
%   outputs
%       - latcenter: latitude of the center of the computed watch circle.
%       - loncenter: corresponding longitude.
%       - radius: radius of the computed circle.
%
%
% SPOTTER_FITWATCHCIRCLE.m computes circle, least-squares.
%
% Since the calculation has to be done on a cartesian grid, geographical
% coordinates are first transformed into UTM coordinates, and then
% back to (latitude, longitude) at the end. These conversions are
% computed by separate functions.
%
% For a great reference on the (simple) linear algebra
% of fitting a circle to (x, y) data, see
% https://lucidar.me/en/mathematics/least-squares-fitting-of-circle/
%
%
% Olavo Badaro Marques.


%% Convert latitude/longitude to UTM coordinates
% to do calculations on a cartesian grid

[UTMEasting, UTMNorthing, UTMzone] = lltoUTM(latitude, longitude);


%% Remove the mean to make the linear algebra work

%
UTMEasting_mean = mean(UTMEasting);
UTMNorthing_mean = mean(UTMNorthing);

%
x_data = UTMEasting - UTMEasting_mean;
y_data = UTMNorthing - UTMNorthing_mean;


%% Create matrix and solve for least-squares coefficients

%
A_matrix = [x_data, y_data, ones(length(y_data), 1)];
B_matrix = x_data.^2 + y_data.^2;

%
X_vector = pinv(A_matrix)*B_matrix;


%% Calculate circle parameter

%
xc = X_vector(1)/2;
yc = X_vector(2)/2;
r = sqrt(4*X_vector(3) + X_vector(1)*X_vector(1) + X_vector(2)*X_vector(2))/2;


%% Convert center from UTM anomaly to geographical coordinate

%
eastingcenter = xc + UTMEasting_mean;
northingcenter = yc + UTMNorthing_mean;

%
[latcenter, loncenter] = UTMtoll(northingcenter, eastingcenter, UTMzone);

%
radius = r;

