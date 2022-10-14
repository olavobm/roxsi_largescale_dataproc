function [ux, uy, angle_rot] = ROXSI_uv_ENtoXY(ue, vn, site, lmagdec, roxsigrid)
%%  [ux, uy, angle_rot] = ROXSI_uv_ENtoXY(ue, vn, site, magdec, roxsigrid)
%
%   inputs
%       - ue: horizontal velocity component along east.
%       - vn: horizontal velocity component along north.
%       - site: string indicating ROXSI site (must match ROXSI_xygrids.mat).
%       - lmagdec (optional): set to false (default) to convert from
%                             geographic north. True to convert from
%                             magnetic coordinate system.
%       - roxsigrid (optional): structure from ROXSI_xygrids.mat.
%
%   outputs
%       - ux: velocity along local X axis.
%       - uy: velocity along local Y axis.
%       - angle_rot: angle of rotation (in degrees).
%
%
% ROXSI_UV_ENTOXY.m is a high-level function that rotates horizontal
% velocity components from ENU (ue and vn) to components (ux, uv)
% in the local coordinate system of a ROXSI array. These local coordinate
% systems (one for each ROXSI array) are defined by the structure in the
% file ROXSI_xygrids.mat.
%
% The input site (a string or character array) must match the fields in
% the structure saved in ROXSI_xygrids.mat. This file is loaded in the
% function, or alternatively you can pass the structure variable as input.
% (the function loads the mat from the current working directory).
%
% For lmagdec==false (default), no correction for magnetic declination is
% applied, i.e. the rotation is from the GEOGRAPHIC ENU to the local XY
% coordinate system. Set lmagdec to true to rotate from the MAGNETIC
% ENU to the local XY. As a reminder, magnetic declination is the angle
% of the magnetic north from the geographic north (positive in the
% CLOCKWISE direction).
% 
% The angle of rotation applied can be output as angle_rot. This is the
% angle (in degrees) of the rotation of the coordinate system, where
% positive in the COUNTERCLOCKWISE sense (a.k.a. trigonometric convention).


%% Check inputs

%
if ~isequal(size(ue), size(vn))
    error('Input velocity components do not have the same dimensions.')
end

%
if ~exist('lmagdec', 'var')
    lmagdec = false;
end

%
if (~iscolumn(ue)) && size(ue, 1) > nbinsTH
    larrayopposite = true;
    %
    ue = ue.';
    vn = vn.';
else
    larrayopposite = false;
end


%%
%
if ~exist('roxsigrid', 'var')
    %
    roxsigrid = load('ROXSI_xygrids.mat');
    %
    fnameaux = fieldnames(roxsigrid);
    roxsigrid = roxsigrid.(fnameaux{1});
end

%
if lmagdec
    magdec_angle = 12.86;
else
    magdec_angle = 0;
end


%%

%
angle_rot = 270 - roxsigrid.(site).angleref + magdec_angle;

%
matrix_rot = [ cosd(angle_rot), sind(angle_rot); ...
              -sind(angle_rot), cosd(angle_rot) ];


%%

%
if iscolumn(ue)
    ue = ue.';
    vn = vn.';
    lcolumn = true;
else
    lcolumn = false;
end


%%

%
ux = NaN(size(ue));
uy = ux;

%
for i = 1:size(ue, 1)

    %
    uv_rot_aux = matrix_rot * [ue(i, :); vn(i, :)];

    %
    ux(i, :) = uv_rot_aux(1, :);
    uy(i, :) = uv_rot_aux(2, :);

end

%%
%
if lcolumn
    ux = ux.';
    uy = uy.';
end

%
if larrayopposite
    ux = ux.';
    uy = uy.';
end


    