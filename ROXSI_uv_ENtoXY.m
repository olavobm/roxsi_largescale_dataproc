function [ux, uy] = ROXSI_uv_ENtoXY(ue, vn, site, roxsigrid)
%%  [ux, uy] = ROXSI_uv_ENtoXY(u, v, site, roxsigrid)
%
%
%
%
%
%
%
%
%
%
%
%

%%

%
if ~exist('roxsigrid', 'var')
    %
    roxsigrid = load('ROXSI_xygrids.mat');
    %
    fnameaux = fieldnames(roxsigrid);
    roxsigrid = roxsigrid.(fnameaux{1});
end


%%
if ~isequal(size(ue), size(vn))
    error('!!!!')
end


%%
%
angle_rot = 270 - roxsigrid.(site).angleref;

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

%
if lcolumn
    ux = ux.';
    uy = uy.';
end
    