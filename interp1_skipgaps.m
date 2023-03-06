function xinterp = interp1_skipgaps(t, x, gapTH, tgrid)
%% xinterp = INTERP1_SKIPGAPS(t, x, gapTH, tgrid)
%
%   inputs
%       - t: independent variable.
%       - x: dependent variable to interpolate
%       - gapTH: threshold of the gap
%       - tgrid:
%
%   outputs
%       - xinterp: interpolated variable.
%
%
%
%

% MAKE SURE IT WORKS IN UNDESIRED SITUATIONS

%%


% Break data apart in continuos segments (where sampling
% time never goes larger than threshold). Then interpolate
% data without interpolating over gaps

%
[Nsegs, indsegments] = find_continuous_segments(t, gapTH);


%% Assign time grid to data output structure
% and pre-allocate space for gridded variables

%
xinterp = NaN(length(tgrid), 1);

        
%% Now interpolate variables to time grid

% Interpolate over continuous segments
for i3 = 1:Nsegs
    %
    inds_data = indsegments(i3, 1) : 1 : indsegments(i3, 2);

    % Get grid points within the segment
    linseg_aux = (tgrid >= t(inds_data(1))) & (tgrid <= t(inds_data(end)));
    
    % Interpolate
    xinterp(linseg_aux) = interp1(t(inds_data), x(inds_data), tgrid(linseg_aux));

             
% %     % Interpolate over variables
% %     for i2 = 1:length(list_variables_aux)
% %         %
% %         if ~strcmp(list_variables_aux{i2}, 'dtime')
% %             %
% %             spotterL1.(list_subfields{i2}).(list_variables_aux{i2})(linseg_aux) = ...
% %                     interp1(dataSpotter.(list_subfields{i2}).dtime(inds_data_aux), ...
% %                             dataSpotter.(list_subfields{i2}).(list_variables_aux{i4})(inds_data_aux), ...
% %                             time_grid_aux(linseg_aux));
% % 
% %         end
% %     end

end