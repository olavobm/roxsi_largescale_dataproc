function [Nsegs, indsegments] = find_continuous_segments(x, gapTH)
%% [Nsegs, indsegments] = FIND_CONTINUOUS_SEGMENTS(x, gapTH)
%
%   inputs
%       -
%       -
%
%   outputs
%       -
%       -
%
%
%
%
%


%%

%
x_diff = diff(x);

%
Npts = length(x);


%%

% Find gaps
indaboveTH = find(x_diff > gapTH);

%
Nsegs = length(indaboveTH) + 1;


%%

% --------------
%
if Nsegs == 1
    indsegments = [1, Npts];
    
    return
end

% --------------
%
indsegments = NaN(Nsegs, 2);
indsegments(1, 1) = 1;
indsegments(end, 2) = Npts;

%
for i = 1:(Nsegs-1)
    indsegments(i, 2) = indaboveTH(i);
    indsegments(i+1, 1) = indsegments(i, 2) + 1;
end



