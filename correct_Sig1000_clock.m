function [dtime_out, dtimegaps] = correct_Sig1000_clock(dtime)
%% [dtime_out, dtimegaps] = CORRECT_SIG1000_CLOCK(dtime)
% 
%   inputs
%       -
% 
%   outputs
%       -
%       -
% 
% 
% CORRECT_SIG1000_CLOCK.m corrects the clock issues found on the
% Signature1000 SN 100231 (deployed at X05at). Looking at the
% diff(time) shows several clock inversions (a total of ???
% identified after running this function). These inversions are
% also accompanied by clock skipping forward (order of 2 seconds).
% Looking at a few examples, the problem is that the ADCP first
% stores a few data points in the future, then it goes back to a
% few data points and then it repeats the data points after the
% forward skip. This clock issue leads to a typical gap of ~0.8
% seconds.
% 
% The correction is to simply eliminate the first instance
% of the few data points and eliminate the clock inversion.
% The output of this function can then be used to eliminate
% the bad data points in the dependent variables.
% 
% 
% 


%% Threshold of number of points to scan around the clock
% inversion to find the points that need to be deleted.
% From a couple of examples, I expect that the bad/repeated
% points will come before the clock goes backwards, but just
% in case, I will write the code to look for those bad points
% around the clock inversion

% 
npts_TH = 30;


%%

%
ind_clock_backwards = find(diff(dtime) < 0);


%%

for i = 1:length(ind_clock_backwards)

    %
    indbadpt_aux = ind_clock_backwards(i);

% %     %
% %     indbadpt_aux = (ind_clock_backwards(i) - npts_TH) : (ind_clock_backwards(i) + npts_TH);

    %
    ind_match_aux = find(dtime((indbadpt_aux+1):(indbadpt_aux + 1 + npts_TH)) == dtime(indbadpt_aux));

    keyboard
end