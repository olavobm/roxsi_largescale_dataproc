function [lkeepdata, Nclockinvs, timelims_clockinvs] = correct_Sig1000_clock(dtime)
%% [lkeepdata, Nclockinvs, timelims_clockinvs] = CORRECT_SIG1000_CLOCK(dtime)
% 
%   inputs
%       -
% 
%   outputs
%       -
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

%
Nclockinvs = length(ind_clock_backwards);


%%

%
lkeepdata = true(size(dtime));

%
timelims_clockinvs = NaT(Nclockinvs, 2);


%%

%
for i = 1:Nclockinvs

    % The index of the (supposedly) last bad timestamp
    % that is repeated
    indlastbadpt_aux = ind_clock_backwards(i);

    %
    lkeepdata(indlastbadpt_aux) = false;
    
% %     %
% %     indbadpt_aux = (ind_clock_backwards(i) - npts_TH) : (ind_clock_backwards(i) + npts_TH);

    %
    ind_match_aux = find(dtime((indlastbadpt_aux+1) : ...
                               (indlastbadpt_aux + 1 + npts_TH)) == ...
                         dtime(indlastbadpt_aux));

    % The index (from the beginning the array dtime) of the timestamp
    % that is the same as the repeated one. But this is the one that
    % will be kept
    indmatch_lastbad_aux = indlastbadpt_aux + ind_match_aux;
    
    %
    indloop = 1;
    
    %
    while indloop < npts_TH
        
        %
        ind_nextbad_aux = indlastbadpt_aux - indloop;
        ind_nextmatch_aux = indmatch_lastbad_aux - indloop;
        
        % If there is a match
        if (dtime(ind_nextbad_aux) == dtime(ind_nextmatch_aux))
            
            %
            lkeepdata(ind_nextbad_aux) = false;
            
            %
            indloop = indloop + 1;
            
        % 
        else
            % When there is no match, all bad points
            % have been found and the while loop should
            % be terminated.
            %
            % Add points to terminate the while loop
            indloop = indloop + npts_TH;
            
        end
    end
    
    % These are the times of the good timestamps bracketing
    % the bad timestamp segment that will be removed
    timelims_clockinvs(i, 1) = dtime(ind_nextbad_aux);
    timelims_clockinvs(i, 2) = dtime(indlastbadpt_aux + 1);
    
    
    keyboard
end