function [lkeepdata, Nclockinvs, timelims_clockinvs] = correct_Sig1000_clock(dtime)
%% [lkeepdata, Nclockinvs, timelims_clockinvs] = CORRECT_SIG1000_CLOCK(dtime)
% 
%   inputs
%       - dtime: datetime vector.
% 
%   outputs
%       - lkeepdata: logical array, same size as dtime, where false
%                    is for data points in the clock inversion that
%                    need to be removed.
%       - Nclockinvs: number of clock inversions.
%       - timelims_clockinvs: an array of datetimes (Nclockinvs x 2),
%                             with the datetimes of the good timestamps
%                             bracketing each period of clock inversion.
%
% 
% CORRECT_SIG1000_CLOCK.m corrects the clock issues found on the
% Signature1000 SN 100231 (deployed at X05at). Looking at the
% diff(time) shows several clock inversions (a total of 69
% identified after running this function for the entire deployment
% in ROXSI 2022). These inversions are also accompanied by clock
% skipping forward (order of 2 seconds). Looking at a few examples,
% the problem is that the ADCP first stores a few data points in
% the future, then it goes back to a few data points and then it
% repeats the data points after the forward skip. This clock issue
% leads to a typical gap of ~0.8 seconds.
% 
% The correction is to simply eliminate the first instance
% of the few data points and eliminate the clock inversion.
% The output of this function can then be used to eliminate
% the bad data points in the dependent variables.
%
% Looking at the corrected data, there seems to be a another handful
% of clock issues where the clock does not reverse in time, it just
% goes forward and backward by a small amount about the sampling period.
%
% In addition to SN 100231 (X05at), there is another Signature1000
% (SN ???) with a handful of clock inversions.


%% Threshold of number of points to scan around the clock
% inversion to find the points that need to be deleted.

% 
npts_TH = 30;


%% Find the clock inversions and get the total number of them

%
ind_clock_backwards = find(diff(dtime) < 0);

%
Nclockinvs = length(ind_clock_backwards);


%% Initialize output variables

%
lkeepdata = true(size(dtime));

%
timelims_clockinvs = NaT(Nclockinvs, 2);
timelims_clockinvs.TimeZone = dtime.TimeZone;


%% Find the lcoation of the clock inversions

% Loop over clock inversions
for i = 1:Nclockinvs

    % The index of the last bad timestamp
    % that is repeated in the i'th instance
    % of clock inversion
    indlastbadpt_aux = ind_clock_backwards(i);

    % Add the location of this last bad data point
    lkeepdata(indlastbadpt_aux) = false;
    
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

end