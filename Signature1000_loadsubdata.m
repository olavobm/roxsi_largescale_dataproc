function list_files = Signature1000_filesintimelims(sigSN, time_lims)
%% list_files = SIGNATURE1000_FILESINTIMELIMS(dirfiles, time_lims)
%
%   inputs
%       - dirfiles:
%       - time_lims: in datenum.
%
%   outputs
%       - list_files:
%
%
%
%
%
% Loads Signature1000_time_limits.mat, created by ????
%
%


%%

%
sigTimeSplit = load('Signature1000_time_limits.mat');
sigTimeSplit = sigTimeSplit.sig1000timelims;


%%

%
list_fields = fieldnames(sigTimeSplit);

%
lmatchSN = contains(list_fields, sigSN);

%
if ~any(lmatchSN)
    error('SN not valid!')
end


%% Apply clock drift correction to time limits


% % time_lims.(list_fields{lmatchSN}).timelims


%%

%
lafter_begin = (sigTimeSplit.(list_fields{lmatchSN}).timelims.time_begin >= time_lims(1));
lbefore_end = (sigTimeSplit.(list_fields{lmatchSN}).timelims.time_begin <= time_lims(2));

%
ind_first = find(lafter_begin, 1, 'first');
ind_last = find(lbefore_end, 1, 'last');


%%


list_files = sigTimeSplit.(list_fields{lmatchSN}).timelims.filenames(ind_first:ind_last);




