function [list_files, sigTimeSplit] = Signature1000_filesintimelims(sigSN, time_lims)
%% list_files = SIGNATURE1000_FILESINTIMELIMS(dirfiles, time_lims)
%
%   inputs
%       - dirfiles:
%       - time_lims: in datenum.
%
%   outputs
%       - list_files:
%       - sigTimeSplit: the structure variable
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
lfiles_edge_1_inlims = (time_lims(1) >= sigTimeSplit.(list_fields{lmatchSN}).timelims.time_begin) & ...
                       (time_lims(1) <  sigTimeSplit.(list_fields{lmatchSN}).timelims.time_end);
%
lfiles_edge_2_inlims = (time_lims(2) >= sigTimeSplit.(list_fields{lmatchSN}).timelims.time_begin) & ...
                       (time_lims(2) <  sigTimeSplit.(list_fields{lmatchSN}).timelims.time_end);
%
lfiles_fully_inlims = (time_lims(1) <= sigTimeSplit.(list_fields{lmatchSN}).timelims.time_begin) & ...
                      (time_lims(2) >  sigTimeSplit.(list_fields{lmatchSN}).timelims.time_end);

%
lfiles_withdata_inlims = lfiles_edge_1_inlims | lfiles_edge_2_inlims | lfiles_fully_inlims;

%
ind_first = find(lfiles_withdata_inlims, 1, 'first');
ind_last = find(lfiles_withdata_inlims, 1, 'last');


%%

list_files = sigTimeSplit.(list_fields{lmatchSN}).timelims.filenames(ind_first:ind_last);




