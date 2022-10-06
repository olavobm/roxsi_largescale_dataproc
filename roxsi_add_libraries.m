function roxsi_add_libraries()
%% roxsi_add_libraries()
%
% ROXSI_ADD_LIBRARIES.m is a function to add other
% libraries that are used by the ROXSI data processing
% repository. This function has no inputs or outputs.


%%

%
dir_thisfunction = pwd();
dir_otherlibraries = fullfile(dir_thisfunction, '..');


%%

addpath(genpath(fullfile(dir_otherlibraries, 'cmocean')))
addpath(genpath(fullfile(dir_otherlibraries, 'ADCPtools')))



