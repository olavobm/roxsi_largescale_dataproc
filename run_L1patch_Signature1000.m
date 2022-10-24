%% Script that (if necessary) patches Signatur1000 data

%
clear
close all


%%

%
dir_parent_data = ['/home/omarques/Documents/obm_ROXSI' ...
                   '/obm_DataLocal/Level1_Data/split_data/'];

%
list_Signature = {'X11_101941'};

%
diroutput = dir_parent_data;

%% Start diary


%% Patch data

Signature1000_patch_L1data(dir_parent_data, list_Signature, diroutput)


%% Close diary

% % toc(runPatchScript)


