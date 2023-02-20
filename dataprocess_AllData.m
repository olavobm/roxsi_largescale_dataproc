%% Script that runs all of the processing scripts for all
% instruments that Olavo has wrote code for (or you can
% run just specific blocks to process data from only certain instruments)

clear
close all


%% 

%
run('create_data_emptydirtree.m')

% If a directory tree already exists, it will be
% renamed and a new one will be created.


%% Process Spotters (note that the above is currently
% required for L2 processing of smart mooring data)

%
run('Spotters_proc_lvl_1.m'), close all
run('bathymetry_around_Spotters.m'), close all
% % run('Spotters_proc_lvl_2.m'), close all

% Merge Spotter data from B01 and Smart Mooring E07 and E09.
run('Spotter_merge_L1data.m')


%%
return

%% Aquadopps

%
run('Aquadopp_proc_lvl_1.m'), close all
% % run('Aquadopp_proc_lvl_2.m'), close all


%% Process Smart Moorings

%
run('SpottersSmart_proc_lvl_1.m'), close all
run('SpottersSmart_timegrid_pressure.m'), close all
% % run('SmartMooring_proc_lvl_2.m'), close all

% Merge Smart Mooring data
run('SmartMooring_merge_L1data.m')


%% Signature1000

%
run('Signature1000_proc_lvl_1_scalars.m'), close all
run('Signature1000_proc_lvl_1.m'), close all
% % run('Signature1000_proc_lvl_2.m'), close all
