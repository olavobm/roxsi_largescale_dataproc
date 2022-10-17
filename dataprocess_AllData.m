%% Script that runs all of the processing scripts for all
% instruments that Olavo has wrote code for (or you can
% run just specific blocks to process data from only certain instruments)

clear
close all



%% Process  Smart Moorings

%
run('SpottersSmart_proc_lvl_1.m')
%
run('SpottersSmart_timegrid_pressure.m')
%
run('SpottersSmart_proc_lvl_2.m')



%% Process Spotters (note that the above is currently
% required for L2 processing of smart mooring data)

%
run('Spotters_proc_lvl_1.m')
run('Spotters_proc_lvl_2.m')


%% Aquadopps



%% Signature1000