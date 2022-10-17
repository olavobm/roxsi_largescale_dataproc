%% Script that runs all of the processing scripts for all
% instruments that Olavo has wrote code for (or you can
% run just specific blocks to process data from only certain instruments)

clear
close all



%% Process Smart Moorings

%
run('Spotters_proc_lvl_1.m')


%% Process Spotters

%
run('SpottersSmart_proc_lvl_1.m')
run('SpottersSmart_timegrid_pressure.m')