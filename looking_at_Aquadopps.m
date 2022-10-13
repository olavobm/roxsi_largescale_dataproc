%% Looking at aquadopp and low-backscatter
% and coming up with L2 data

clear
close all


%%

%

%
dir_AQDP_L1 = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/Aquadopp_Level1/';

%
dataE04 = load(fullfile(dir_AQDP_L1, 'roxsi_aquadopp_L1_13172_E04.mat'));
dataA03 = load(fullfile(dir_AQDP_L1, 'roxsi_aquadopp_L1_5380_A03.mat'));

%
dataE04 = dataE04.aquadoppL1;
dataA03 = dataA03.aquadoppL1;

%%

dataE04_L2 = dataE04;
dataA03_L2 = dataA03;


%% Compute depth-averaged currents

%
dataE04_L2.Uezavg = mean(dataE04_L2.Ue, 1, 'omitnan');
dataE04_L2.Vnzavg = mean(dataE04_L2.Vn, 1, 'omitnan');
dataE04_L2.Wupzavg = mean(dataE04_L2.Wup, 1, 'omitnan');
%
dataA03_L2.Uezavg = mean(dataA03_L2.Ue, 1, 'omitnan');
dataA03_L2.Vnzavg = mean(dataA03_L2.Vn, 1, 'omitnan');
dataA03_L2.Wupzavg = mean(dataA03_L2.Wup, 1, 'omitnan');

%
dataE04_L2.averaged.Uezavg = mean(dataE04_L2.averaged.Ue, 1, 'omitnan');
dataE04_L2.averaged.Vnzavg = mean(dataE04_L2.averaged.Vn, 1, 'omitnan');
dataE04_L2.averaged.Wupzavg = mean(dataE04_L2.averaged.Wup, 1, 'omitnan');
%
dataA03_L2.averaged.Uezavg = mean(dataA03_L2.averaged.Ue, 1, 'omitnan');
dataA03_L2.averaged.Vnzavg = mean(dataA03_L2.averaged.Vn, 1, 'omitnan');
dataA03_L2.averaged.Wupzavg = mean(dataA03_L2.averaged.Wup, 1, 'omitnan');



%%