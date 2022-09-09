%% Get Spotter data for directional spectra calculations

clear
close all


%%
% ---------------------------------------------
% ---------------------------------------------
% ---------------------------------------------

%%

addpath('/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/code_proc/')

%% Load NOAA tide

%
noaatides = load(['/Volumes/GoogleDrive/Shared drives/ROXSI' ...
                  '/LargeScale_Data_2022/RAW/noaa_mry_tides/tides_NOAA_MRY.mat']);
noaatides = noaatides.noaaTides;


%% Spotter data directory

%
dir_data = ['/Volumes/GoogleDrive/Shared drives/ROXSI' ...
            '/LargeScale_Data_2022/Level1_Data/Spotter_Level1/'];

% Or
% % load('/Volumes/LaCie/RAW/Spotters/SDcards/B01_spot1150.mat')


%% Take B03s data

%
list_files = {'B01_spot1150.mat', 'B01_spot1158.mat', 'B03_spot1152.mat', ...
              'B05_spot1153.mat', 'X01_spot1151.mat', 'X03_spot1157.mat', ...
              'X04_spot1155.mat'};

%%

%
spotterData = load(fullfile(dir_data, 'B03_spot1152.mat'));
spotterData = spotterData.s;

%
spotterData.bulkparameters.time.TimeZone = 'America/Los_Angeles';
spotterData.location.time.TimeZone = 'America/Los_Angeles';
spotterData.displacement.time.TimeZone = 'America/Los_Angeles';


%% Get 5 days of data to test the calculations

%
time_lims = [datetime(2022, 07, 11, 0, 0, 0), ...
             datetime(2022, 07, 16, 0, 0, 0)];
time_lims.TimeZone = 'America/Los_Angeles';


% Plot significant wave height and the period that I have chosen
figure
    plot(spotterData.bulkparameters.time, ...
         spotterData.bulkparameters.("Significant Wave Height"), '.-')

    %
    overlayline('v', time_lims, '-k')


%% Get location/coordinates

%
linlims_location = (spotterData.location.time >= time_lims(1)) & ...
                   (spotterData.location.time <= time_lims(2));

%
data4dirspec.location.dtime = spotterData.location.time(linlims_location);
data4dirspec.location.latitude = spotterData.location.("latitude (decimal degrees)")(linlims_location);
data4dirspec.location.longitude = spotterData.location.("longitude (decimal degrees)")(linlims_location); 


%% Get displacement

%
linlims_displacement = (spotterData.displacement.time >= time_lims(1)) & ...
                       (spotterData.displacement.time <= time_lims(2));

%
data4dirspec.displacement.dtime = spotterData.displacement.time(linlims_displacement);
data4dirspec.displacement.x = spotterData.displacement.("x (m)")(linlims_displacement);
data4dirspec.displacement.y = spotterData.displacement.("y(m)")(linlims_displacement);
data4dirspec.displacement.z = spotterData.displacement.("z(m)")(linlims_displacement);



%% Get and idealized bathymetry timeseries (just a placeholder
% constant value + tides)

%
zmean_ref = 22;

%
data4dirspec.location.waterdepth = zmean_ref + ...
                                          interp1(noaatides.dtime, ...
                                                  noaatides.zMSL, ...
                                                  data4dirspec.location.dtime);

%
data4dirspec.location.z_bathy = -data4dirspec.location.waterdepth;


%% Put in the format that the rest of the code wants

%
spotloc.dtime = data4dirspec.location.dtime;
spotloc.lat = data4dirspec.location.latitude;
spotloc.lon = data4dirspec.location.longitude;
spotloc.depth = data4dirspec.location.waterdepth;

%
spotdisp.dtime = data4dirspec.displacement.dtime;
spotdisp.x = data4dirspec.displacement.x;
spotdisp.y = data4dirspec.displacement.y;
spotdisp.z = data4dirspec.displacement.z;


%%


save("demo_displacement_deployed.mat", "spotdisp")
save("demo_location_deployed.mat", "spotloc")

