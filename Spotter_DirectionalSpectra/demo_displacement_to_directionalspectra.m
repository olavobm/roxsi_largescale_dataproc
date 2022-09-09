% demo_Process_displacement_to_directional_spectra_wafo.m
% adapted by Olavo on 09 Sep 2022
% 
% 
% Patrick Collins
% 25 Apr 22
% TITLE: process_displacement_to_directional_spectra_wafo
%
% PURPOSE: take displacement data and calculate directional spectra using
% the wafo scripts
%
% DATA REQUIREMENTS:
% - displacement data as a .mat data file
% - location data covering the displacement domain
% - the location data must have been processes to include the depth data
%
% NOTES: the displacement data will be trimmed to start and end on the
% hour, without regard to the analysis period.

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

%
displacement_data_path = pwd;
location_data_path = pwd;
%
workingpath = pwd;


%%

addpath(genpath(fullfile(workingpath, 'wafo')))


%%
%
dt = 0.4; %the sampling interval in seconds
nfft = 256; %the number of points in each window (nfft)
dtheta = 1; %the angle bin width for the spectral analysis
analysis_period_hours = 0.5; %the analysis period in hours. 


%%
% % dspec_method = ["EMEM","IMLM","MLM","MEM"]; % the directional algorithms to use

% dspec_method = ["EMEM","IMLM","MLM"]; % the directional algorithms to use

dspec_method = "IMLM";

% CHECK that the index for the displacement and location files are
% consistent (line 64)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The directional algorithms available are listed in dat2spec.m 
% 'BDM'  Bayesian Directional Spectrum Estimation Method 
% 'MLM'  Maximum Likelihood Method (default)
% 'IMLM' Iterative  Maximum Likelihood Method 
% 'MEM'  Maximum Entropy Method   (slow)
% 'EMEM' Extended Maximum Entropy Method

%%

% index the displacement and depth files
% The location data must have been processed to include depth 
% % cd(displacement_data_path)
disp = dir(fullfile(displacement_data_path, "*displacement_deployed.mat"));
loc = dir(fullfile(displacement_data_path, "*location_deployed.mat"));


%%

%
analysis_period_seconds = analysis_period_hours * 3600; % hourly estimates are being generated
N = analysis_period_seconds/dt;
t = 0:dt:(N-1)*dt; % a time vector that will be passed to wafo 
t = t(:);

%
Nt = 360/dtheta + 1; % the number of anglular bins, dtheta makes more sense

%
T = nfft.*dt; % the window length in seconds
df = 1/T; % frequency resolution, a consequence of T
fN = 0.5/dt; % nyquist frequency
f = 0:df:fN; f = f(:);

%
pos = [0,  0,  0; ...
       0,  0,  0; ...
       0,  0,  0; ...
       1, 16, 17; ...
       1,  0,  0].';    % data to the wafo, look at the
                        % dat2dspec.m and tran.m for more details 

%%
% now loop through each displacement file to generate directional spectra

for ii = 1: length(disp)
    
    load(disp(ii).name);
    %%%%%%%%% This needs to correct for the displacement file name 
    site_name = disp(ii).name(1:end-26);
    load(loc(ii).name); % check that the index is the same for the displacement and location file... it should ne

    % pull out an analysis period of data
    % ind is the index of the start of the first full hour
    ind_start = find(minute(spotdisp.dtime) == 0, 1);

    dtime = spotdisp.dtime(ind_start) : hours(analysis_period_hours) : ...
        spotdisp.dtime(end) - hours(analysis_period_hours);

    analysis_periods = length(dtime); % this is the number of analysis period in the data

    % preallocate for speed
    S_f_theta_temp = nan(length(f), Nt, analysis_periods, length(dspec_method));
    D_f_theta_temp = nan(length(f), Nt, analysis_periods, length(dspec_method));
    S_f_temp = nan(length(f), analysis_periods);
    depth = nan(1, analysis_periods);
    lat = nan(1, analysis_periods);
    lon = nan(1, analysis_periods);

    % Print message to the screen
    %
    disp(' ')
    disp(' ')
    disp(['----- The total number of analysis ' ...
          'periods is: ' num2str(analysis_periods) ' -----'])
    
    %
    for sample = 1 : analysis_periods
        
        %
        disp(' ')
        disp(' ')

        %
        data_index = ind_start + (sample -1) *N : ...
                        ind_start + sample *N -1;
        xt = spotdisp.x(data_index);
        yt = spotdisp.y(data_index);
        zt = spotdisp.z(data_index);
        dtime_sample  = spotdisp.dtime(data_index); 
     

        %the index of the depth record with the same dtime as the
        %displacement 

        good = spotloc.dtime >= dtime_sample(1) & spotloc.dtime <= dtime_sample(end);
        lat(sample) = mean(spotloc.lat(good),'omitnan');
        lon(sample) = mean(spotloc.lon(good),'omitnan');
        depth(sample) = mean(spotloc.depth(good), 'omitnan');
        h = abs(depth(sample));

        Data = [t, zt, xt, yt];
        
        for jj = 1 : length(dspec_method)
            if jj == 1
                [Sd,D,Sw] = dat2dspec(Data,pos,h,nfft,Nt,dspec_method(jj));
            else
                [Sd,D] = dat2dspec(Data,pos,h,nfft,Nt,dspec_method(jj));
            end
            S_f_theta_temp(:,:,sample,jj) = Sd.S'.*(2*pi); %scaled because I want direction in degrees
            D_f_theta_temp(:,:,sample,jj) = D.S';
        end
%       
        S_f_temp(:,sample) = Sw.S.*(2*pi);    

        %
        disp(['----- Done with iteration/period ' num2str(sample) ' ' ...
              'out of ' num2str(analysis_periods) ' -----'])
       
    end

    % from the wafo documentation
%             theta  = angle vector -pi..pi of length Nt 
%                      (theta = 0 -> + x-axis, theta = pi/2 -> + y-axis) 
% so in this angle definition, theta = 0 mean waves heading towards the
% east, (from the west) 

    % change the angular definition from cartesian to nautical
    dir = rad2deg(Sd.theta(1:end));
    dir_naut_unsorted= mod(270 - rad2deg(Sd.theta(1:end-1)),360);
    [dir_naut, ind] = sort(dir_naut_unsorted);

    dspec.site = site_name;
%     save each dspec method to the dspec structure
    for jj = 1 : length(dspec_method)
        eval("dspec.S_f_theta_" + dspec_method(jj) + "= S_f_theta_temp(:,ind,:,jj);")
        eval("dspec.D_f_theta_" + dspec_method(jj) + "= D_f_theta_temp(:,ind,:,jj);")
    end
    dspec.S_f = S_f_temp;
    dspec.f = f;
    dspec.direction_nautical = dir_naut;
    dspec.dtime = dtime(:);
    dspec.depth = depth(:);
    dspec.lon   = lon(:);
    dspec.lat   = lat(:);
    dspec.nfft = nfft;
    dspec.analysis_period_hours = analysis_period_hours;

    fname = site_name + "_dspec.mat";

    save(fname, "dspec" ,'-v7.3')

end



