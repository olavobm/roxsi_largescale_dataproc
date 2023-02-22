function adcpdirspec = dirspectra_ADCP_uvp(dtsegment, cutoff, dtimegrid, dtwindow, zhab, dtimedata, udata, vdata, pdata, bottomdepth, dspecmethod)
%% 
%
%   inputs
%       - dtsegment:
%       - cutoff:
%       - dtimegrid:
%       - dtstats:
%       - dtimedata:
%       - u:
%       - v:
%       - p:
%       - dspecmethod:
%
%   outputs
%       - adcpdirspectra
%
%
% DIRSPECTRA_ADCP_UVP.m is a high-level function that calls WAFO's
% toolbox to compute a directional spectrum from ADCP + pressure data.
%
% The direction vector in the output structure is in degrees,
% rotating CCW, and wrapping at -180/+180 (i.e. the convention
% of atan2). This direction coordinate is for where the waves
% are GOING TO.
%
% Olavo Badaro Marques


%%

%
if ~exist('dspecmethod', 'var')
    %
    % dspecmethod = ["EMEM","IMLM","MLM","MEM"];    % the directional algorithms to use
    dspecmethod = "EMEM";
end


%%

dtheta = 1;
Nangbins = 360/dtheta + 1; % the number of anglular bins, dtheta makes more sense


%% Get sampling period from the data time vector

% Should make sure it's gridded!!!!

%
dtsampling = seconds(dtimedata(2) - dtimedata(1));    % this should have double class!

%
Nptswindow = dtwindow/dtsampling;


%% Get number of points where fft will be computed over

%
nfft = (1/dtsampling) * dtsegment;    % probably needs to be even.

%
df = (1/dtsegment);

%
noverlap = nfft/2;

%
f = 0:df:cutoff;
f = f(:);


%% Parameters

%
% % dt = seconds(diff(dtstats(;
% % nfft = 512;    % 512 points is 128 seconds with the 4 Hz sampling of B13
% % 
% % dtheta = 1;
% % 
% % %
% % %dspec_method = ["EMEM","IMLM","MLM","MEM"]; % the directional algorithms to use
% % dspec_method = "EMEM"; % the directional algorithms to use
% % 
% % %
% % % cutoff = 0.20;    % the cut-off frequency to use because the transfer function blows up 
% % % cutoff = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = specoptset('dat2dspec', 'noverlap', noverlap, ...
                                  'window', hanning(nfft), 'igam', 1);
%%% the 'igam' parameter is how you define the z coordinate system, see
%%% specoptset for more details. 1 is z=0 at mslm with z increasing
%%% upwards, this becomes important when setting the 'pos' variable

%% More parameters

% 
% % analysis_period_seconds = analysis_period_hours * 3600; % hourly estimates are being generated
% N = dtsegment/dt;


% A time vector from 0 used by WAFO
t = 0:dtsampling:((Nptswindow-1)*dtsampling);
t = t(:);


%%

%
adcpdirspec.dtime = dtimegrid;
adcpdirspec.dtwindow = dtwindow;

%
adcpdirspec.df = df;
adcpdirspec.frequency = f;

%
adcpdirspec.cutoff = cutoff;


%% Pre-allocate variables

% % vel_range = adcp.vel_range;
% % % depth = adcp.depth; 

% % ind_start = find(minute(adcp.dtime) == 0, 1);
% % dtime = adcp.dtime(ind_start) : hours(analysis_period_hours) : ...
% %         adcp.dtime(end);
% % analysis_periods = length(dtime);


% preallocate for speed
S_f_theta_temp = NaN(length(f), Nangbins, length(adcpdirspec.dtime), length(dspecmethod));
D_f_theta_temp = S_f_theta_temp;

%
S_f_temp = NaN(length(f), length(adcpdirspec.dtime));
% % depth_avg = NaN(1, length(adcpdirspec.dtime));


%% Break the data apart in columns do the calculation more efficiently
% (this way doesn't allow for overlaps and the window is the same
% as the grid spacing)

%
[indsub, reshapeNdims, indgridlims] = reshapeintoWindows(dtimedata, dtimegrid);


% Get data
u_array = reshape(udata(indsub), reshapeNdims);
v_array = reshape(vdata(indsub), reshapeNdims);
%
p_array = reshape(pdata(indsub), reshapeNdims);

%
depth_avg = mean(reshape(bottomdepth(indsub), reshapeNdims), 1, 'omitnan');
h = depth_avg;     % h must be a POSITIVE 

%
lgoodcolumns = ~any(isnan(u_array), 1) & ~any(isnan(p_array), 1);

%
indfilloutput = indgridlims(1) : 1 : indgridlims(2);


%% Loop over the windows and compute directional spectra

% Loop over columns of the data arrays (i.e. the windows)
for i1 = 1:reshapeNdims(2)
    
    % Skip calculation if there is a NaN
    if ~lgoodcolumns(i1)
        continue
    end
    
    %%
    
    % the z coordinate for pos is set by the 'igam' option sent to
    % dat2dspec, igam = 1 is z=0 at surface, positive z up
    pos = [0 0 0;
           0 0 0;
           1.*[(-h(i1) + zhab(1)), (-h(i1) + zhab(1)), (-h(i1) + zhab(2))];
           10 11 9;
           0 0 1]';     % line 4, 10 = U, 11 = V, 9 = p 12 = W

    % %     t = 0:dt:(dt*(length(u)-1));    % the way this was defined, the
                                        % data is assumed to be trimmed /
                                        % or you can get different DOF

    %
    Data = [t(:), u_array(:, i1), v_array(:, i1), (1e4.*p_array(:, i1))];    % convert pressure to Pa from dbar 
   
    
    %% Call Pat's high-level function, that calls WAFO's function
    % that computes directional spectra
    
    %
    warning('off', 'WAFO:W2K')    % turn warning message off. It's a
                                  % warning because the wavenumber computed
                                  % from frequency doesn't converge. But
                                  % the result I get looks Ok, and it's
                                  % likely some silly problem with the WAFO
                                  % code
    
    % Loop over methods
    for i2 = 1 : length(dspecmethod)
        
        %
        if i2 == 1
            [Sd, D, Sw] = dat2dspec_frequency_cut(Data, pos, h(i1), ...
                                                  nfft, Nangbins, ...
                                                  dspecmethod(i2), options, cutoff);
        else
            [Sd, D] = dat2dspec_frequency_cut(Data, pos, h(i1), ...
                                              nfft, Nangbins, ...
                                              dspecmethod(i2), options, cutoff);
        end
        
        %
        S_f_theta_temp(:, :, indfilloutput(i1), i2) = Sd.S.' .*(2*pi);   % dimension: direction x frequency
        D_f_theta_temp(:, :, indfilloutput(i1), i2) = D.S.';
        
    end
    
    %       
    S_f_temp(:, indfilloutput(i1)) = Sw.S.*(2*pi);    
    
end

 

%% Get the data in a way that's not very efficient

% % %
% % for i1 = 1:length(adcpdirspec.dtime)
% % 
% % % %             % 
% % % %         linanalysis_disp = (spotterL1.displacement.dtime >= (dtime_proc_aux(sample) - (hours(analysis_period_hours)/2))) & ...
% % % %                            (spotterL1.displacement.dtime  < (dtime_proc_aux(sample) + (hours(analysis_period_hours)/2)));
% % % %         % This is slower than the index approach, but shouldn't
% % % %         % make much of a difference because the directional 
% % % %         % spectra takes a lot more time
% % % %         %
% % % %         xt = spotterL1.displacement.x(linanalysis_disp);
% % % %         yt = spotterL1.displacement.y(linanalysis_disp);
% % % %         zt = spotterL1.displacement.z(linanalysis_disp);
% % % %         dtime_sample = spotterL1.displacement.dtime(linanalysis_disp); 
% % % %      
% % % % 
% % % %         % the index of the depth record with the same dtime as the displacement 
% % % %         linanalysis_location = (spotterL1.location.dtime >= (dtime_proc_aux(sample) - (hours(analysis_period_hours)/2))) & ...
% % % %                                (spotterL1.location.dtime  < (dtime_proc_aux(sample) + (hours(analysis_period_hours)/2)));
% % % %         %
% % % %         lat(sample) = mean(spotterL1.location.latitude(linanalysis_location), 'omitnan');
% % % %         lon(sample) = mean(spotterL1.location.longitude(linanalysis_location), 'omitnan');
% % % %         depth(sample) = mean(spotterL1.location.z_msl(linanalysis_location), 'omitnan');
% % % %         %
% % % %         h = abs(depth(sample));
% %     
% %     %% Select data
% %     
% %     %
% %     ind_match_time = find(dtimedata == adcpdirspec.dtime(i1));
% %     
% %     %
% %     if isempty(ind_match_time)
% %         continue
% %     end
% %     
% %     %
% %     inds_getdata_aux = (-Nptswindow/2 : 1 : Nptswindow/2) + ind_match_time;
% %     inds_getdata_aux = inds_getdata_aux(1:end-1);    % such that we have the correct number of points
% %     
% %     %
% %     if any(inds_getdata_aux < 1) || any(inds_getdata_aux > length(udata))
% %         continue
% %     end
% % 
% % %     %
% % %     linsample_aux = dataADCP.dtime >= 
% % %     good = adcp.dtime>= dtime(sample) & ...
% % %            adcp.dtime < dtime(sample) + hours(analysis_period_hours);
% %       
% %     % Get data
% %     u = udata(inds_getdata_aux);
% %     v = vdata(inds_getdata_aux);
% %     %
% %     p = pdata(inds_getdata_aux);     % p = detrend(p, 2);
% %     
% %     % -----------
% %     if any(isnan(u))
% %         warning('NaN found in data')
% %         continue
% %     end
% %     % -----------
% %     
% %     %
% %     u = u(:);
% %     v = v(:);
% %     p = p(:);
% %     
% %     %
% %     depth_avg(i1) = mean(bottomdepth(inds_getdata_aux));
% %     h = depth_avg(i1);    % h should be a POSITIVE number 
% % 
% %     
% %     %%
% %     
% %     % the z coordinate for pos is set by the 'igam' option sent to
% %     % dat2dspec, igam = 1 is z=0 at surface, positive z up
% %     pos = [0 0 0;
% %            0 0 0;
% %            1.*[(-h + zhab(1)), (-h + zhab(1)), (-h + zhab(2))];
% %            10 11 9;
% %            0 0 1]';     % line 4, 10 = U, 11 = V, 9 = p 12 = W
% %     
% % % %     t = 0:dt:(dt*(length(u)-1));    % the way this was defined, the
% %                                         % data is assumed to be trimmed /
% %                                         % or you can get different DOF
% %     
% %     %
% %     Data = [t(:), u(:), v(:), (1e4.*p(:))];    % convert pressure to Pa from dbar 
% %     
% % 
% %     %% Call Pat's high-level function, that calls WAFO's function
% %     % that computes directional spectra
% %     
% %     %
% %     warning('off', 'WAFO:W2K')    % turn warning message off. It's a
% %                                   % warning because the wavenumber computed
% %                                   % from frequency doesn't converge. But
% %                                   % the result I get looks Ok, and it's
% %                                   % likely some silly problem with the WAFO
% %                                   % code
% %     
% %     % Loop over methods
% %     for i2 = 1 : length(dspecmethod)
% %         
% %         %
% %         if i2 == 1
% %             [Sd, D, Sw] = dat2dspec_frequency_cut(Data, pos, h, ...
% %                                                   nfft, Nangbins, ...
% %                                                   dspecmethod(i2), options, cutoff);
% %         else
% %             [Sd, D] = dat2dspec_frequency_cut(Data, pos, h, ...
% %                                               nfft, Nangbins, ...
% %                                               dspecmethod(i2), options, cutoff);
% %         end
% %         
% %         %
% %         S_f_theta_temp(:, :, i1, i2) = Sd.S.' .*(2*pi);   % dimension: direction x frequency
% %         D_f_theta_temp(:, :, i1, i2) = D.S.';
% %         
% %     end
% %     
% %     %       
% %     S_f_temp(:, i1) = Sw.S.*(2*pi);    
% %     
% % end


%% 

%
adcpdirspec.See = S_f_temp;
adcpdirspec.bottomdepth = depth_avg(:);

%%

% Remove repeated direction bin (i.e. -pi and +pi)
S_f_theta_temp = S_f_theta_temp(:, 1:end-1, :, :);
D_f_theta_temp = D_f_theta_temp(:, 1:end-1, :, :);

% Remove repetition and go from radians to degrees
adcpdirspec.direction = (180/pi) * Sd.theta(1:end-1);


%%

% % %
% % dir_naut_unsorted = mod(270 - rad2deg(Sd.theta(1:end-1)), 360);   % 1:end-1 to remove repeated direction bin
% % 
% % %
% % [dir_naut, ind] = sort(dir_naut_unsorted);
% % %
% % S_f_theta_temp = S_f_theta_temp(:, ind, :, :);
% % D_f_theta_temp = D_f_theta_temp(:, ind, :, :);


%% Add variables to output structure

%
adcpdirspec.list_methods_dirspec = dspecmethod;

%
for i = 1 : length(dspecmethod)
	%
    adcpdirspec.(dspecmethod(i)).See = S_f_theta_temp(:, :, :, i);
    adcpdirspec.(dspecmethod(i)).D = D_f_theta_temp(:, :, :, i);
    
    % Average in time
    adcpdirspec.(dspecmethod(i)).See_avg = mean(adcpdirspec.(dspecmethod(i)).See, 3, 'omitnan');
    
end





