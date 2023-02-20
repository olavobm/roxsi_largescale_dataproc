function dirspec = dirspectra_xyz(dtsegment, cutoff, dtimegrid, dtwindow, dtimedata, xdata, ydata, zdata, bottomdepth, dspecmethod)
%% dirspec = dirspectra_XYZ(dtsegment, cutoff, dtimegrid, dtwindow, dtimedata, xdata, ydata, zdata, bottomdepth, dspecmethod)
%
%   inputs
%       - dtsegment:
%       - cutoff:
%       - dtimegrid:
%       - dtstats:
%       - dtimedata:
%       - xdata:
%       - ydata:
%       - zdata:
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
dirspec.dtime = dtimegrid;
dirspec.dtwindow = dtwindow;

%
dirspec.df = df;
dirspec.frequency = f;

%
dirspec.cutoff = cutoff;


%% Pre-allocate variables

%
S_f_theta_temp = NaN(length(f), Nangbins, length(dirspec.dtime), length(dspecmethod));
D_f_theta_temp = S_f_theta_temp;

%
S_f_temp = NaN(length(f), length(dirspec.dtime));
depth_avg = NaN(1, length(dirspec.dtime));


%%

%
for i1 = 1:length(dirspec.dtime)


    %% Select data
    
    %
    ind_match_time = find(dtimedata == dirspec.dtime(i1));
    
    %
    if isempty(ind_match_time)
        continue
    end
    
    %
    inds_getdata_aux = (-Nptswindow/2 : 1 : Nptswindow/2) + ind_match_time;
    inds_getdata_aux = inds_getdata_aux(1:end-1);    % such that we have the correct number of points
    
    %
    if any(inds_getdata_aux < 1) || any(inds_getdata_aux > length(xdata))
        continue
    end
    
    % Get data
    x = xdata(inds_getdata_aux);
    y = ydata(inds_getdata_aux);
    z = zdata(inds_getdata_aux);
    
    % ---------------------
    if any(isnan(x))
        %
        warning('NaN found in data. Skipping the corresponding timestamp for directional spectrum.')
        %
        continue
    end
    % ---------------------
    
    %
    x = x(:);
    y = y(:);
    z = z(:);
    
    %
    depth_avg(i1) = mean(bottomdepth(inds_getdata_aux));
    h = depth_avg(i1);    % h should be a POSITIVE number 

    
    %%
    
    % the z coordinate for pos is set by the 'igam' option sent to
    % dat2dspec, igam = 1 is z=0 at surface, positive z up
    pos = [0,  0,  0;
           0,  0,  0;
           0,  0,  0;
           1, 16, 17;
           1,  0,  0].';
       
    %
    Data = [t(:), z(:), x(:), y(:)];
    

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
            [Sd, D, Sw] = dat2dspec_frequency_cut(Data, pos, h, ...
                                                  nfft, Nangbins, ...
                                                  dspecmethod(i2), options, cutoff);

        else
            [Sd, D] = dat2dspec_frequency_cut(Data, pos, h, ...
                                              nfft, Nangbins, ...
                                              dspecmethod(i2), options, cutoff);
        end
        
        %
        S_f_theta_temp(:, :, i1, i2) = Sd.S.' .*(2*pi);   % dimension: direction x frequency
        D_f_theta_temp(:, :, i1, i2) = D.S.';
        
    end
    
    %       
    S_f_temp(:, i1) = Sw.S.*(2*pi);    
    
end


%% Add variables to output structure

%
dirspec.See = S_f_temp;
dirspec.bottomdepth = depth_avg(:);


%%

% Remove repeated direction bin (i.e. -pi and +pi)
S_f_theta_temp = S_f_theta_temp(:, 1:end-1, :, :);
D_f_theta_temp = D_f_theta_temp(:, 1:end-1, :, :);

% Remove repetition and go from radians to degrees
dirspec.direction = (180/pi) * Sd.theta(1:end-1);


%% Add variables to output structure

%
dirspec.list_methods_dirspec = dspecmethod;

%
for i = 1 : length(dspecmethod)
	%
    dirspec.(dspecmethod(i)).See = S_f_theta_temp(:, :, :, i);
    dirspec.(dspecmethod(i)).D = D_f_theta_temp(:, :, :, i);
    
    % Average in time
    dirspec.(dspecmethod(i)).See_avg = mean(dirspec.(dspecmethod(i)).See, 3, 'omitnan');
    
end





