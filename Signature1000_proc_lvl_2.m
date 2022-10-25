%% Signature1000_proc_lvl_2
%
% Since:
%   - This L2 processing for Signatures is basically
%     the same as for Aquadopps and;
%   - there is a lot of Signature data that may require
%     doing analyses on different files;
% I will write this script and a higher level function
% to do the L2 processing on Signature (or maybe any ADCP) data.
%
% If I take depth-averaged velocity, or in the first bin,
% I can probably patch all of the timeseries together
% and then do the analysis.


clear
close all


%%
% --------------------------------------------
% ----------------- PREAMBLE -----------------
% --------------------------------------------

%% Directory of L1 data

%
dirparent_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/';
dir_dataL1 = fullfile(dirparent_data, 'Level1_Data', 'Signature_Level1');


%% Directory to output L2 data and figures

%
dir_output_parent = dirparent_data;
dir_dataL2 = fullfile(dir_output_parent, 'Level2_Data', 'Signature_Level2');


%%
% --------------------------------------------
% ---------- PROCESSING PARAMETERS -----------
% --------------------------------------------

%%

% All Signatures
list_Signature = {'B10_103045', ...    % Do B first, because they have less
                  'B13_103046', ...    % data than A01, so if something is
                  'B15_103056', ...    % wrong with the code, the error
                  'B17_101923', ...    % pops up sooner
                  'A01_103043', ...
                  'C01_102128', ...
                  'X05_100231', ...
                  'X11_101941'};

Nsignatures = length(list_Signature);


%%
% --------------------------------------------
% ------------- DEFINE TIME GRID -------------
% --------------------------------------------

%% Define this instead of using the data
% (it's simpler this way)

% All ADCPs programmed to start at 2022/06/21 18:00:00
% (B10 not in the water at start time)
time_lims_L2 = datetime(2022, 06, 21, 19, 00, 00) : hours(1) : ...
               datetime(2022, 07, 25, 10, 00, 00);    % most ADCPs recovered by the 21st at 9:00AM.
time_lims_L2.TimeZone = 'America/Los_Angeles';


%%
% --------------------------------------------
% --------------------------------------------
% --------------------------------------------

%%

%
for i1 = 1:Nsignatures

    %%

    %% Load data and get reduced data
    
    %
    list_dirsegments = dir(fullfile(dir_dataL1, 'segment_*'));
    
    % Loop over segments
    for i2 = 1:length(list_dirsegments)
        %
        data_aux = load(fullfile(list_dirsegments(i2).folder, ...
                                 list_dirsegments(i2).name, ...
                                 ['roxsi_signature_L1_' list_Signature{i2} '.mat']));

        %
        if i2==1
            %
            sigL2.SN = sigL1.SN;
            sigL2.mooringID = sigL1.mooringID;
            % ... others ...
            %
            sigL2.zhab = sigL1.zhab;
        end

        %
        sigL2.dtime = sigL1.dtime;
        %
        sigL2.pressure = sigL1.pressure;
        %
        sigL2.udepthavg = mean(sigL1.u, 1, 'omitnan');
        sigL2.vdepthavg = mean(sigL1.v, 1, 'omitnan');
        sigL2.wdepthavg = mean(sigL1.w, 1, 'omitnan');
        %
        sigL2.ubin1 = sigL1.u(1, :);
        sigL2.vbin1 = sigL1.v(1, :);
        sigL2.wbin1 = sigL1.w(1, :);

        %
        sigL2.udepthavg = sigL2.udepthavg(:);
        sigL2.vdepthavg = sigL2.vdepthavg(:);
        sigL2.wdepthavg = sigL2.wdepthavg(:);
        %
        sigL2.ubin1 = sigL2.ubin1(:);
        sigL2.vbin1 = sigL2.vbin1(:);
        sigL2.wbin1 = sigL2.wbin1(:);
    end
    
    % Compute depth-averaged current
    
    
    %% Compute hourly depth-averaged current
    
    
    %% Compute spectra
    
    
    
    %% Compute See spectra from pressure (and vel??)

end