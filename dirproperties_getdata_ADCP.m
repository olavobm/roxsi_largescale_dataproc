%% Get ADCP + pressure data for directional spectra and
% directional property calculations

clear
close all

%
addpath(genpath('/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc/'))


%%
% ----------------------------------------
% ----------------------------------------
% ----------------------------------------

%%

%
dirL1data = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/';

%
diroutputdata = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results_2/data_for_dirspectra/';


%%

% All Signatures
list_Signatures = {'B10', 'B13', 'B15', 'B17', 'A01', 'C01', 'X05', 'X11'};
% % list_Signatures = {};


%%

%
list_Aquadopps = {'B11', 'X13', ...    % 2 MHz
                  'B08', 'E03', 'E04', 'E06', 'E12', 'X06'};     % 1 MHz



% % list_Aquadopps = {'B08', '?????'};
% % list_Aquadopps = {'B08', 'B11'};



%% Time limits

% % % Small segment
% % dtimeload = [datetime(2022, 06, 27, 00, 00, 00), ...
% %              datetime(2022, 07, 01, 00, 00, 00)];

% Time limits for the whole experiment
dtimeload = [datetime(2022, 06, 21, 17, 30, 00), ...
             datetime(2022, 07, 25, 18, 00, 00)];

%
dtimeload.TimeZone = 'America/Los_Angeles';


%% List of metadata variables

list_metadatafields = {'SN', 'mooringID', ...
                       'latitude', 'longitude', 'site', 'X', 'Y', 'coordsystem', ...
                       'transducerHAB', 'binsize', 'zhab', 'samplingrateHz'};

                   
%% Define a threshold for interpolating over gaps.
% (note that there are slightly different definitions
% of gap -- e.g. time between good data points bracketing
% a gap; time between bad data points at the edges of a gap.)

gapTH = 5;    % in seconds


%%
% ----------------------------------------
% -------- GET SIGNATURE1000 DATA --------
% ----------------------------------------

%% Print progress message

%
disp(' '), disp(' ')
%
disp('----- Starting to get Signature1000 data for directional estimates. Loading data from:  -----')
for i = 1:length(list_Signatures)
    disp([num2str(i) ') ' list_Signatures{i}])
end


%% Now get the data

%
for i1 = 1:length(list_Signatures)

    %% Print message

    disp(['---- Getting data from mooring ' list_Signatures{i1} ' ----'])

    %% First load the scalars data to get the bottom depth

    %
    filescalar_aux = dir(fullfile(dirL1data, 'Signature_Level1', ...
                                  ['roxsi_signature_L1_' list_Signatures{i1} '_*_scalars.mat']));
    %
    sigL1scalars = load(fullfile(filescalar_aux.folder, filescalar_aux.name));
    sigL1scalars = sigL1scalars.sigL1;


    %% Copy relevant metadata

    %
    for i2 = 1:length(list_metadatafields)
        %
        dataADCP.(list_metadatafields{i2}) = sigL1scalars.(list_metadatafields{i2});
        %
        if strcmp(list_metadatafields{i2}, 'mooringID')
            dataADCP.instrument = "Signature1000";
        end
    end
    

    %% Select the ADCP cell that will be taken

    %
    min_bottomdepth = min(sigL1scalars.bottomdepthfrompres);

    %
    ind_binselect = find((sigL1scalars.zhab + (sigL1scalars.binsize/2)) <= min_bottomdepth, 1, 'last');

    % If ind_binselect is empty, then this should be a shallow
    % ADCP and the first bin will be taken instead
    if isempty(ind_binselect)
        ind_binselect = 1;
    end

    %
    dataADCP.zhab = dataADCP.zhab(1:ind_binselect);
    %
    dataADCP.binN = ind_binselect;

    %
    disp(['---- Selecting velocity bin number ' num2str(dataADCP.binN) ...
          ', where hab = ' num2str(dataADCP.zhab(end), '%.1f') ' m ----'])


    %% Get the data

    %
    sigL1 = Signature_load_L1fulldata(list_Signatures{i1}, dtimeload, dataADCP.binN);

    %
    dataADCP.dtime = sigL1.dtime;
    dataADCP.pressure = sigL1.pressure;
    %
    dataADCP.bottomdepthfrompres = sigL1.bottomdepthfrompres(:);
    %
    dataADCP.u = sigL1.u(:);
    dataADCP.v = sigL1.v(:);
    dataADCP.w = sigL1.w(:);


    %% Check for NaNs

    %
    inds_nan = find(isnan(dataADCP.u));
    Npts_velNaN = length(inds_nan);

    %
    disp(['Velocity data have ' num2str(Npts_velNaN) ' NaNs.'])
    %
    if Npts_velNaN>0
        warning('NaNs found in L1 data!')
    end

    
    %% Interpolate over short gaps (short as defined by
    % the threshold constant on top of this code)
    
    %
    figure, plot(isnan(dataADCP.u)), hold on, grid on
    
    
    %
    if Npts_velNaN>0

        %
        dt_samplingrate = seconds(diff(dataADCP.dtime(1:2)));
        
        % If there is a (normal) number of gaps
        if length(inds_nan)>1
            
            %
            ind_space_between_gaps = find(diff(inds_nan) > 1);

            % indboundgaps ASSUMES THAT THE TIMESERIES DOESN'T
            % START OR END WITH NANS!!!
            %
            indboundgaps = NaN(length(ind_space_between_gaps), 2);
            %
            indboundgaps(1, 1) = inds_nan(1);
            %
            for i2 = 1:(length(ind_space_between_gaps)-1)
                indboundgaps(i2, 2) = inds_nan(ind_space_between_gaps(i2));
                indboundgaps(i2+1, 1) = inds_nan(ind_space_between_gaps(i2) + 1);
            end
            %
            indboundgaps(end, 1) = inds_nan(ind_space_between_gaps(end) + 1);
            indboundgaps(end, 2) = inds_nan(end);
            
            
        % If there is only one NaN point (that has happened for B13)
        else
            
            indboundgaps = inds_nan.*[1, 1];
        end
        

        % Now loop over continuos gaps, and interpolate
        % if their length is equal or shorter than the threshold
        length_gaps_aux = (indboundgaps(:, 2) - indboundgaps(:, 1)) + 1;
        for i2 = 1:size(indboundgaps, 1)
            %
            if (dt_samplingrate*length_gaps_aux(i2)) <= gapTH
                %
                inds_goodbracket_aux = [(indboundgaps(i2, 1) - 1), ...
                                        (indboundgaps(i2, 2) + 1)];
                inds_gap_aux = indboundgaps(i2, 1) : 1 : indboundgaps(i2, 2);
                %
                dataADCP.u(inds_gap_aux) = interp1(dataADCP.dtime(inds_goodbracket_aux), dataADCP.u(inds_goodbracket_aux), dataADCP.dtime(inds_gap_aux));
                dataADCP.v(inds_gap_aux) = interp1(dataADCP.dtime(inds_goodbracket_aux), dataADCP.v(inds_goodbracket_aux), dataADCP.dtime(inds_gap_aux));
                dataADCP.w(inds_gap_aux) = interp1(dataADCP.dtime(inds_goodbracket_aux), dataADCP.u(inds_goodbracket_aux), dataADCP.dtime(inds_gap_aux));
                    
            end
        end
        
    end
    
    %
    disp(['Velocity data have ' num2str(length(find(isnan(dataADCP.u)))) ' NaNs.'])

    
    %% Save the data

    %
    save(fullfile(diroutputdata, ['adcpdata_' list_Signatures{i1} '.mat']), ...
         'dataADCP', '-v7.3')

    %%

    clear dataADCP sigL1
end


%%
% ----------------------------------------
% ---------- GET AQUADOPP DATA -----------
% ----------------------------------------

%% Print progress message

%
disp(' '), disp(' ')
%
disp('----- Starting to get Aquadopp data for directional estimates. Loading data from:  -----')
for i = 1:length(list_Aquadopps)
    disp([num2str(i) ') ' list_Aquadopps{i}])
end


%% Now get the data

%
for i1 = 1:length(list_Aquadopps)

    %% Print message

    disp(['---- Getting data from mooring ' list_Aquadopps{i1} ' ----'])

    %% Load data

% %     % Use different folder until I reprocess Aquadopp
% %     % and remove the fake gap because of the atm pressure issue
% %     file_aux = dir(fullfile(diroutputdata, ...
% %                             ['roxsi_aquadopp_L1_' list_Aquadopps{i1} '_*.mat']));
    %
    file_aux = dir(fullfile(dirL1data, 'Aquadopp_Level1', ...
                            ['roxsi_aquadopp_L1_' list_Aquadopps{i1} '_*.mat']));
    %
    aquadoppL1 = load(fullfile(file_aux.folder, file_aux.name));
    aquadoppL1 = aquadoppL1.aquadoppL1;


    %% Copy relevant metadata

    %
    for i2 = 1:length(list_metadatafields)
        %
        if isfield(aquadoppL1, list_metadatafields{i2})
            %
            dataADCP.(list_metadatafields{i2}) = aquadoppL1.(list_metadatafields{i2});
            %
            if strcmp(list_metadatafields{i2}, 'mooringID')
                dataADCP.instrument = "Aquadopp";
            end
        end
    end
    

    %% Select the ADCP cell that will be taken

    %
    min_bottomdepth = min(aquadoppL1.bottomdepthfrompres);

    %
    ind_binselect = find((aquadoppL1.zhab + (aquadoppL1.binsize/2)) <= min_bottomdepth, 1, 'last');
    
% %     ind_binselect = 1;    % or select the bottom-most bin for Aquadopps because it less noisy

    % If ind_binselect is empty, then this should be a shallow
    % ADCP and the first bin will be taken instead
    if isempty(ind_binselect)
        ind_binselect = 1;
    end

    %
    dataADCP.zhab = dataADCP.zhab(1:ind_binselect);
    %
    dataADCP.binN = ind_binselect;

    %
    disp(['---- Selecting velocity bin number ' num2str(dataADCP.binN) ...
          ', where hab = ' num2str(dataADCP.zhab(end), '%.1f') ' m ----'])


    %% Get the data

    %
    dataADCP.dtime = aquadoppL1.dtime;
    dataADCP.pressure = aquadoppL1.pressure;
    %
    dataADCP.bottomdepthfrompres = aquadoppL1.bottomdepthfrompres(:);
    %
    dataADCP.u = aquadoppL1.u(ind_binselect, :).';
    dataADCP.v = aquadoppL1.v(ind_binselect, :).';
    dataADCP.w = aquadoppL1.w(ind_binselect, :).';


    %% Check for NaNs

    %
    Npts_velNaN = length(find(isnan(dataADCP.u)));

    %
    disp(['Velocity data have ' num2str(Npts_velNaN) ' NaNs.'])
    %
    if Npts_velNaN>0
        warning('NaNs found in L1 data!')
% %         keyboard
    end


    %% Save the data

    %
    save(fullfile(diroutputdata, ['adcpdata_' list_Aquadopps{i1} '.mat']), ...
         'dataADCP', '-v7.3')

    %%

    clear dataADCP sigL1
end


%%
% ----------------------------------------
% ---------- GET SPOTTER??? DATA -----------
% ----------------------------------------