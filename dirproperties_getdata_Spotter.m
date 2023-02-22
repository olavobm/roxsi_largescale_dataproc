%% Get Spotter data and organize for computing
% directional spectra and properties in the
% same way as how I did for ADCPs (do it for
% Smart Mooring Spotters)

clear
close all


%%
% ---------------------------------------
% ---------------------------------------
% ---------------------------------------

%%

% % %
% % dir_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Level1/';
% % %
% % dir_output = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results/data_for_dirspectra/';

%%
%
dir_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Level1/gridded/';
%
dir_output = '/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results_3/data_for_dirspectra/';



%%

%
% % list_files = {'roxsi_spotter_L1_B01.mat', 'roxsi_spotter_L1_B03_1152.mat', 'roxsi_spotter_L1_B05_1153.mat', ...
% %               'roxsi_spotter_L1_X01_1151.mat', 'roxsi_spotter_L1_X03_1157.mat', 'roxsi_spotter_L1_X04_1155.mat', ...
% % ...%
% %               'roxsi_spotter_L1_E01_1851.mat', 'roxsi_spotter_L1_E02_1859.mat', 'roxsi_spotter_L1_E05_1853.mat', ...
% %               'roxsi_spotter_L1_E07.mat', 'roxsi_spotter_L1_E08_1852.mat', 'roxsi_spotter_L1_E09.mat', ...
% %               'roxsi_spotter_L1_E10_1848.mat', 'roxsi_spotter_L1_E11_1860.mat', 'roxsi_spotter_L1_E13_1849.mat'};

%
list_files = {'roxsi_spotter_L1_B01.mat', 'roxsi_spotter_L1_B03.mat', 'roxsi_spotter_L1_B05.mat', ...
              'roxsi_spotter_L1_X01.mat', 'roxsi_spotter_L1_X03.mat', 'roxsi_spotter_L1_X04.mat', ...
...%
              'roxsi_spotter_L1_E01.mat', 'roxsi_spotter_L1_E02.mat', 'roxsi_spotter_L1_E05.mat', ...
              'roxsi_spotter_L1_E07.mat', 'roxsi_spotter_L1_E08.mat', 'roxsi_spotter_L1_E09.mat', ...
              'roxsi_spotter_L1_E10.mat', 'roxsi_spotter_L1_E11.mat', 'roxsi_spotter_L1_E13.mat'};

%
list_mooringID = {'B01', 'B03', 'B05', 'X01', 'X03', 'X04', ...
                  'E01', 'E02', 'E05', 'E07', 'E08', 'E09', 'E10', 'E11', 'E13'};
              
              
%
if length(list_files)~=length(list_mooringID)
	error('Data pointers have not been defined consistently.') 
end


%%
% ---------------------------------------
% ---------------------------------------
% ---------------------------------------

%%

% % % % dt_spectra
% % dt_displacement = 0.4;    
% % 
% % %
% % gapTH.displacement = 0.5;
% % gapTH.location = 0.5;
% % gapTH.SST = 0.5;
% % 
% % % in seconds
% % 
% % %
% % timepars.displacement.dt = 0.4;
% % timepars.location.dt = 60;
% % timepars.SST.dt = 60;
% % 
% % %
% % timepars.displacement.gapTH = 0.5;
% % timepars.location.gapTH = 70;
% % % % timepars.SST.gapTH = 70;
% % timepars.SST.gapTH = 5*60;    % temperature on B05 did not sample as intended, so I'll increase the gap TH
% % 
% % %
% % list_subfields = fieldnames(timepars);


%%

%
for i1 = 1:length(list_files)
    
    %%
    
    %
    dataSpotter = load(fullfile(dir_data, list_files{i1}));
    dataSpotter = dataSpotter.spotterL1;
    
    % Remove existing spectra (I'll recompute)
    if isfield(dataSpotter, 'spectra')
        dataSpotter = rmfield(dataSpotter, 'spectra');
    end
    
    %
    spotterL1 = dataSpotter;

    
    %% If it's as Smart Mooring, then may need to get
    % the depth from pressure
    
    % SO FAR ONLY WORKS WITH ONE SMART MOORING FILE PER MOORING.
    % WOULD NEED TO MERGE DATA FOR SMART MOORINGS THAT WERE REPLACED!!!
    
    %
    if strcmp(list_mooringID{i1}(1), 'E')
       
% %         %
% %         dirfile_aux = dir(['/home/omarques/Documents/obm_ROXSI/obm_DataLocal' ...
% %                            '/Level1_Data/Spotter_Smart_Level1/gridded' ...
% %                            '/roxsi_smartmooring_L1_' list_mooringID{i1} 'sp_*_gridded.mat']);
        %
        dirfile_aux = dir(['/home/omarques/Documents/obm_ROXSI/obm_DataLocal' ...
                           '/Level1_Data/Spotter_Smart_Level1/merged' ...
                           '/roxsi_smartmooring_L1_' list_mooringID{i1} '.mat']);           
                       
        %
        smartmooring_aux = load(fullfile(dirfile_aux.folder, dirfile_aux.name));
        smartmooring_aux = smartmooring_aux.spotsmartL1;
        
        %
        time_window_aux = seconds(10*60);
        time_smooth_aux = spotterL1.location.dtime(1) : time_window_aux: ...
                          spotterL1.location.dtime(end);
        %
        time_smooth_aux = time_smooth_aux(2:end-1);    % trim edges
        
        %
        npts_inwindow = seconds(time_window_aux) / seconds(diff(smartmooring_aux.dtime(1:2)));
        Ndata_aux = length(smartmooring_aux.dtime);
        newdims_aux = [npts_inwindow, floor(Ndata_aux / npts_inwindow)];
        %
        indsub_aux = 1:(newdims_aux(1)*newdims_aux(2));
        
        %
        timearray_aux = reshape(smartmooring_aux.dtime(indsub_aux), newdims_aux);
        pressurearray_aux = reshape(smartmooring_aux.pressure(indsub_aux), newdims_aux);

        %
        avgtime_aux = mean(timearray_aux, 1);
        avgpressure_aux = mean(pressurearray_aux, 1, 'omitnan');


        % Computes bottom depth from average pressure
        bottomdepth_aux = (1e4*avgpressure_aux./(9.8*1025)) + smartmooring_aux.zhab;

        %
        bottomdepthavg_interp = interp1(avgtime_aux, ...
                                        bottomdepth_aux, ...
                                        spotterL1.location.dtime);
        bottomdepthavg_interp = bottomdepthavg_interp(:);
        
        %
        ind_1 = find(~isnan(bottomdepthavg_interp), 1, 'first');
        ind_2 = find(~isnan(bottomdepthavg_interp), 1, 'last');
        
        %
        spotterL1.location.z_msl = -bottomdepthavg_interp(ind_1:ind_2);
        %
        spotterL1.location.dtime = spotterL1.location.dtime(ind_1:ind_2);
        spotterL1.location.latitude = spotterL1.location.latitude(ind_1:ind_2);
        spotterL1.location.longitude = spotterL1.location.longitude(ind_1:ind_2);
        
        
    end
    
    
    %% Save data structure
    
    %
    disp(['--- Saving gridded Spotter data from mooring ' list_mooringID{i1} ' ---'])
    
    %
    save(fullfile(dir_output, ['spotterdata_' list_mooringID{i1} '.mat']), 'spotterL1', '-v7.3')

    
    %%
    
    clear spotterL1 dataSpotter;
     
end

