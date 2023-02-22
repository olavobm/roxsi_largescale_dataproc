%% Interpolate Spotter data to a regular time grids

clear
close all


%%
% ---------------------------------------
% ---------------------------------------
% ---------------------------------------

%%

%
dir_data = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Level1/not_gridded/';

%
dir_output = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level1_Data/Spotter_Level1/gridded/';


%%

%
list_files = {'roxsi_spotter_L1_B01.mat', ...
              'roxsi_spotter_L1_B03_1152.mat', 'roxsi_spotter_L1_B05_1153.mat', ...
              'roxsi_spotter_L1_X01_1151.mat', 'roxsi_spotter_L1_X03_1157.mat', 'roxsi_spotter_L1_X04_1155.mat', ...
...%
              'roxsi_spotter_L1_E01_1851.mat', 'roxsi_spotter_L1_E02_1859.mat', 'roxsi_spotter_L1_E05_1853.mat', ...
              'roxsi_spotter_L1_E07.mat', ...
              'roxsi_spotter_L1_E08_1852.mat', ...
              'roxsi_spotter_L1_E09.mat', ...
              'roxsi_spotter_L1_E10_1848.mat', 'roxsi_spotter_L1_E11_1860.mat', 'roxsi_spotter_L1_E13_1849.mat'};
          
% Corresponding mooring IDs
list_mooringID = {'B01', ...
                  'B03', 'B05', ...
                  'X01', 'X03', 'X04', ...
                  'E01', 'E02', 'E05', ...
                  'E07', ...
                  'E08', ...
                  'E09', ...
                  'E10', 'E11', 'E13'};
              

%
if length(list_files)~=length(list_mooringID)
	error('Data pointers have not been defined consistently.') 
end


%%
% ---------------------------------------
% ---- DEFINE TIME GRID/GAP CRITERIA ----
% ---------------------------------------

%%

% in seconds

%
gapTH.displacement = 0.5;
gapTH.location = 0.5;
gapTH.SST = 0.5;

%
timepars.displacement.dt = 0.4;
timepars.location.dt = 60;
timepars.SST.dt = 60;

%
timepars.displacement.gapTH = 0.5;
timepars.location.gapTH = 70;
% % timepars.SST.gapTH = 70;
timepars.SST.gapTH = 5*60;    % temperature on B05 did not sample as intended, so I'll increase the gap TH

%
list_subfields = fieldnames(timepars);


%%
% ---------------------------------------
% ---------- INTERPOLATE DATA -----------
% ---------------------------------------


%%

%
for i1 = 1:length(list_files)
    
    %%
    %
    dataSpotter = load(fullfile(dir_data, list_files{i1}));
    dataSpotter = dataSpotter.spotterL1;
    
% %     % Remove existing spectra (I'll recompute)
% %     dataSpotter = rmfield(dataSpotter, 'spectra');
    
    
    %%
    
    spotterL1.mooringID = dataSpotter.mooringID;
    spotterL1.SN = dataSpotter.SN;
    spotterL1.site = dataSpotter.site;
    spotterL1.latitude = dataSpotter.latitude;
    spotterL1.longitude = dataSpotter.longitude;
    spotterL1.X = dataSpotter.X;
    spotterL1.Y = dataSpotter.Y;
        
        
    %%
    
    % Loop over substructures where data will be interpolated
    for i2 = 1:length(list_subfields)
        
        %
        if ~isfield(dataSpotter, list_subfields{i2})
            continue
        end
        
        %%
        
        %
        time_grid_aux = dataSpotter.(list_subfields{i2}).dtime(1) : ...
                        seconds(timepars.(list_subfields{i2}).dt) : ...
                        dataSpotter.(list_subfields{i2}).dtime(end);


        % Break data apart in continuos segments (where sampling
        % time never goes larger than threshold). Then interpolate
        % data without interpolating over gaps

        %
        diff_time_seconds = seconds(diff(dataSpotter.(list_subfields{i2}).dtime));

        %
        Ndata_aux = length(dataSpotter.(list_subfields{i2}).dtime);
        
        % Find gaps
        indaboveTH = find(diff_time_seconds > timepars.(list_subfields{i2}).gapTH);

        %
        Nsegs_aux = length(indaboveTH) + 1;
        %
        indsegments_aux = NaN(Nsegs_aux, 2);
        indsegments_aux(1, 1) = 1;
        indsegments_aux(end, 2) = Ndata_aux;

        %
        for i3 = 1:(Nsegs_aux-1)
            indsegments_aux(i3, 2) = indaboveTH(i3);
            indsegments_aux(i3+1, 1) = indsegments_aux(i3, 2) + 1;
        end

        %% Assign time grid to data output structure
        % and pre-allocate space for gridded variables
        
        %
        spotterL1.(list_subfields{i2}).dtime = time_grid_aux(:);
        
        %
        list_variables_aux = fieldnames(dataSpotter.(list_subfields{i2}));
        
        %
        for i3 = 1:length(list_variables_aux)
            if ~strcmp(list_variables_aux{i3}, 'dtime')
                spotterL1.(list_subfields{i2}).(list_variables_aux{i3}) = NaN(length(time_grid_aux), 1);
            end
        end

        
        %% Now interpolate variables to time grid
        
        % Interpolate over continuous segments
        for i3 = 1:Nsegs_aux
            %
            inds_data_aux = indsegments_aux(i3, 1) : 1 : indsegments_aux(i3, 2);

            % Get grid points within the segment
            linseg_aux = (time_grid_aux >= dataSpotter.(list_subfields{i2}).dtime(inds_data_aux(1))) & ...
                         (time_grid_aux <= dataSpotter.(list_subfields{i2}).dtime(inds_data_aux(end)));

            % Interpolate over variables
            for i4 = 1:length(list_variables_aux)
                %
                if ~strcmp(list_variables_aux{i4}, 'dtime')
                    %
                    spotterL1.(list_subfields{i2}).(list_variables_aux{i4})(linseg_aux) = ...
                            interp1(dataSpotter.(list_subfields{i2}).dtime(inds_data_aux), ...
                                    dataSpotter.(list_subfields{i2}).(list_variables_aux{i4})(inds_data_aux), ...
                                    time_grid_aux(linseg_aux));

                end
            end
            
        end
    end
    
    
    %% Save output data structure
    
    %
    disp(['--- Saving gridded Spotter data from mooring ' list_mooringID{i1} ' ---'])
    
    %
    save(fullfile(dir_output, ['roxsi_spotter_L1_' list_mooringID{i1} '.mat']), 'spotterL1', '-v7.3')

    
    %%
    
    clear dataSpotter spotterL1
end