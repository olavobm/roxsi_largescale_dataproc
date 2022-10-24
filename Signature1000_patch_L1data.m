function Signature1000_patch_L1data(dirsplitdata, list_Signature, diroutput)
%% SIGNATURE1000_PATCH_L1DATA(dirsplitdata, list_Signature, diroutput)
%
%   inputs
%       -
%       -
%       -
%
%   outputs
%       -
%


%%

data_file_names = {'roxsi_signature_L1_', '.mat'; ...
                   'roxsi_signature_L1_', '_scalars.mat'; ...
                   'roxsi_signature_L1_', '_alongbeamvel.mat'};

%%

list_segments = dir(fullfile(dirsplitdata, 'segment_*'));

% this assumes the numbers are sequential in filename. That is,
% if there are, e.g., 15 segments, folders should be segment_01,
% segment_02, etc, instead of segment_1 and segment_2.


%%

% Loop over Signatures in the list
for i1 = 1:length(list_Signature)

    % Loop over types of files
    for i2 = 1:size(data_file_names, 1)

        % First load data in the first segment, and other segments
        % will be patched to it
        data_seg_all = load(fullfile(list_segments(1).folder, ...
                          list_segments(1).name, ...
                          [data_file_names{i2, 1}, list_Signature{i1}, data_file_names{i2, 2}]));
        %
        list_variables = fieldnames(data_seg_all);

        % If it's the beam data, find the structure variable that has the datetime in it
        if (i2==1) || (i2==2)
            varL1struct = 'sigL1';
        elseif (i2==3)
            varL1struct = 'sigL1metadata';
        else

        end

        keyboard
        % Loop over segment/folders (starting at the second segment)
        for i3 = 2:length(list_segments)

            %
            data_seg_next = load(fullfile(list_segments(i3).folder, ...
                                 list_segments(i3).name, ...
                                 [data_file_names{i2, 1}, list_Signature{i1}, data_file_names{i2, 2}]));

            % Find the timestamp on the next segment that is the
            % last one in the patched data segment
            ind_match_aux = find(data_seg_next.(varL1struct).dtime == ...
                                 data_seg_all.(varL1struct).dtime(end));

            %
            Nlengthtime = length(data_seg_next.(varL1struct).dtime);

            %
            for i4 = 1:length(list_variables)
                %
                if strcmp(varL1struct, list_variables{i4})

                    %
                    list_fields_in_struct = fieldnames(data_seg_next.(varL1struct));

                    %
                    for i5 = 1:length(list_fields_in_struct)
                        
                        % If it's a variable with time in the ROW dimension
                        if size(data_seg_next.(varL1struct).(list_fields_in_struct{i5}), 1) == Nlengthtime
% %                             cat_dim = 1;

                            % Concatenate averaged variables variables
                            data_seg_all.(varL1struct).(list_fields_in_struct{i5}) = ...
                                        [data_seg_all.(varL1struct).(list_fields_in_struct{i5}); ...
                                         data_seg_next.(varL1struct).(list_fields_in_struct{i5})((ind_match_aux+1):end, :)];

                        % If it's a variable with time in the COLUMN dimension
                        elseif size(data_seg_next.(varL1struct).(list_fields_in_struct{i5}), 2) == Nlengthtime
%                             cat_dim = 2;

                            data_seg_all.(varL1struct).(list_fields_in_struct{i5}) = ...
                                        [data_seg_all.(varL1struct).(list_fields_in_struct{i5}), ...
                                         data_seg_next.(varL1struct).(list_fields_in_struct{i5})(:, (ind_match_aux+1):end)];

                        end

                    

                    end

                %
                else

                    % If it's a variable with time in the ROW dimension
                    if size(data_seg_next.(list_variables{i4}), 1) == Nlengthtime
                        %
                        data_seg_all.(list_variables{i4}) = ...
                                            [data_seg_all.(list_variables{i4}); ...
                                             data_seg_next.(list_variables{i4})((ind_match_aux+1):end, :)];

                    % If it's a variable with time in the COLUMN dimension
                    elseif size(data_seg_next.(list_variables{i4}), 2) == Nlengthtime

                        data_seg_all.(list_variables{i4}) = ...
                                            [data_seg_all.(list_variables{i4}), ...
                                             data_seg_next.(list_variables{i4})(:, (ind_match_aux+1):end)];
                    end
                    
                end

            end


            % Do the averaged fields now
            if (i2==1) || (i2==2)
                %
                ind_match_avg_aux = find(data_seg_next.sigL1.averaged.dtime == ...
                                         data_seg_all.sigL1.averaged.dtime(end));
                %
                list_fields_in_struct = fieldnames(data_seg_all.sigL1.averaged);

                %
                Nlengthavgdtime = length(data_seg_next.sigL1.averaged.dtime);

                %
                for i4 = 1:length(list_fields_in_struct)

                    % If it's a variable with time in the ROW dimension
                    if size(data_seg_next.sigL1.averaged.(list_fields_in_struct{i4}), 1) == Nlengthavgdtime
                        
                        %
                        data_seg_all.sigL1.averaged.(list_fields_in_struct{i4}) = ...
                                    [data_seg_all.sigL1.averaged.(list_fields_in_struct{i4}); ...
                                     data_seg_next.sigL1.averaged.(list_fields_in_struct{i4})((ind_match_avg_aux+1):end, :)];


                    % If it's a variable with time in the COLUMN dimension
                    elseif size(data_seg_next.sigL1.averaged.(list_fields_in_struct{i4}), 2) == Nlengthavgdtime
                        

                        % Concatenate averaged variables variables
                        data_seg_all.sigL1.averaged.(list_fields_in_struct{i4}) = ...
                                    [data_seg_all.sigL1.averaged.(list_fields_in_struct{i4}), ...
                                     data_seg_next.sigL1.averaged.(list_fields_in_struct{i4})(:, (ind_match_avg_aux+1):end)];

                    end

                    
                end
                

            end

            %
            clear data_seg_next
        end   

    
        %% Now save patched data
    
        %
        save(fullfile(diroutput, ...
                      [data_file_names{i2, 1}, list_Signature{i1}, data_file_names{i2, 2}]), ...
                      '-struct', 'data_seg_all')
        
        
        %%
        clear data_seg_all

    end



    
end



