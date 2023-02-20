function sigL1 = Signature_load_L1fulldata(mooringID, dtimelims, indbins, list_vars, dirdata)
%% sigL1 = Signature_load_L1fulldata(mooringID, dtimelims, list_vars, indbins, dirdata)
%
%   inputs
%       - mooringID:
%       - dtimelims:
%       - list_vars (optional):
%       - dirdata (optional):
%
%   outputs
%       - sigL1:
%
%
% SIGNATURE_LOADL1DATA.m ......
%
% Basic functionality is operational. It does NOT load data from
% alongbeam velocity data files or backscatter.



%%

if ~exist('dirdata', 'var')
    %
    dirdata = '/project/CSIDE/ROXSI/LargeScale_Data_2022/Level1_Data/Signature_Level1/';
end


%%

if ~exist('list_vars', 'var') || isempty(list_vars)
    list_vars = {'dtime', 'pressure', 'temperature', ...
                 'heading', 'pitch', 'roll', ...
                 'bottomdepthfrompres', 'u', 'v', 'w'};
end
                      
%%

if ~exist('indbins', 'var')
    %
    lallbins = true;
else
    lallbins = false;
end


%%

list_metadata_vars = {'SN', 'mooringID', 'latitude', 'longitude', 'site', ...
                      'X', 'Y', 'coordsystem', 'magdec', ...
                      'l5beams', 'Config', ...
                      'transducerHAB', 'binsize', 'cellcenter', 'zhab', ...
                      'samplingrateHz'};


%% Get list of segment folders

%
list_segments = dir(fullfile(dirdata, 'segment*'));


%% Get time limits for each segment and whether
% there is a signature data file in that segment

%
datapointer.dirdata = dirdata;
%
datapointer.Nsegments = length(list_segments);
datapointer.mooringID = mooringID;
datapointer.lwithdata = false(datapointer.Nsegments, 1);
datapointer.dtimelims = NaT(datapointer.Nsegments, 2);

%
for i = 1:datapointer.Nsegments
    
    %
    list_fileswithID = dir(fullfile(dirdata, list_segments(i).name, ...
                                    ['roxsi_signature_L1_' mooringID '_*']));
                       
	%
    if isempty(list_fileswithID)
        continue
    end
    
    %
    datapointer.lwithdata(i) = true;
    
    %
    file_period = fullfile(dirdata, list_segments(i).name, 'segment_period.txt');
    
    %
    fid_aux = fopen(file_period, 'r');
    line_period = fgetl(fid_aux);
    
    %
    datapointer.dtimelims(i, 1) = datetime(datenum(line_period(24:42), 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum');
    datapointer.dtimelims(i, 2) = datetime(datenum(line_period(48:66), 'yyyy/mm/dd HH:MM:SS'), 'ConvertFrom', 'datenum');
    
    %
    fclose(fid_aux);
    
end

%
datapointer.dtimelims.TimeZone = 'America/Los_Angeles';


%% Return the function with datapointer if only 1 input

if nargin==1
    sigL1 = datapointer;
    return
end


%%

%
if isempty(dtimelims.TimeZone)
    dtimelims.TimeZone = 'America/Los_Angeles';
end

%
ind_first_file = find(datapointer.dtimelims(:, 1) <= dtimelims(1), 1, 'last');
% Need to deal with beginning separately
if isempty(ind_first_file)
    ind_first_file = find(datapointer.lwithdata, 1, 'first');
end

%
ind_last_file = find(datapointer.dtimelims(:, 2) >= dtimelims(2), 1, 'first');
% Need to deal with end separately
if isempty(ind_last_file)
    ind_last_file = find(datapointer.lwithdata, 1, 'last');
end

%
datapointer.indloadfiles = ind_first_file:ind_last_file;
datapointer.Nloadfiles = length(datapointer.indloadfiles);


%% Now load the data within dtimelims from the input

%
disp(['--- Loading L1 Signature data from ' num2str(datapointer.Nloadfiles) ' files ---'])

%
for i1 = 1:datapointer.Nloadfiles

    
    %%
    
    % Progress message
    disp(['Loading file ' num2str(i1) ' out of ' ...
          num2str(datapointer.Nloadfiles) ' (segment ' ...
          num2str(datapointer.indloadfiles(i1), '%.2d') ')'])
    
    %
    list_fileswithID = dir(fullfile(dirdata, ...
                                    list_segments(datapointer.indloadfiles(i1)).name, ...
                                    ['roxsi_signature_L1_' mooringID '_*.mat']));
         
	% It's assumed that the first one is the primary file (with X-Y-up
	% velocity components)
    data_aux = load(fullfile(list_fileswithID(1).folder, list_fileswithID(1).name));
    
    
% %     %% Get fields that don't change between different data files -- not really necessary
% %     
% %     %
% %     if i==1
% %         %
% %         for i2 = 1:length(list_metadata_vars)
% %             %
% %             sigL1.(list_metadata_vars{i2}) = data_aux.sigL1.(list_metadata_vars{i2});
% %         end
% %     end
    
	%%
    
    if i1==1 && lallbins
        nbins = size(data_aux.sigL1.u, 1);
        indbins = 1:nbins;
    end
    
    
    %%
    
    %
    if i1==1
        
        %
        sigL1 = data_aux.sigL1;
        sigL1 = rmfield(sigL1, 'averaged');
        
        %
        lafterstart = (sigL1.dtime >= dtimelims(1));
        
        %
        for i2 = 1:length(list_vars)
            if isvector(data_aux.sigL1.(list_vars{i2}))
                sigL1.(list_vars{i2}) = sigL1.(list_vars{i2})(lafterstart);
            else
                sigL1.(list_vars{i2}) = sigL1.(list_vars{i2})(indbins, lafterstart);
            end
        end
        
    %
    else
        
        %
        ind_firstnew_aux = 1 + find(data_aux.sigL1.dtime == sigL1.dtime(end));

        %
        for i2 = 1:length(list_vars)
            %
            if isvector(data_aux.sigL1.(list_vars{i2}))
                %
                sigL1.(list_vars{i2}) = [sigL1.(list_vars{i2}); data_aux.sigL1.(list_vars{i2})(ind_firstnew_aux:end)];
            else
                %
                sigL1.(list_vars{i2}) = [sigL1.(list_vars{i2}), data_aux.sigL1.(list_vars{i2})(indbins, ind_firstnew_aux:end)];
            end

        end
        
    end
    
    %
    if i1==datapointer.Nloadfiles
        
        %
        lbeforeend = (sigL1.dtime <= dtimelims(2));
        
        %
        for i2 = 1:length(list_vars)
            if isvector(sigL1.(list_vars{i2}))
                sigL1.(list_vars{i2}) = sigL1.(list_vars{i2})(lbeforeend);
            else
                sigL1.(list_vars{i2}) = sigL1.(list_vars{i2})(:, lbeforeend);
            end
        end
        
    end


end




