%% Group of sea-swell L2 results from all moorings
% in the same data structure. Not all results are
% selected, so read this script to check what variables
% are used


clear
close all

%%

%
% dir_data = ['/home/omarques/Documents/obm_ROXSI/Analyses/dirspectra_product/new_results/data_with_dirspectra/'];
dir_data =['/home/omarques/Documents/obm_ROXSI/Analyses' ...
           '/dirspectra_product/new_results_3/data_with_dirspectra/'];

%
dir_output = '/home/omarques/Documents/obm_ROXSI/obm_DataLocal/Level3_Data/';


%%

%
list_dir_L2_files = dir(fullfile(dir_data, 'roxsi_*_nodirspec_*.mat'));
list_dir_soloD_L2_files = dir(fullfile(dir_data, 'roxsi_solodL2_*.mat'));

%
list_allfiles = [list_dir_L2_files; list_dir_soloD_L2_files];


%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------


%%

for i = 1:length(list_allfiles)
    
	%
    data_aux = load(fullfile(list_allfiles(i).folder, list_allfiles(i).name));
    fieldname_aux = fieldnames(data_aux);
    if length(fieldname_aux)~=1
        error('!!!!')
    end
    
    %
    mooringID_aux = list_allfiles(i).name(end-6:end-4);
    
    %
    dataAll.(mooringID_aux) = data_aux.(fieldname_aux{1});
    
end

%
dataAll = orderfields(dataAll);

%
list_allmoorings = fieldnames(dataAll);


%% Define ADCPs (which have shorter time vectors)

list_ADCPs = {'A01', 'B08', 'B10', 'B11', 'B13', 'B15', ...
              'B17', 'C01', 'E03', 'E04', 'E06', 'E12', ...
              'X05', 'X06', 'X11', 'X13'};


%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------


%
list_freqbands = {'swell', 'sea', 'seaswell'};


%%

% Could incorporate this in the processing script ...

%
for i1 = 1:length(list_allmoorings)

    %
    for i2 = 1:length(list_freqbands)

        % To avoid SoloDs
        if isfield(dataAll.(list_allmoorings{i1}), list_freqbands{i2})
            continue
        end

        %
        freq_lims_aux = dataAll.(list_allmoorings{i1}).moments.(list_freqbands{i2}).freqlims;

        %
        [meanfreq_aux, peakfreq_aux] = wave_bulk_frequency(dataAll.(list_allmoorings{i1}).frequency, ...
                                                           dataAll.(list_allmoorings{i1}).See, ...
                                                           freq_lims_aux, 2);

        %
        dataAll.(list_allmoorings{i1}).moments.(list_freqbands{i2}).meanfreq = meanfreq_aux;
        dataAll.(list_allmoorings{i1}).moments.(list_freqbands{i2}).peakfreq = peakfreq_aux(:);
    end
end


%%
% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------

%%

%
dataSS.Nmoorings = length(list_allmoorings);

% -----------------
% Pre-allocate
%
dataSS.mooringID = strings(1, dataSS.Nmoorings);
dataSS.site = dataSS.mooringID;
dataSS.instrument = dataSS.mooringID;
%
dataSS.latitude = NaN(1, dataSS.Nmoorings);
dataSS.longitude = dataSS.latitude;
dataSS.X = dataSS.latitude;
dataSS.Y = dataSS.latitude;
% -----------------

%
dataSS.fftparameters.fftwindow = dataAll.B01.fftwindow;
dataSS.fftparameters.fftavgwindow = dataAll.B01.fftavgwindow;
dataSS.fftparameters.freqcutoff = dataAll.B01.freqcutoff;
dataSS.fftparameters.df = diff(dataAll.B01.frequency(1:2));

%
dataSS.dt = dataAll.B01.dt;
dataSS.dtime = dataAll.B01.dtime;
%
dataSS.frequency = dataAll.B01.frequency;


%% Add metadata and bottom depth

%
Nt = length(dataSS.dtime);

%
dataSS.bottomdepth = NaN(Nt, dataSS.Nmoorings);
% % dataSS.cp = NaN(Nt, dataSS.Nmoorings);
% % dataSS.cg = NaN(Nt, dataSS.Nmoorings);

%
for i = 1:dataSS.Nmoorings

    %
    dataSS.mooringID(i) = convertCharsToStrings(dataAll.(list_allmoorings{i}).mooringID);
    dataSS.site(i) = dataAll.(list_allmoorings{i}).site;
    %
    if ~isfield(dataAll.(list_allmoorings{i}), 'instrument')
        dataSS.instrument(i) = "Spotter";
    else
        dataSS.instrument(i) = dataAll.(list_allmoorings{i}).instrument;
    end

    %
    dataSS.latitude(i) = dataAll.(list_allmoorings{i}).latitude;
    dataSS.longitude(i) = dataAll.(list_allmoorings{i}).longitude;
    dataSS.X(i) = dataAll.(list_allmoorings{i}).X;
    dataSS.Y(i) = dataAll.(list_allmoorings{i}).Y;

    %
    if any(strcmp(list_ADCPs, list_allmoorings{i}))
        indsget = 1:length(dataAll.(list_allmoorings{i}).dtime);
        %
        ind_match = find(dataSS.dtime == dataAll.(list_allmoorings{i}).dtime(1));
        indsfill = indsget + ind_match - 1;
    else
        indsget = 1:Nt;
        indsfill = indsget;
    end
    %
    dataSS.bottomdepth(indsfill, i) = dataAll.(list_allmoorings{i}).bottomdepth(:);
    %
% %     dataSS.cp(indsfill, i) = dataAll.(list_allmoorings{i}).cp;
% %     dataSS.cg(indsfill, i) = dataAll.(list_allmoorings{i}).cg;
end


%% Add data from all moorings

%
prealloc_aux = NaN(Nt, dataSS.Nmoorings);

%
for i1 = 1:length(list_freqbands)
    %
    dataSS.(list_freqbands{i1}).freqlims = dataAll.B01.moments.(list_freqbands{i1}).freqlims;

    %
    dataSS.(list_freqbands{i1}).Hs = prealloc_aux;
    dataSS.(list_freqbands{i1}).Flux = prealloc_aux;
    dataSS.(list_freqbands{i1}).meanfreq = prealloc_aux;
    dataSS.(list_freqbands{i1}).peakfreq = prealloc_aux;
    %
    dataSS.(list_freqbands{i1}).meank = prealloc_aux;
% %     %
% %     dataSS.(list_freqbands{i1}).cp = prealloc_aux;
% %     dataSS.(list_freqbands{i1}).cg = prealloc_aux;

    %
    for i2 = 1:dataSS.Nmoorings

        %
        if any(strcmp(list_ADCPs, list_allmoorings{i2}))
            indsget = 1:length(dataAll.(list_allmoorings{i2}).dtime);
            %
            ind_match = find(dataSS.dtime == dataAll.(list_allmoorings{i2}).dtime(1));
            indsfill = indsget + ind_match - 1;
        else
            indsget = 1:Nt;
            indsfill = indsget;
        end
        

        % Add variables to output data structure
        if ~isfield(dataAll.(list_allmoorings{i2}), 'moments')
            % For SoloDs
            dataSS.(list_freqbands{i1}).Hs(indsfill, i2) = dataAll.(list_allmoorings{i2}).(list_freqbands{i1}).Hs(indsget, 1);
            dataSS.(list_freqbands{i1}).Flux(indsfill, i2) = dataAll.(list_allmoorings{i2}).(list_freqbands{i1}).Flux(indsget, 1);
            %
            dataSS.(list_freqbands{i1}).meanfreq(indsfill, i2) = dataAll.(list_allmoorings{i2}).(list_freqbands{i1}).meanfreq(indsget);
            dataSS.(list_freqbands{i1}).peakfreq(indsfill, i2) = dataAll.(list_allmoorings{i2}).(list_freqbands{i1}).peakfreq(indsget);
        else
            % For the others
            %
            dataSS.(list_freqbands{i1}).Hs(indsfill, i2) = dataAll.(list_allmoorings{i2}).EMEM.(list_freqbands{i1}).Hs(indsget, 1);
            dataSS.(list_freqbands{i1}).Flux(indsfill, i2) = dataAll.(list_allmoorings{i2}).nodirection.(list_freqbands{i1}).Flux(indsget, 1);
            %
            dataSS.(list_freqbands{i1}).meanfreq(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).meanfreq(indsget);
            dataSS.(list_freqbands{i1}).peakfreq(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).peakfreq(indsget);
        end

        %
% %         dataSS.(list_freqbands{i1}).meank = prealloc_aux;
    end
end


%% Add directional quantities

%
prealloc_aux = NaN(Nt, dataSS.Nmoorings);

%
for i1 = 1:length(list_freqbands)
    %
    dataSS.(list_freqbands{i1}).freqlims = dataAll.B01.moments.(list_freqbands{i1}).freqlims;

    %
    dataSS.(list_freqbands{i1}).a1 = prealloc_aux;
    dataSS.(list_freqbands{i1}).b1 = prealloc_aux;
    dataSS.(list_freqbands{i1}).a2 = prealloc_aux;
    dataSS.(list_freqbands{i1}).b2 = prealloc_aux;
    %
    dataSS.(list_freqbands{i1}).meandir = prealloc_aux;
    dataSS.(list_freqbands{i1}).dirspread = prealloc_aux;
    dataSS.(list_freqbands{i1}).Fx = prealloc_aux;
    dataSS.(list_freqbands{i1}).Fy = prealloc_aux;

    %
    for i2 = 1:dataSS.Nmoorings

        % Skip for SoloD
        if ~isfield(dataAll.(list_allmoorings{i2}), 'moments')
            continue
        end

        %
        if any(strcmp(list_ADCPs, list_allmoorings{i2}))
            indsget = 1:length(dataAll.(list_allmoorings{i2}).dtime);
            %
            ind_match = find(dataSS.dtime == dataAll.(list_allmoorings{i2}).dtime(1));
            indsfill = indsget + ind_match - 1;
        else
            indsget = 1:Nt;
            indsfill = indsget;
        end
        
        %
        dataSS.(list_freqbands{i1}).a1(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).a1(indsget);
        dataSS.(list_freqbands{i1}).b1(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).b1(indsget);
        dataSS.(list_freqbands{i1}).a2(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).a2(indsget);
        dataSS.(list_freqbands{i1}).b2(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).b2(indsget);
        %
        dataSS.(list_freqbands{i1}).meandir(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).meandir(indsget);
        dataSS.(list_freqbands{i1}).dirspread(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).dirspread(indsget);
        dataSS.(list_freqbands{i1}).Fx(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).Fx(indsget);
        dataSS.(list_freqbands{i1}).Fy(indsfill, i2) = dataAll.(list_allmoorings{i2}).moments.(list_freqbands{i1}).Fy(indsget);


    end
end


%%

%
save(fullfile(dir_output, 'largescale_seaswell_ROXSI2022.mat'), 'dataSS', '-v7.3')

%
disp('--- DONE GROUPING ALL OF THE L2 DATA ---')




%%
% -------------------------------------------------------
% -------------------------------------------------------
% -------------------------------------------------------



%%

% % INCORPORTATE SOMETHING ABOUT SELECTING ONLY THE CONVIENT DATA???

% % 
% % %
% % % % % % list_moorings = ["B01", "B03", "B05", "E08", "B08", "B10", "B11", "B13", "B15"];
% % % % list_moorings = ["B01", "B03", "B05", "E08", "B10", "B13", "B15"];    % no Aquadopps
% % 
% % %
% % % list_moorings = ["B01", "B03", "B05", "E05", "E08", "B10", "B13", "B15", "A01", "C01"];
% % list_moorings = ["B01", "B03", "B05", ...
% %                  "E01", "E02", "E05", "E07", "E08", "E09", "E10", "E11", "E13", ...
% %                  "B10", "B13", "B15", "A01", "C01"];
% % 
% % list_moorings = ["B01", "B03", "B05", ...
% %                  "E01", "E02", "E05", "E07", "E08", "E09", "E10", "E11", "E13", ...
% %                  "B10", "B11", "B13", "B15", "A01", "C01"];
% % 
% % % % % Asilomar only
% % % % list_moorings = ["X01", "X03", "X04", "X05", "X06", "X11", "X13"];
% % % % list_moorings = ["X01", "X03", "X04", "X05"];
% % % % list_moorings = ["X04", "X05"];
% % 
% % 
% % %
% % % listband_pick = 'swell';
% % listband_pick = 'sea';
% % % listband_pick = 'seaswell';
% % 
% % % -----------
% % %
% % lpick = strcmp(dataSS.mooringID, "X05");
% % %
% % dataSS.(listband_pick).meandir(:, lpick) = dataSS.(listband_pick).meandir(:, lpick) - (58*pi/180);
% % % -----------
% % 
% % 
% % %
% % lgetmoor = false(dataSS.Nmoorings, 1);
% % %
% % for i = 1:dataSS.Nmoorings
% %     %
% %     if any(strcmp(list_moorings, dataSS.mooringID(i)))
% %         lgetmoor(i) = true;
% %     end
% % end
% % 
% % %
% % Hs_TH = 0.5;
% % 
% % %
% % lcommondata = ~isnan(mean(dataSS.(listband_pick).Hs(:, lgetmoor), 2));
% % ltimeavg = dataSS.(listband_pick).Hs(:, strcmp(dataSS.mooringID, "B01")) > Hs_TH;
% % %
% % lplot = lcommondata & ltimeavg;

