clear
close all


%%
% ----------------------------------------
% -------------- LOAD DATA ---------------
% ----------------------------------------

%%

%
% % dir_parent_data = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/figures_bydate/2022_10_26/Signature_L2/Signature_Level2/';


%% Load SoloD L2

%
dir_soloD_L2 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/figures_bydate/2022_10_26/SoloD_L2/';

% % %
% % dataspec.soloD_B14 = load(fullfile(dir_soloD, 'roxsi_soloD_L2_B14.mat'));
% % dataspec.soloD_B14 = dataspec.soloD_B14.SOLOD_L2;

%
list_soloDL2 = dir(fullfile(dir_soloD_L2, 'roxsi_soloD_L2_*.mat'));

%
for i = 1:length(list_soloDL2)

    %
    data_aux = load(fullfile(list_soloDL2(i).folder, list_soloDL2(i).name));
    dataspec.(['soloD_' list_soloDL2(i).name(16:18)]) = data_aux.SOLOD_L2;

end

% Remove C09 from the data
dataspec = rmfield(dataspec, 'soloD_C09');

% return

%% Load Aquadopp L2

%
dir_aquadopp_L2 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/figures_bydate/2022_10_26/Aquadopp_L2/Aquadopp_Level2/';
list_aquadoppL2 = dir(fullfile(dir_aquadopp_L2, 'roxsi_aquadopp_L2_*.mat'));

%
for i = 1:length(list_aquadoppL2)

    %
    data_aux = load(fullfile(list_aquadoppL2(i).folder, list_aquadoppL2(i).name));

    % Wave statistics was only computed for Aquadopps sampling at 1 Hz
    if (data_aux.aquadoppL2.samplingtime == 1)
        %
        dataspec.(['aquadopp_' list_aquadoppL2(i).name(19:21)]) = data_aux.aquadoppL2;
    end

end


%% Load Signature L2

%
dir_signature_L2 = '/Volumes/LaCie/ROXSI/LargeScale_Data_2022/Level2_Data/Signature_Level2/';
%
list_signatureL2 = dir(fullfile(dir_signature_L2, 'roxsi_signature_L2_X*.mat'));

%
for i = 1:length(list_signatureL2)

    %
    data_aux = load(fullfile(list_signatureL2(i).folder, list_signatureL2(i).name));
    dataspec.(['signature_' list_signatureL2(i).name(20:22)]) = data_aux.sigL2;

end


%% Load Spotter L1

%
% % dir_spotter_L1 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/figures_bydate/2022_10_26/';
dir_spotter_L1 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/figures_bydate/2022_10_26/Spotter_L2/Spotter_Level2/';

%
list_spotterL1 = dir(fullfile(dir_spotter_L1, 'roxsi_spotter_L2_*.mat'));

%
for i = 3:length(list_spotterL1)

    %
    data_aux = load(fullfile(list_spotterL1(i).folder, list_spotterL1(i).name));
% %     dataspec.(['spotter_' list_spotterL1(i).name(18:20)]) = data_aux.spotterL1;
    dataspec.(['spotter_' list_spotterL1(i).name(18:20)]) = data_aux.spotterL2;

end

%%

% % %
% % list_signature.B10.filename = 'roxsi_signature_L2_B10_103045.mat';
% % list_signature.B13.filename = 'roxsi_signature_L2_B13_103046.mat';
% % list_signature.B15.filename = 'roxsi_signature_L2_B15_103056.mat';
% % list_signature.B17.filename = 'roxsi_signature_L2_B17_101923.mat';
% % 
% % %
% % list_fiels = fieldnames(list_signature);
% % 
% % %
% % for i = 1:length(list_fiels)
% %     dataspec.(list_fiels{i}) = load(fullfile(dir_sigL2, list_signature.(list_fiels{i}).filename), 'sigL2');
% %     dataspec.(list_fiels{i}) = dataspec.(list_fiels{i}).sigL2;
% % end

%%
% % keyboard


%% Load Duet at Asilomar(???) L2

%% Coastline

%
coastlineMP = load('/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/code_proc/coastline_Monterey.mat');
coastlineMP = coastlineMP.coastLineMP;

% Trim

% X/Y coordinates China Rock

% X/Y coordinates Asilomar

%% Bathymetry???


%%

%
list_alldata = fieldnames(dataspec);

%%
% ----------------------------------------
% ---- ........................ ----------
% ----------------------------------------

%%
%
for i1 = 1:length(list_alldata)

    % Identify type of instrument
    ind_underscore_aux = strfind(list_alldata{i1}, '_');
    %
    strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    %
    if strcmp(strinstr_aux, 'soloD')

        dataspec.(list_alldata{i1}).dtime = ...
                        datetime(dataspec.(list_alldata{i1}).time_dnum, 'ConvertFrom', 'datenum');
        dataspec.(list_alldata{i1}).dtime.TimeZone = 'America/Los_Angeles';

        dataspec.(list_alldata{i1}).See = dataspec.(list_alldata{i1}).See.';
    end
   
end

%%

%
for i1 = 1:length(list_alldata)

    % Identify type of instrument
    ind_underscore_aux = strfind(list_alldata{i1}, '_');
    %
    strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    %
    if strcmp(strinstr_aux, 'aquadopp')

        dataspec.(list_alldata{i1}).meanfreqSS = 1./ dataspec.(list_alldata{i1}).TmeanSS;

        %
        dataspec.(list_alldata{i1}).meanfreqSS = dataspec.(list_alldata{i1}).meanfreqSS(:);
        dataspec.(list_alldata{i1}).HsigSS = dataspec.(list_alldata{i1}).HsigSS(:);

        %
        dataspec.(list_alldata{i1}).bottomdepthfrompres = dataspec.(list_alldata{i1}).pressuremean;
        dataspec.(list_alldata{i1}).bottomdepthfrompres = dataspec.(list_alldata{i1}).bottomdepthfrompres(:);
        dataspec.(list_alldata{i1}).bottomdepthfrompres = 1e4*dataspec.(list_alldata{i1}).bottomdepthfrompres/(1030*9.8);
    end

    %
    if strcmp(strinstr_aux, 'signature')

        dataspec.(list_alldata{i1}).meanfreqSS = 1./ dataspec.(list_alldata{i1}).TmeanSS;

        %
        dataspec.(list_alldata{i1}).meanfreqSS = dataspec.(list_alldata{i1}).meanfreqSS(:);
        dataspec.(list_alldata{i1}).HsigSS = dataspec.(list_alldata{i1}).HsigSS(:);

        %
        dataspec.(list_alldata{i1}).bottomdepthmean = dataspec.(list_alldata{i1}).bottomdepthmean(:);

    end
   
end


%%
% ----------------------------------------
% --- LOAD AND PATCH RESULTS FROM THE ----
% ----------- TWO B01 SPOTTERS -----------
% ----------------------------------------

%%


% % dir_spotter_L1 = '/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/figures_bydate/2022_10_26/Spotter_L2/Spotter_Level2/';
% % 
% % %
% % spotter_B01_1150 = load(fullfile(dir_spotter_L1, 'roxsi_spotter_L2_B01_1150_reduced.mat'));
% % spotter_B01_1158 = load(fullfile(dir_spotter_L1, 'roxsi_spotter_L2_B01_1158_reduced.mat'));
% % %
% % spotter_B01_1150 = spotter_B01_1150.spotterL2;
% % spotter_B01_1158 = spotter_B01_1158.spotterL2;
% % 
% % %
% % dataspec.spotter_B01 = spotter_B01_1158;
% % %
% % time_mid_aux = datetime((datenum(spotter_B01_1158.dtime(end)) + ...
% %                         datenum(spotter_B01_1150.dtime(1)))./2, ...
% %                         'ConvertFrom', 'datenum', 'TimeZone', 'America/Los_Angeles');
% % %
% % dataspec.spotter_B01.dtime = [dataspec.spotter_B01.dtime; time_mid_aux; spotter_B01_1150.dtime];
% % dataspec.spotter_B01.See = [dataspec.spotter_B01.See, NaN(length(spotter_B01_1158.frequency), 1), spotter_B01_1150.See];
% % dataspec.spotter_B01.meanfreq = [dataspec.spotter_B01.meanfreq; NaN; spotter_B01_1150.meanfreq];
% % dataspec.spotter_B01.Hsig = [dataspec.spotter_B01.Hsig; NaN; spotter_B01_1150.Hsig];
% % dataspec.spotter_B01.bottomdepth = [dataspec.spotter_B01.bottomdepth, NaN, spotter_B01_1150.bottomdepth];
% % 
% % % Update this
% % list_alldata = fieldnames(dataspec);

%%
% ----------------------------------------
% --- LOAD AND PATCH RESULTS FROM THE ----
% ----------- TWO B01 SPOTTERS -----------
% ----------------------------------------



%% Load Smart mooring as Spotter

% % %
% % data_aux = load(fullfile(dir_spotter_L1, 'roxsi_spotter_L1_E08_1852.mat'));
% % dataspec.spotter_E08 = data_aux.spotterL1;
% % 
% % %
% % dataspec.spotter_E08.bottomdepth


%%
% ----------------------------------------
% ---- ........................ ----------
% ----------------------------------------

%%

%
varfieldnames.soloD.bottomdepth = 'mean_depth';
varfieldnames.aquadopp.bottomdepth = 'bottomdepthfrompres';
varfieldnames.signature.bottomdepth = 'bottomdepthmean';
varfieldnames.spotter.bottomdepth = 'bottomdepth';

%
varfieldnames.soloD.timespec = 'dtime';
varfieldnames.aquadopp.timespec = 'dtime';
varfieldnames.signature.timespec = 'dtime';
varfieldnames.spotter.timespec = 'dtime';

%
varfieldnames.soloD.frequency = 'freq';
varfieldnames.aquadopp.frequency = 'frequency';
varfieldnames.signature.frequency = 'frequency';
varfieldnames.spotter.frequency = 'frequency';

%
varfieldnames.soloD.HsigSS = 'HsigSS';
varfieldnames.aquadopp.HsigSS = 'HsigSS';
varfieldnames.signature.HsigSS = 'HsigSS';
% varfieldnames.spotter.HsigSS = 'HsigSS';
varfieldnames.spotter.HsigSS = 'Hsig';

%
varfieldnames.soloD.meanfreqSS = 'meanfreqSS';
varfieldnames.aquadopp.meanfreqSS = 'meanfreqSS';
varfieldnames.signature.meanfreqSS = 'meanfreqSS';
% varfieldnames.spotter.meanfreqSS = 'meanfreqSS';
varfieldnames.spotter.meanfreqSS = 'meanfreq';


%%
% ----------------------------------------
% ---- COMPUTE MEAN WAVENUMBER USING -----
% ---------- LINEAR WAVE THEORY ----------
% ----------------------------------------

%%

%
g = 9.8;
rho0 = 1030;

%%

% (PS: frequency in Hz and k in radians per meter)
disp_rel = @(k, freq, H) g*k*tanh(k*H) - (2*pi*freq)^2;


%% Compute mean wavenumber correspondent to the mean frequency

%
for i1 = 1:length(list_alldata)
% 
    % Identify type of instrument
    ind_underscore_aux = strfind(list_alldata{i1}, '_');
    %
    strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    % Make bottom depth a column vector and positive
    dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth) = ...
            abs(dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth)(:));
    %
    bottomdepth_aux = dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth);
    meanfreq_aux = dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).meanfreqSS);

    %
    dataspec.(list_alldata{i1}).kmeanSS = NaN(length(bottomdepth_aux), 1);

    % Loop over time (i.e. over varying bottom depth)
    for i2 = 1:length(bottomdepth_aux)

        %
        if (bottomdepth_aux(i2) > 1) && ~isnan(meanfreq_aux(i2))
            %
            disp_rel_eval = @(k) disp_rel(k, meanfreq_aux(i2), bottomdepth_aux(i2));
            %
            meank_aux = fzero(disp_rel_eval, [(2*pi/(5000)), (2*pi/(1))]);

            %
            dataspec.(list_alldata{i1}).kmeanSS(i2) = meank_aux;

        end
    end
end


%%
% ----------------------------------------
% -- COMPUTE cg, ENERGY, and ENERGY FLUX -
% ----------------------------------------

%% Compute group velocity

%
for i1 = 1:length(list_alldata)

    ind_underscore_aux = strfind(list_alldata{i1}, '_');
    strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    %
    k_aux = dataspec.(list_alldata{i1}).kmeanSS;

    %
    dataspec.(list_alldata{i1}).kH = k_aux .* ...
                                   dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth);

    %
    dataspec.(list_alldata{i1}).cp = sqrt(g * tanh(dataspec.(list_alldata{i1}).kH) ./ k_aux);
    
    %
    dataspec.(list_alldata{i1}).cg = dataspec.(list_alldata{i1}).cp .* ...
                                   0.5 .* (1 + (2*dataspec.(list_alldata{i1}).kH./sinh(2*dataspec.(list_alldata{i1}).kH)));

end


%% Recompute Hsig and mean frequency for a noncontaminated frequency band

%
new_freqband = [0.045, 0.2];

%
for i1 = 1:length(list_alldata)

    %
    ind_underscore_aux = strfind(list_alldata{i1}, '_');
    strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    %
    dtime_aux = dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).timespec);
    frequency_aux = dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).frequency);
    frequency_aux = frequency_aux(:);
    See_aux = dataspec.(list_alldata{i1}).See;
    %
    linfreq_aux = (frequency_aux >= new_freqband(1)) & ...
                  (frequency_aux <= new_freqband(2));
    %
    m0_aux = trapz(frequency_aux(linfreq_aux), See_aux(linfreq_aux, :), 1);
    m1_aux = trapz(frequency_aux(linfreq_aux), ...
                                (See_aux(linfreq_aux, :) .* ...
                                 repmat(frequency_aux(linfreq_aux), 1, length(dtime_aux))), 1);

    %
    dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).meanfreqSS) = 1./(m0_aux./m1_aux);
    dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).HsigSS) = 4*sqrt(m0_aux);
end



%% Compute depth-integrated total energy

%
for i1 = 1:length(list_alldata)
   
    ind_underscore_aux = strfind(list_alldata{i1}, '_');
    strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    %
    dataspec.(list_alldata{i1}).E = (1/16) * rho0 * g * ...
                        (dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).HsigSS)).^2;
    dataspec.(list_alldata{i1}).E = dataspec.(list_alldata{i1}).E(:);

end


%% Compute depth-integrated energy flux

%
lgooddims = true(length(list_alldata), 1);

%
for i = 1:length(list_alldata)
    
    %
    if ~isequal(size(dataspec.(list_alldata{i}).E), ...
                size(dataspec.(list_alldata{i}).cg))
        
        %
        lgooddims(i) = false;
    else

        %
        dataspec.(list_alldata{i}).FluxSS = dataspec.(list_alldata{i}).E .* ...
                                            dataspec.(list_alldata{i}).cg;
    end
end


% % 
% % %%
% % % ----------------------------------------
% % % ------------- CALCULATIONS -------------
% % % ----------------------------------------
% % %
% % % Depth-integrated energy 
% % 
% % %%
% % 
% % for i = 1:length(list_fiels)
% %     %
% %     dataspec.(list_fiels{i}).Energy = (1/8)*1030*9.8 * ...
% %                                     (dataspec.(list_fiels{i}).HsigSS).^2;
% %     %
% %     dataspec.(list_fiels{i}).Flux = sqrt(9.8*dataspec.(list_fiels{i}).bottomdepthmean) .* ...
% %                                     dataspec.(list_fiels{i}).Energy;
% % end
% % 

%%
% ----------------------------------------------------------------
% ---- TIMEGRID??? AVERAGE??? SORT RESULTS BY BOTTOM DEPTH??? ----
% ----------------------------------------------------------------

%%

%
datagridded.dtime = datetime(2022, 06, 17, 20, 00, 00) : hours(1) : ...
                    datetime(2022, 07, 20, 02, 00, 00);
datagridded.dtime.TimeZone = 'America/Los_Angeles';
%
datagridded.dtime = datagridded.dtime(:);


%% Results along B line

%
ind_Bline = 0;

%
for i1 = 1:length(list_alldata)

    % Identify type of instrument
    ind_underscore_aux = strfind(list_alldata{i1}, '_');

    %
    if strcmp(list_alldata{i1}(ind_underscore_aux+1), 'B')
        ind_Bline = ind_Bline + 1;
    end
end

%
datagridded.Bline.mooringID = strings(ind_Bline, 1);
datagridded.Bline.instrument = strings(ind_Bline, 1);
%
datagridded.Bline.X = NaN(ind_Bline, 1);
datagridded.Bline.Y = NaN(ind_Bline, 1);
datagridded.Bline.meandepth = NaN(ind_Bline, 1);
%
datagridded.Bline.FluxSS = NaN(ind_Bline, length(datagridded.dtime));


%
ind_Bline = 0;
%
for i = 1:length(list_alldata)

    % Identify type of instrument
    ind_underscore_aux = strfind(list_alldata{i}, '_');

    %
    if strcmp(list_alldata{i}(ind_underscore_aux+1), 'B')

        %
        ind_Bline = ind_Bline + 1;

        %
        datagridded.Bline.mooringID(ind_Bline) = convertCharsToStrings(list_alldata{i}((end-2):end));
        instrument_char_aux = list_alldata{i}(1:(ind_underscore_aux-1));
        datagridded.Bline.instrument(ind_Bline) = convertCharsToStrings(instrument_char_aux);

        %
        datagridded.Bline.X(ind_Bline) = dataspec.(list_alldata{i}).X;
        datagridded.Bline.Y(ind_Bline) = dataspec.(list_alldata{i}).Y;
        datagridded.Bline.meandepth(ind_Bline) = mean(dataspec.(list_alldata{i}).(varfieldnames.(instrument_char_aux).bottomdepth), 'omitnan');

        %
        datagridded.Bline.HsigSS(ind_Bline, :) = ...
                            interp1(dataspec.(list_alldata{i}).dtime, ...
                                    dataspec.(list_alldata{i}).(varfieldnames.(instrument_char_aux).HsigSS), ...   % same name for all??
                                    datagridded.dtime);

        %
        datagridded.Bline.FluxSS(ind_Bline, :) = ...
                            interp1(dataspec.(list_alldata{i}).dtime, ...
                                    dataspec.(list_alldata{i}).FluxSS, ...
                                    datagridded.dtime);
    end
end


%% Interpolate spectra and get average spectrum for each instrument


%
lavg_spec = true(1, length(datagridded.dtime));

%
for i = 1:length(list_alldata)

    % Identify type of instrument
    ind_underscore_aux = strfind(list_alldata{i}, '_');

    %
    if strcmp(list_alldata{i}(ind_underscore_aux+1), 'B')

        %
        instrument_char_aux = list_alldata{i}(1:(ind_underscore_aux-1));

        %
        freq_aux = dataspec.(list_alldata{i}).(varfieldnames.(instrument_char_aux).frequency);
        %
        See_aux = NaN(length(freq_aux), length(datagridded.dtime));

        %
        tic
        for i2 = 1:length(freq_aux)

            %
            See_aux(i2, :) = interp1(dataspec.(list_alldata{i}).dtime, ...
                                     dataspec.(list_alldata{i}).See(i2, :), ...
                                     datagridded.dtime);

            

        end
        toc

        %
        lavg_spec = lavg_spec & ~isnan(mean(See_aux, 1, 'omitnan'));

        %
        datagridded.Bline.avgspec.(list_alldata{i}).frequency = freq_aux(:);
        datagridded.Bline.avgspec.(list_alldata{i}).See = See_aux;
        
    end
end

% ------------------------------------------------
%
for i = 1:length(list_alldata)

    % Identify type of instrument
    ind_underscore_aux = strfind(list_alldata{i}, '_');

    %
    if strcmp(list_alldata{i}(ind_underscore_aux+1), 'X')

        %
        instrument_char_aux = list_alldata{i}(1:(ind_underscore_aux-1));

        %
        datagridded.Bline.avgspec.(list_alldata{i}).Seeavg = ...
            mean(datagridded.Bline.avgspec.(list_alldata{i}).See(:, lavg_spec), 2);
    end
end



%%
% ----------------------------------------
% --- SORT THE BY CROSS-SHORE DISTANCE ---
% ------ AND ORGANIZE FINAL RESULTS ------
% ----------------------------------------

%%


%% Sort by cross-shore location, but not really necessary apparently

%
[~, ind_sortB] = sort(datagridded.Bline.X);

%
datagridded.Bline.mooringID = datagridded.Bline.mooringID(ind_sortB);
datagridded.Bline.instrument = datagridded.Bline.instrument(ind_sortB);
%
datagridded.Bline.X = datagridded.Bline.X(ind_sortB, :);
datagridded.Bline.Y = datagridded.Bline.Y(ind_sortB, :);
datagridded.Bline.meandepth = datagridded.Bline.meandepth(ind_sortB, :);
%
datagridded.Bline.HsigSS = datagridded.Bline.HsigSS(ind_sortB, :);
datagridded.Bline.FluxSS = datagridded.Bline.FluxSS(ind_sortB, :);


%%
%
datagridded.Bline.timeavg.timelims = [datetime(2022, 06, 23), datetime(2022, 07, 20)];
datagridded.Bline.timeavg.timelims.TimeZone = 'America/Los_Angeles';
%
% % datagridded.Bline.timeavg.lintimeavg = (datagridded.dtime >= datagridded.Bline.timeavg.timelims(1)) & ...
% %                                        (datagridded.dtime <= datagridded.Bline.timeavg.timelims(2));
datagridded.Bline.timeavg.lintimeavg = ~isnan(mean(datagridded.Bline.FluxSS, 1));

datagridded.Bline.timeavg.HsigSS = mean(datagridded.Bline.HsigSS(:, datagridded.Bline.timeavg.lintimeavg), 2, 'omitnan');
datagridded.Bline.timeavg.FluxSS = mean(datagridded.Bline.FluxSS(:, datagridded.Bline.timeavg.lintimeavg), 2, 'omitnan');


%%

%
datagridded.Bline.FluxSSfraction = (datagridded.Bline.FluxSS ./ ...
                                    repmat(datagridded.Bline.FluxSS(1, :), length(datagridded.Bline.mooringID), 1));
%
datagridded.Bline.timeavg.FluxSSfraction = (datagridded.Bline.timeavg.FluxSS ./ ...
                                             datagridded.Bline.timeavg.FluxSS(1));


% % %%
% % return


%% Results along X line

%%
% ------------------------------------------------------
% ------------------------------------------------------
% -------------------- MAKE FIGURES --------------------
% ------------------------------------------------------
% ------------------------------------------------------


% % % %% Plot pressure spectra with contaminations from SoloDs
% % % 
% % % %
% % % figure
% % %     %
% % %     hold on
% % %     %
% % %     plot(dataspec.soloD_B04.freq, ...
% % %          mean(dataspec.soloD_B04.See, 2, 'omitnan'))
% % %     plot(dataspec.soloD_B06.freq, ...
% % %          mean(dataspec.soloD_B06.See, 2, 'omitnan'))
% % %     plot(dataspec.soloD_B07.freq, ...
% % %          mean(dataspec.soloD_B07.See, 2, 'omitnan'))
% % %     plot(dataspec.soloD_B09.freq, ...
% % %          mean(dataspec.soloD_B09.See, 2, 'omitnan'))
% % %     plot(dataspec.soloD_B12.freq, ...
% % %          mean(dataspec.soloD_B12.See, 2, 'omitnan'))
% % % 
% % %     %
% % %     set(gca, 'FontSize', 16, 'Box', 'on', ...
% % %              'XGrid', 'on', 'YGrid', 'on', ...
% % %              'YScale', 'log', 'YLim', [6e-2, 3])
% % % 
% % %     %
% % %     overlayline('v', 0.3', '-k', 'LineWidth', 2)
% % % 
% % %     %
% % %     xlabel('Frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     ylabel('Elevation variance [m$^2$ Hz$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     %
% % %     title('Average spectrum: SoloD''s B04, B06, B07, B09, B12', 'Interpreter', 'Latex', 'FontSize', 16)
% % % 
% % % 
% % % %% Plot average spectrum from all instruments on B line
% % % 
% % % %
% % % list_fileds_instr = fieldnames(datagridded.Bline.avgspec);
% % % 
% % % %
% % % newFigDims([11.6, 7.32])
% % %     %
% % %     hold on
% % %     %
% % %     for i = 1:length(list_fileds_instr)
% % %         plot(datagridded.Bline.avgspec.(list_fileds_instr{i}).frequency, ...
% % %              datagridded.Bline.avgspec.(list_fileds_instr{i}).Seeavg)
% % %     end
% % %     %
% % %     plot(datagridded.Bline.avgspec.soloD_B04.frequency, ...
% % %          datagridded.Bline.avgspec.soloD_B04.Seeavg, '-k', 'LineWidth', 2)
% % %     plot(datagridded.Bline.avgspec.soloD_B07.frequency, ...
% % %          datagridded.Bline.avgspec.soloD_B07.Seeavg, '-k', 'LineWidth', 2)
% % % 
% % %     %
% % %     set(gca, 'FontSize', 16, 'Box', 'on', ...
% % %              'XGrid', 'on', 'YGrid', 'on', ...
% % %              'YScale', 'log', 'YLim', [1e-2, 1.5], ...
% % %              'XLim', [1e-10, 0.31])
% % % 
% % %     %
% % %     overlayline('v', 0.3, '--k', 'LineWidth', 1)
% % %     overlayline('v', 0.2, '--r', 'LineWidth', 2)
% % % 
% % %     %
% % %     xlabel('Frequency [Hz]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     ylabel('Elevation variance [m$^2$ Hz$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', 16)
% % %     %
% % %     title('Average spectrum all instruments on B line', 'Interpreter', 'Latex', 'FontSize', 16)


%% Timeseries from all instruments, everywhere

% % %
% % newFigDims([14.5, 5.3])
% %     hold on
% %     for i = 1:length(list_alldata)
% %         if lgooddims(i)
% %             plot(dataspec.(list_alldata{i}).dtime, ...
% %                  dataspec.(list_alldata{i}).FluxSS, '.-')
% %         end
% %     end
% % 
% %     %
% %     grid on
% %     set(gca, 'FontSize', 16, 'Box', 'on', ...
% %              'XLim', dataspec.(list_alldata{1}).dtime([1, end]))
% %     

%% Timeseries from all instruments on B line

%
newFigDims([14.5, 5.3])
    hold on
    plot(datagridded.dtime, datagridded.Bline.FluxSS.', '.-')

    %
    grid on
    set(gca, 'FontSize', 16, 'Box', 'on')
    %
    overlayline('v', datagridded.Bline.timeavg.timelims, '--r', 'LineWidth', 2)

    %
    ylabel('Flux [W m$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', 18)
    %
    title('Depth-integrated flux on B line', 'Interpreter', 'Latex', 'FontSize', 18)


%% draft -- plot mean normalized flux (by offshore most estimate)
% as a function of x (distance from shore)
%
% (overlay mean bathymetry)

% % %
% % newFigDims([9.8, 5])
% %     %
% %     plot(datagridded.Bline.X, ...
% %          datagridded.Bline.timeavg.FluxSSfraction, ...
% %          '.-k', 'MarkerSize', 24, 'LineWidth', 2)
% %     %
% %     grid on
% %     set(gca, 'FontSize', 16, 'Box', 'on', 'XLim', [-800, 0])
% %     %
% %     overlayline('h', 1, '-k')
% % 
% %     %
% %     xlabel('X [m]', 'Interpreter', 'Latex', 'FontSize', 18)
% %     ylabel('Flux divided by offshore flux', 'Interpreter', 'Latex', 'FontSize', 18)
% %     title('B line', 'Interpreter', 'Latex', 'FontSize', 22)


%% Plot mean normalized flux (by offshore most estimate)
% as a function of x (distance from shore)
%
% (overlay mean bathymetry)

%
l_spotter = strcmp(datagridded.Bline.instrument, "spotter");
l_soloD = strcmp(datagridded.Bline.instrument, "soloD");
l_aquadopp = strcmp(datagridded.Bline.instrument, "aquadopp");
l_signature = strcmp(datagridded.Bline.instrument, "signature");

%
newFigDims([10.7, 6.85])
    %
    axes('Position', [0.15, 0.15, 0.7, 0.6]);
    hold on
    %
    yyaxis left
    %
    plot(datagridded.Bline.X, ...
         datagridded.Bline.timeavg.FluxSSfraction, ...
         '-b', 'LineWidth', 2)
    %
    plot(datagridded.Bline.X(l_spotter), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_spotter), ...
         '.k', 'MarkerSize', 60)
    plot(datagridded.Bline.X(l_soloD), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_soloD), ...
         'sb', 'MarkerSize', 24, 'MarkerFaceColor', 'b')
    plot(datagridded.Bline.X(l_aquadopp), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_aquadopp), ...
         'pr', 'MarkerSize', 32, 'MarkerFaceColor', 'r')
    plot(datagridded.Bline.X(l_signature), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_signature), ...
         'dg', 'MarkerSize', 20, 'MarkerFaceColor', 'g')


    %
    grid on
    set(gca, 'FontSize', 16, 'Box', 'on')
    overlayline('h', 1, '-k')
    %
    xlabel('X [m]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Flux divided by offshore flux', 'Interpreter', 'Latex', 'FontSize', 18)

    %
    yyaxis right
    %
    plot(datagridded.Bline.X, -datagridded.Bline.meandepth, '.-k', 'LineWidth', 2) 
    ylabel('bottom depth [m]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylim([-35, 5])

    %
    haxs = gca;
    haxs.YAxis(1).Color = 'b';
    haxs.YAxis(2).Color = 'k';
    %
    haxs.YAxis(1).TickValues = 0:0.25:2;
    haxs.YAxis(2).TickValues = -35:5:5;

    %
    xlim([-850, 0])
%     haxs.XAxis(2).Limits = [-850, 0];
    % With B01 spotter
    haxs.YAxis(1).Limits = [0, 1.7];
    haxs.YAxis(2).Limits = [-25, 5];

% %     % Without B01 spotter
% %     haxs.YAxis(1).Limits = [0, 1.7];
% %     haxs.YAxis(2).Limits = [-35, 3.2];

    %
    title({['Energy flux decay on B line (flux at offshore'];['boundary = ' num2str(datagridded.Bline.timeavg.FluxSS(1)/1e3, '%.1f') ' kW m$^{-1}$, $H_{sig} =$ ' num2str(datagridded.Bline.timeavg.HsigSS(1), '%.1f') ' m)']}, ...
          'Interpreter', 'Latex', 'FontSize', 22)
    

%% With legend

%
l_spotter = strcmp(datagridded.Bline.instrument, "spotter");
l_soloD = strcmp(datagridded.Bline.instrument, "soloD");
l_aquadopp = strcmp(datagridded.Bline.instrument, "aquadopp");
l_signature = strcmp(datagridded.Bline.instrument, "signature");

%
newFigDims([10.7, 6.85])
    %
    axes('Position', [0.1, 0.15, 0.625, 0.6]);
    hold on
    %
    yyaxis left
    %
    plot(datagridded.Bline.X, ...
         datagridded.Bline.timeavg.FluxSSfraction, ...
         '-b', 'LineWidth', 2)
    %
    hp_Spotter = plot(datagridded.Bline.X(l_spotter), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_spotter), ...
         '.k', 'MarkerSize', 60);
    hp_SoloD = plot(datagridded.Bline.X(l_soloD), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_soloD), ...
         'sb', 'MarkerSize', 24, 'MarkerFaceColor', 'b');
    hp_Aquadopp = plot(datagridded.Bline.X(l_aquadopp), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_aquadopp), ...
         'pr', 'MarkerSize', 32, 'MarkerFaceColor', 'r');
    hp_Sig = plot(datagridded.Bline.X(l_signature), ...
         datagridded.Bline.timeavg.FluxSSfraction(l_signature), ...
         'dg', 'MarkerSize', 20, 'MarkerFaceColor', 'g');


    %
    grid on
    set(gca, 'FontSize', 16, 'Box', 'on')
    overlayline('h', 1, '-k')
    %
    xlabel('X [m]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Flux divided by offshore flux', 'Interpreter', 'Latex', 'FontSize', 18)

    %
    yyaxis right
    %
    plot(datagridded.Bline.X, -datagridded.Bline.meandepth, '.-k', 'LineWidth', 2) 
    ylabel('bottom depth [m]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylim([-35, 5])

    %
    haxs = gca;
    haxs.YAxis(1).Color = 'b';
    haxs.YAxis(2).Color = 'k';
    %
    haxs.YAxis(1).TickValues = 0:0.25:2;
    haxs.YAxis(2).TickValues = -35:5:5;

    %
    xlim([-900, 0])
%     haxs.XAxis(2).Limits = [-850, 0];
    % With B01 spotter
    haxs.YAxis(1).Limits = [0, 1.2];
    haxs.YAxis(2).Limits = [-20, 0];

% %     % Without B01 spotter
% %     haxs.YAxis(1).Limits = [0, 1.7];
% %     haxs.YAxis(2).Limits = [-35, 3.2];

    %
    title({['Energy flux on X line (flux at offshore'];['boundary = ' num2str(datagridded.Bline.timeavg.FluxSS(1)/1e3, '%.1f') ' kW m$^{-1}$, $H_{sig} =$ ' num2str(datagridded.Bline.timeavg.HsigSS(1), '%.1f') ' m)']}, ...
          'Interpreter', 'Latex', 'FontSize', 22)
    
    %
    hleg = legend([hp_Spotter, hp_SoloD, hp_Aquadopp, hp_Sig], ...
                  'Spotter', 'SoloD', 'Aquadopp', 'Signature', ...
                  'Location', 'EastOutside');
        hleg.Position = [0.825, 0.325, 0.15, 0.35];
        


