%% Run and energy flux calculation

clear
close all


%%

% % %One of each instrument
% % list_mooring_IDs = {'B03', 'E08', 'B10', 'B11', 'B12'};
%list_mooring_IDs = {'B10'};

% Along B line
list_mooring_IDs = {'B03', 'B04', 'B05', 'B06', 'B07', ...
                    'E08', ...
                    'B08', 'B09', 'B10', 'B11', 'B12', ...
                    'B13', 'B14', 'B15', 'B16', 'B17', 'B18'};

%
list_mooring_IDs = {'E08', 'E08pres'};
list_mooring_IDs = {'E02', 'E02pres'};

    
%%

freq_lims = [0.06, 0.1];
% freq_lims = [0.1, 0.15];


%%

%
clear waveStats

%
for i = 1:length(list_mooring_IDs)

    %
    if length(list_mooring_IDs{i})==3
        %
        waveStats.(list_mooring_IDs{i}) = wavestats_roxsi2022(list_mooring_IDs{i}(1:3), freq_lims);
    else
        %
        waveStats.(list_mooring_IDs{i}) = wavestats_roxsi2022(list_mooring_IDs{i}(1:3), freq_lims, true);
    end
end



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
for i1 = 1:length(list_mooring_IDs)
% % % 
% %     % Identify type of instrument
% %     ind_underscore_aux = strfind(list_alldata{i1}, '_');
% %     %
% %     strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));
% % 
% %     % Make bottom depth a column vector and positive
% %     dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth) = ...
% %             abs(dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth)(:));
    

    %
    bottomdepth_aux = waveStats.(list_mooring_IDs{i1}).bottomdepth;
    meanfreq_aux = 1./waveStats.(list_mooring_IDs{i1}).Tmean;
    %
    waveStats.(list_mooring_IDs{i1}).kmean = NaN(length(meanfreq_aux), 1);


    % Loop over time (i.e. over varying bottom depth)
    for i2 = 1:length(bottomdepth_aux)

        %
        if (bottomdepth_aux(i2) > 1) && ~isnan(meanfreq_aux(i2))
            %
            disp_rel_eval = @(k) disp_rel(k, meanfreq_aux(i2), bottomdepth_aux(i2));
            %
            meank_aux = fzero(disp_rel_eval, [(2*pi/(5000)), (2*pi/(1))]);

            %
            waveStats.(list_mooring_IDs{i1}).kmean(i2) = meank_aux;

        end
    end
end


%%
% ----------------------------------------
% -- COMPUTE cg, ENERGY, and ENERGY FLUX -
% ----------------------------------------

%% Compute group velocity

%
for i1 = 1:length(list_mooring_IDs)

% %     ind_underscore_aux = strfind(list_alldata{i1}, '_');
% %     strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    %
    k_aux = waveStats.(list_mooring_IDs{i1}).kmean;

    %
    waveStats.(list_mooring_IDs{i1}).kH = k_aux .* waveStats.(list_mooring_IDs{i1}).bottomdepth;

    %
    waveStats.(list_mooring_IDs{i1}).cp = sqrt(g * tanh(waveStats.(list_mooring_IDs{i1}).kH) ./ k_aux);
    
    %
    waveStats.(list_mooring_IDs{i1}).cg = waveStats.(list_mooring_IDs{i1}).cp .* ...
                                   0.5 .* (1 + (2*waveStats.(list_mooring_IDs{i1}).kH./sinh(2*waveStats.(list_mooring_IDs{i1}).kH)));

end



%% Compute depth-integrated total energy

%
for i1 = 1:length(list_mooring_IDs)
   
% %     ind_underscore_aux = strfind(list_alldata{i1}, '_');
% %     strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));

    %
    waveStats.(list_mooring_IDs{i1}).E = (1/16) * rho0 * g * (waveStats.(list_mooring_IDs{i1}).Hsig).^2;
    waveStats.(list_mooring_IDs{i1}).E = waveStats.(list_mooring_IDs{i1}).E(:);

end


%% Compute depth-integrated energy flux


%
for i1 = 1:length(list_mooring_IDs)
   

    %
    waveStats.(list_mooring_IDs{i1}).Flux = waveStats.(list_mooring_IDs{i1}).E .* ...
                                            waveStats.(list_mooring_IDs{i1}).cg;

end


%%
% -------------------------------------
% ----- DO ENERGETICS SPECTRALLY ------
% -------------------------------------

%%

% (PS: frequency in Hz and k in radians per meter)
disp_rel = @(k, freq, H) g*k*tanh(k*H) - (2*pi*freq)^2;


%% Compute wavenumber AT EACH FREQUENCY

%
tic
for i1 = 1:length(list_mooring_IDs)
    
    %
    waveStats.(list_mooring_IDs{i1}).spectrally.k = NaN(length(waveStats.(list_mooring_IDs{i1}).frequency), ...
                                                        length(waveStats.(list_mooring_IDs{i1}).dtime));
    waveStats.(list_mooring_IDs{i1}).spectrally.cg = waveStats.(list_mooring_IDs{i1}).spectrally.k;
    waveStats.(list_mooring_IDs{i1}).spectrally.E = waveStats.(list_mooring_IDs{i1}).spectrally.cg;
    waveStats.(list_mooring_IDs{i1}).spectrally.Flux = waveStats.(list_mooring_IDs{i1}).spectrally.cg;
% % % 
% %     % Identify type of instrument
% %     ind_underscore_aux = strfind(list_alldata{i1}, '_');
% %     %
% %     strinstr_aux = list_alldata{i1}(1:(ind_underscore_aux-1));
% % 
% %     % Make bottom depth a column vector and positive
% %     dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth) = ...
% %             abs(dataspec.(list_alldata{i1}).(varfieldnames.(strinstr_aux).bottomdepth)(:));
    

    %
    bottomdepth_aux = waveStats.(list_mooring_IDs{i1}).bottomdepth;
    
    %
    for i2 = 1:length(waveStats.(list_mooring_IDs{i1}).dtime)
        %
        if (bottomdepth_aux(i2) > 1)
            %
            for i3 = 1:length(waveStats.(list_mooring_IDs{i1}).frequency)
                %
                if (waveStats.(list_mooring_IDs{i1}).frequency(i3) >= 0.05) && ...
                   (waveStats.(list_mooring_IDs{i1}).frequency(i3) < 0.35)    
                    
                    %
                    disp_rel_eval = @(k) disp_rel(k, waveStats.(list_mooring_IDs{i1}).frequency(i3), bottomdepth_aux(i2));
                    %
                    k_aux = fzero(disp_rel_eval, [(2*pi/(5000)), (2*pi/(1))]);

                    %
                    waveStats.(list_mooring_IDs{i1}).spectrally.k(i3, i2) = k_aux;
        
                end
            end
        end
    end
    

    toc
    
    % Compute group velocity
    waveStats.(list_mooring_IDs{i1}).spectrally.kH = waveStats.(list_mooring_IDs{i1}).spectrally.k .* ...
                                                     repmat(waveStats.(list_mooring_IDs{i1}).bottomdepth(:).', length(waveStats.(list_mooring_IDs{i1}).frequency), 1);

    %
    waveStats.(list_mooring_IDs{i1}).spectrally.cp = sqrt(g * tanh(waveStats.(list_mooring_IDs{i1}).spectrally.kH) ./ waveStats.(list_mooring_IDs{i1}).spectrally.k);
    
    %
    waveStats.(list_mooring_IDs{i1}).spectrally.cg = waveStats.(list_mooring_IDs{i1}).spectrally.cp .* ...
                                   0.5 .* (1 + (2*waveStats.(list_mooring_IDs{i1}).spectrally.kH./sinh(2*waveStats.(list_mooring_IDs{i1}).spectrally.kH)));
                               
	% Compute Energy
    waveStats.(list_mooring_IDs{i1}).spectrally.E = rho0 * g * (waveStats.(list_mooring_IDs{i1}).See);
    %waveStats.(list_mooring_IDs{i1}).spectrally.E = waveStats.(list_mooring_IDs{i1}).spectrally.E(:);
    
    %
    linlims_aux = waveStats.(list_mooring_IDs{i1}).linlims;
    %
    waveStats.(list_mooring_IDs{i1}).spectrally.inband.freq_lims = waveStats.(list_mooring_IDs{i1}).freq_lims;
    waveStats.(list_mooring_IDs{i1}).spectrally.inband.linlims = linlims_aux;
    
    %
    waveStats.(list_mooring_IDs{i1}).spectrally.inband.Hsig = trapz(waveStats.(list_mooring_IDs{i1}).frequency(linlims_aux), ...
                                                                    waveStats.(list_mooring_IDs{i1}).See(linlims_aux, :), 1);
	waveStats.(list_mooring_IDs{i1}).spectrally.inband.Hsig = 4*sqrt(waveStats.(list_mooring_IDs{i1}).spectrally.inband.Hsig);
    
    % Compute Energy Flux in one frequency band
    waveStats.(list_mooring_IDs{i1}).spectrally.inband.Flux = trapz(waveStats.(list_mooring_IDs{i1}).frequency(linlims_aux), ...
                                                                    waveStats.(list_mooring_IDs{i1}).spectrally.E(linlims_aux, :) .* ...
                                                                    waveStats.(list_mooring_IDs{i1}).spectrally.cg(linlims_aux, :), 1);
    
	%
    waveStats.(list_mooring_IDs{i1}).spectrally.cumulativeFlux = NaN(size(waveStats.(list_mooring_IDs{i1}).spectrally.E));
	%
    waveStats.(list_mooring_IDs{i1}).spectrally.cumulativeFlux(linlims_aux, :) = cumsum(waveStats.(list_mooring_IDs{i1}).spectrally.E(linlims_aux, :) .* ...
                                                                                        waveStats.(list_mooring_IDs{i1}).spectrally.cg(linlims_aux, :), 1) .* ...
                                                                                 waveStats.(list_mooring_IDs{i1}).df;
                                                         
	%
    disp(['--- Done with spectral energy calculation for mooring ' list_mooring_IDs{i1} ' ---'])
	%
    toc
end


%%
% -------------------------------------
% ------------- TIME GRID -------------
% -------------------------------------

%% Load mooring table

%
addpath(genpath('/home/omarques/Documents/MATLAB/roxsi_largescale_dataproc/'))

%
mooringtable = load(['/home/omarques/Documents/MATLAB' ...
                     '/roxsi_largescale_dataproc/ROXSI2022_mooringtable.mat']);
mooringtable = mooringtable.mooringtable;


%%

%
datagridded.Nmoorings = length(list_mooring_IDs);

%
datagridded.mooringID = list_mooring_IDs;

%
Nmoorings = length(list_mooring_IDs);


%%

%
datagridded.instrument = strings(1, Nmoorings);

%
datagridded.latitude = NaN(Nmoorings, 1);
datagridded.longitude = NaN(Nmoorings, 1);
% % datagridded.X = NaN(Nmoorings, 1);
% % datagridded.Y = NaN(Nmoorings, 1);

%
for i1 = 1:Nmoorings

    %
    datagridded.instrument(i1) = convertCharsToStrings(waveStats.(list_mooring_IDs{i1}).instrument);

    %
    indmatch_table = 0;
    lloop_aux = true;
    %
    while (indmatch_table < size(mooringtable, 1)) && lloop_aux
        %
        indmatch_table = indmatch_table + 1;
        %
        moorID_ontable_aux = char(mooringtable.mooringID(indmatch_table));
        %
        if strcmp(moorID_ontable_aux(1:3), list_mooring_IDs{i1}(1:3))
            lloop_aux = false;
        end
    end

    %
    datagridded.latitude(i1) = mooringtable.latitude(indmatch_table);
    datagridded.longitude(i1) = mooringtable.longitude(indmatch_table);
end

% Compute (x, y) coordinate
[datagridded.X, datagridded.Y] = ROXSI_lltoxy(datagridded.latitude, datagridded.longitude);
    

%%

%
datagridded.dtime = datetime(2022, 06, 17, 20, 00, 00) : hours(1) : ...
                    datetime(2022, 07, 20, 02, 00, 00);
datagridded.dtime.TimeZone = 'America/Los_Angeles';
%
datagridded.dtime = datagridded.dtime(:);


%%

%
datagridded.meandepth = NaN(length(list_mooring_IDs), 1);
% % datagridded.cg = NaN(length(list_mooring_IDs), length(datagridded.dtime));
datagridded.Flux = NaN(length(list_mooring_IDs), length(datagridded.dtime));
datagridded.Hsig = NaN(length(list_mooring_IDs), length(datagridded.dtime));


%%

%
disp('--- Starting to time grid the data ---')
tic
%
for i = 1:length(list_mooring_IDs)
   
    %
    datagridded.meandepth(i) = mean(waveStats.(list_mooring_IDs{i}).bottomdepth, 'omitnan');

% %     %
% %     datagridded.Flux(i, :) = interp1(waveStats.(list_mooring_IDs{i}).dtime, ...
% %                                      waveStats.(list_mooring_IDs{i}).Flux, ...
% %                                      datagridded.dtime);
                                 
    %
    datagridded.Flux(i, :) = interp1(waveStats.(list_mooring_IDs{i}).dtime, ...
                                     waveStats.(list_mooring_IDs{i}).spectrally.inband.Flux, ...
                                     datagridded.dtime);
                                 
    %
    datagridded.Hsig(i, :) = interp1(waveStats.(list_mooring_IDs{i}).dtime, ...
                                     waveStats.(list_mooring_IDs{i}).spectrally.inband.Hsig, ...
                                     datagridded.dtime);

end

disp('--- Done with interpolation ---')
toc




% % %% Results along B line
% % 
% % %
% % ind_Bline = 0;
% % 
% % %
% % for i1 = 1:length(list_alldata)
% % 
% %     % Identify type of instrument
% %     ind_underscore_aux = strfind(list_alldata{i1}, '_');
% % 
% %     %
% %     if strcmp(list_alldata{i1}(ind_underscore_aux+1), 'B')
% %         ind_Bline = ind_Bline + 1;
% %     end
% % end
% % 
% % %
% % datagridded.Bline.mooringID = strings(ind_Bline, 1);
% % datagridded.Bline.instrument = strings(ind_Bline, 1);
% % %
% % datagridded.Bline.X = NaN(ind_Bline, 1);
% % datagridded.Bline.Y = NaN(ind_Bline, 1);
% % datagridded.Bline.meandepth = NaN(ind_Bline, 1);
% % %
% % datagridded.Bline.FluxSS = NaN(ind_Bline, length(datagridded.dtime));
% % 
% % 
% % %
% % ind_Bline = 0;
% % %
% % for i = 1:length(list_alldata)
% % 
% %     % Identify type of instrument
% %     ind_underscore_aux = strfind(list_alldata{i}, '_');
% % 
% %     %
% %     if strcmp(list_alldata{i}(ind_underscore_aux+1), 'B')
% % 
% %         %
% %         ind_Bline = ind_Bline + 1;
% % 
% %         %
% %         datagridded.Bline.mooringID(ind_Bline) = convertCharsToStrings(list_alldata{i}((end-2):end));
% %         instrument_char_aux = list_alldata{i}(1:(ind_underscore_aux-1));
% %         datagridded.Bline.instrument(ind_Bline) = convertCharsToStrings(instrument_char_aux);
% % 
% %         %
% %         datagridded.Bline.X(ind_Bline) = dataspec.(list_alldata{i}).X;
% %         datagridded.Bline.Y(ind_Bline) = dataspec.(list_alldata{i}).Y;
% %         datagridded.Bline.meandepth(ind_Bline) = mean(dataspec.(list_alldata{i}).(varfieldnames.(instrument_char_aux).bottomdepth), 'omitnan');
% % 
% %         %
% %         datagridded.Bline.HsigSS(ind_Bline, :) = ...
% %                             interp1(dataspec.(list_alldata{i}).dtime, ...
% %                                     dataspec.(list_alldata{i}).(varfieldnames.(instrument_char_aux).HsigSS), ...   % same name for all??
% %                                     datagridded.dtime);
% % 
% %         %
% %         datagridded.Bline.FluxSS(ind_Bline, :) = ...
% %                             interp1(dataspec.(list_alldata{i}).dtime, ...
% %                                     dataspec.(list_alldata{i}).FluxSS, ...
% %                                     datagridded.dtime);
% %     end
% % 


%%
% ------------------------------------------
% ------------------------------------------
% ------------------------------------------
% ------------------------------------------

%%

%
lnogaps = ~isnan(mean(datagridded.Flux, 1));

%%

%
figure
    hold on
    grid on
    %
    plot(datagridded.dtime, datagridded.Flux.', '.-')

    
%%

%
list_instruments = ["Spotter", "SmartMooring", "SoloD", "Signature", "Aquadopp"];

%
figure
    hold on
    %
    plot(datagridded.X, ...
         mean(datagridded.Flux(:, lnogaps), 2), '.-')
    %
    for i = 1:datagridded.Nmoorings
        %
        if strcmp(datagridded.instrument(i), "Spotter")
            %
            plot(datagridded.X(i), ...
                 mean(datagridded.Flux(i, lnogaps), 2), '.y', 'MarkerSize', 42)
        end
        %
        if strcmp(datagridded.instrument(i), "SmartMooring")
            %
            plot(datagridded.X(i), ...
                 mean(datagridded.Flux(i, lnogaps), 2), '.r', 'MarkerSize', 42)
        end
        %
        if strcmp(datagridded.instrument(i), "SoloD")
            %
            plot(datagridded.X(i), ...
                 mean(datagridded.Flux(i, lnogaps), 2), '.b', 'MarkerSize', 42)
        end
        %
        if strcmp(datagridded.instrument(i), "Aquadopp")
            %
            plot(datagridded.X(i), ...
                 mean(datagridded.Flux(i, lnogaps), 2), 'dr', 'MarkerSize', 16, 'MarkerFaceColor', 'r')
        end
        %
        if strcmp(datagridded.instrument(i), "Signature")
            %
            plot(datagridded.X(i), ...
                 mean(datagridded.Flux(i, lnogaps), 2), 'dg', 'MarkerSize', 16, 'MarkerFaceColor', 'g')
        end
    end

    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'XGrid', 'on', 'YGrid', 'on')
    xlabel('Cross-shore [m]', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Energy flux [W m$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', 18)
    
    
%% Plot summary with many panels

% % %
% % figure
% %     %
% %     plot(waveStats.E08.dtime, waveStats.E08.spectrally.inband.Hsig)
% %     hold on
% %     plot(waveStats.E08pres.dtime, waveStats.E08pres.spectrally.inband.Hsig)
% %     
    
%%

%
figure
    %
    plot(datagridded.Hsig(1, :), datagridded.Hsig(2, :), '.')
    %
    hold on
    %
    plot([0, 1.6], [0, 1.6])
    
    
    grid on
    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis([0, 1.6, 0, 1.6])


    