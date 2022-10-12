function L2out = SmartMooring_proc_lvl_2(spotsmartL1, dtspec)
%% L2out = SMARTMOORING_PROC_LVL_2(spotsmartL1, dtspec)
%
%   inputs
%       - spotsmartL1: Smart Mooring level 1 data structure.
%       - dtspec:
%
%   outputs
%       - L2out: Smart Mooring level 2 data structure.
%
%
% SMARTMOORING_PROC_LVL_2.m
%
%
%
%
%
%


%% Copy some metadata from L1 to L2

%
L2out.mooringID = spotsmartL1.mooringID;
L2out.SN = spotsmartL1.SN;

%
mooringLoc = ROXSI_mooringlocation('E01sp');

L2out.site = mooringLoc.roxsiarray;

%
L2out.latitude = spotsmartL1.latitude;
L2out.longitude = spotsmartL1.longitude;


%% Copute x/y coordinates

[L2out.X, L2out.Y] = ROXSI_lltoxy(L2out.latitude, L2out.longitude, L2out.site);

%% Height of the sensor above the bottom

%
L2out.zhab = spotsmartL1.zhab;


%% Parameters for computing spectrum

% Both in seconds
windowfft = 6*60;
windowavg = 60*60;

%
timespeclims = [datetime(2022, 06, 17, 19, 0, 0), ...
                datetime(2022, 07, 20, 05, 0, 0)];

%% Frequency limits of wave bands of interest

freq_bands.IG = [0.005, 0.03];
freq_bands.SS = [0.045, 0.3];


%%
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------



%% Computes pressure spectra

%
tic
[Spp, timespec, freqvec, dof, avgpres] = ...
            spectra_scalar_reg(spotsmart.dtime, spotsmart.pressure, ...
                               windowfft, windowavg, timespeclims);
toc


%% Computes SSH spectra (or do it inside calculation above so that
% the average pressure/depth is used at windowfft instead of windowavg?????).
% Falk's code for SoloDs uses the hourly averaged pressure. The other option
% would be to write a function similar to the code above, but include
% the SSH conversion inside it

%
g = 9.8;
rho0 = 1030;

% Computes bottom depth from average pressure (which is
% hydrostatic over 1 hour)
avgbottomdepth = 1e4*avgpres ./ (rho0*g);

% Convert frequencies to wavenumbers using
% linear wave theory
%
% Function handle to compute wavenumbers
% (PS: frequency in Hz and k in radians per meter)
disp_rel = @(k, H, freq) g*k*tanh(k*H) - (2*pi*freq)^2;


%Compute wavenumbers (in radians per meter) from linear wave theory
k_matrix = NaN(length(freqvec), length(timespec));

% Loop over time
tic
for i1 = 1:length(timespec)

    % Loop over frequencies
    for i2 = 2:length(freqvec)

        %
        disp_rel_solver = @(k) disp_rel(k, avgbottomdepth(i1), freqvec(i2));

        % In radians per meter
        k_matrix(i2, i1) = fzero(disp_rel_solver, [(2*pi/(5000)), (2*pi/(1))]);

    end
end
toc

% Compute transfer function
h_matrix = repmat(avgbottomdepth, length(freqvec), 1);
correction = cosh(k_matrix .* h_matrix) ./ ...
             cosh(k_matrix * spotsmartL1.zhab);
See = Spp .* correction.^2;


%% Now compute significant wave height
% at sea-swell and infragravity wave bands

% Frequency resolution (in Hz)
df = freqvec(2) - freqvec(1);

%
lin_SSband = ((freqvec > freq_bands.SS(1)) & (freqvec < freq_bands.SS(2)));
lin_IGband = ((freqvec > freq_bands.IG(1)) & (freqvec < freq_bands.IG(2)));

%
Hsig_SS = 4*sqrt(sum(See(lin_SSband, :), 1).*df);
Hsig_IG = 4*sqrt(sum(See(lin_IGband, :), 1).*df);

% Plot pcolor of SSH spectra up until upper band
% of sea-swell wave, to check the spectra are not
% blowing up
ltrimspec = (freqvec <= freq_bands.SS(2));

%
figure
    %
    pcolor(timespec, freqvec(ltrimspec), log10(See(ltrimspec, :)))
    shading flat
    %
    hold on

    %
    plot(timespec([1, end]), freq_bands.SS(1).*[1, 1], '--r')
    plot(timespec([1, end]), freq_bands.SS(2).*[1, 1], '--r')
    plot(timespec([1, end]), freq_bands.IG(1).*[1, 1], '--k')
    plot(timespec([1, end]), freq_bands.IG(2).*[1, 1], '--k')

    %
    ylim_aux = ylim;
    ylim([ylim_aux(1), (freq_bands.SS(2) + 3*df)])


    %
    colorbar

%
figure
    %
    plot(freqvec, mean(See, 2), '-k')
    hold on
    plot(freqvec(ltrimspec), mean(See(ltrimspec, :), 2), '.-b')

    %
    set(gca, 'FontSize', 16, 'Box', 'on', ...
             'XGrid', 'on', 'YGrid', 'on', ...
             'XScale', 'log', 'YScale', 'log')
    %
    set(gca, 'XLim', [(df/1.5), (freq_bands.SS(2) + 10*df)])
    set(gca, 'YLim', [1e-3, 1.5])

    %
    ylim_aux = ylim;
    %
    plot(freq_bands.SS(1).*[1, 1], ylim_aux, '--r')
    plot(freq_bands.SS(2).*[1, 1], ylim_aux, '--r')
    plot(freq_bands.IG(1).*[1, 1], ylim_aux, '--k')
    plot(freq_bands.IG(2).*[1, 1], ylim_aux, '--k')


% -------------------------------------------------
% Combine the 2 plots above in the same figure


%% Compute mean and peak frequency in sea-swell band

%
ind_SSband = find(lin_SSband);

% 
[~, ind_peak] = max(See(lin_SSband, :), [], 1);
% Go from indices of See(lin_SSband, :) to indices of See
ind_peak = ind_SSband(ind_peak);

%
peakfreqSS = freqvec(ind_peak);

%
freq_matrix = repmat(freqvec(:), 1, length(timespec));
%
meanfreqSS = sum(See(lin_SSband, :) .* freq_matrix(lin_SSband, :), 1) ./ ...
             sum(See(lin_SSband, :), 1);


%% Add variables to L2 data structure

%
L2out.dtime = timespec;
L2out.frequency = freqvec;

%
L2out.mean_depth = avgbottomdepth;

%
% % L2out.nptsfft = windowfft / samplingfreq

%
L2out.Spp = Spp;
L2out.See = See;

%
L2out.DOF = dof(:);

%
L2out.freqbandSS = freq_bands.SS;
L2out.freqbandIG = freq_bands.IG;

%
L2out.HsigSS = Hsig_SS(:);
L2out.HsigIG = Hsig_IG(:);

%
L2out.meanfreqSS = meanfreqSS(:);
L2out.peakfreqSS = peakfreqSS(:);


%% Add README

%
time_dataproc = datetime('now', 'TimeZone', 'Local');
time_dataproc_char = datestr(time_dataproc, 'yyyy/mm/dd HH:MM:SS');

%
L2out.README = ['ROXSI Level 2 Smart Mooring data. Serial number ' ...
                '(SN) is from the Spotter buoy of the smart mooring ' ...
                'deployed at mooringID (it''s not  the SN from the ' ...
                'pressure sensor). X and Y are cartesian coordinates ' ...
                'in a local coordinate system of the deployment site ' ...
                'Data contains hourly estimates of mean_depth (m), ' ...
                'pressure spectra (dbar^2 / Hz), and sea surface ' ...
                'elevation spectra (m^2 / Hz). A spectrum is computed ' ...
                'every ' num2str(windowfft) ' seconds (with 50% ' ...
                'overlap) and averaged over an hour. Significant ' ...
                'wave height (m) is given for sea-swell (SS) and ' ...
                'infragravity (IG) bands ' ...
                '(and the frequency bands for each are in Hz). Mean ' ...
                'and peak frequencies are also computed for the SS ' ...
                'band. DOF is the number of degrees of freedom for ' ...
                'the spectral estimates (smaller DOF indicates ' ...
                'gaps in the original pressure data). Data processed ' ...
                'by ' mfilename() '.m on ' time_dataproc_char ' (TimeZone ' ...
                time_dataproc.TimeZone ').'];








