function fig_L1_QC_tilt = Aquadopp_scalars_QCplot(aquadoppL1, time_1, time_2)
%% fig_L1_QC_tilt = AQUADOPP_SCALARS_QCPLOT(aquadoppL1, time_1, time_2)
%
%
%
%
%
%
%
%
%
%

%% Make a diagnostic plot of pressure, heading, pitch and roll

% 
fig_L1_QC_tilt = figure;
    %
    set(fig_L1_QC_tilt, 'units', 'normalized')
    set(fig_L1_QC_tilt, 'Position', [0.2, 0.2, 0.4, 0.6])
        %
        haxs_1 = axes(fig_L1_QC_tilt, 'Position', [0.15, 0.73, 0.7, 0.17]);
        haxs_2 = axes(fig_L1_QC_tilt, 'Position', [0.15, 0.52, 0.7, 0.17]);
        haxs_3 = axes(fig_L1_QC_tilt, 'Position', [0.15, 0.31, 0.7, 0.17]);
        haxs_4 = axes(fig_L1_QC_tilt, 'Position', [0.15, 0.10, 0.7, 0.17]);
        %
        haxs_all = [haxs_1, haxs_2, haxs_3, haxs_4];
        hold(haxs_all, 'on')
        %
        plot(haxs_1, aquadoppL1.dtime, aquadoppL1.pressure, '-k')
        plot(haxs_2, aquadoppL1.dtime, aquadoppL1.heading, '-k')
        plot(haxs_3, aquadoppL1.dtime, aquadoppL1.pitch, '-k')
        plot(haxs_4, aquadoppL1.dtime, aquadoppL1.roll, '-k')

    %
    set(haxs_all, 'FontSize', 10, 'Box', 'on', ...
                  'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_all, 'XLim', aquadoppL1.dtime([1, end]) + [-hours(12); hours(12)])
    %
% %     ylim(haxs_2, [0, 360])
    ylim(haxs_3, 14.*[-1, 1])
    ylim(haxs_4, 14.*[-1, 1])


    %
    ylabel(haxs_1, '[dbar]', 'Interpreter', 'Latex', 'FontSize', 12)
    ylabel(haxs_2, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 12)
    ylabel(haxs_3, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 12)
    ylabel(haxs_4, '[degrees]', 'Interpreter', 'Latex', 'FontSize', 12)
    %
    title(haxs_1, ['ROXSI 2022: Aquadopp ' char(aquadoppL1.mooringID) ' - SN ' ...
                   char(aquadoppL1.SN) ': pressure, heading, pitch, and roll'], ...
                  'Interpreter', 'Latex', 'FontSize', 14)
    %
    linkaxes([haxs_1, haxs_2, haxs_3, haxs_4], 'x')


    %
    for i2 = 1:length(haxs_all)
        ylims_aux = ylim(haxs_all(i2));
        %
        plot(haxs_all(i2), [time_1, time_1], ylims_aux, '--r')
        plot(haxs_all(i2), [time_2, time_2], ylims_aux, '--r')
        %
        ylim(haxs_all(i2), ylims_aux)
    end