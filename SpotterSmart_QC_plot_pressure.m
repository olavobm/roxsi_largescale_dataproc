function hfig = SpotterSmart_QC_plot_pressure(time_datenum, pressuredata, time_lims, moorID, SN)
%% hfig = SPOTTERSMART_QC_PLOT_PRESSURE(time_datenum, pressuredata, time_lims, moorID, SN)
%
%   inputs
%       -
%       -
%       -
%       -
%       -
%
%   outputs
%       -
%
%


%%

%
inds_segment = find(time_datenum >= time_lims(1), 1, 'first') : 1 : ...
               find(time_datenum <= time_lims(2), 1, 'last');

%
inds_diff_segment = inds_segment;
if inds_segment(end) == length(time_datenum)
    inds_diff_segment = inds_diff_segment(1:end-1);
end

%
timediff_inseconds = 24*3600*diff(time_datenum);


%%

%
time_datetime = datetime(time_datenum, 'ConvertFrom', 'datenum');

%
lims_plt_datetime = datetime(time_lims, 'ConvertFrom', 'datenum');
lims_plt_pressure = [min(pressuredata(inds_segment)), max(pressuredata(inds_segment))];


%%

% % inds_timeseries = 
% % inds_timediff = (inds_timeseries(1:end-1) + inds_timeseries(2:end))./2;


%% Make the QC figure

%
hfig = figure;
    %
    set(hfig, 'Units', 'normalized')
    set(hfig, 'Position', [0.2, 0.1722, 0.3387, 0.3368])
    %
    haxs_1 = axes('Position', [0.125, 0.725, 0.325, 0.18]);
    haxs_2 = axes('Position', [0.125, 0.445, 0.325, 0.18]);
    haxs_3 = axes('Position', [0.125, 0.12, 0.325, 0.18]);
    %
    haxs_4 = axes('Position', [0.54, 0.25, 0.42, 0.5]);
    %
    hold(haxs_1, 'on')
    hold(haxs_2, 'on')
    hold(haxs_3, 'on')
    hold(haxs_4, 'on')

        %
        plot(haxs_1, time_datetime, '.-k')
        plot(haxs_2, timediff_inseconds, '.-k')
        plot(haxs_3, pressuredata, '.-k')
        %
        plot(haxs_4, time_datetime, pressuredata, '.-k')


    %
    set([haxs_1, haxs_2, haxs_3, haxs_4], ...
                        'FontSize', 16, 'Box', 'on', ...
                        'XGrid', 'on', 'YGrid', 'on')
% %     %
% %     set([haxs_1, haxs_2, haxs_3], 'XLim', [-0.05*length(time_datenum), ...
% %                                            1.05.*length(time_datenum)])
%
    set([haxs_1, haxs_2, haxs_3], 'XLim', inds_segment([1, end]))
    
    %
    set(haxs_1, 'YLim', lims_plt_datetime)
    %
% %     set(haxs_2, 'YLim', [min(timediff_inseconds(inds_diff_segment)), max(timediff_inseconds(inds_diff_segment))])
    ylims_axs2 = get(haxs_2, 'YLim');
    if max(abs(ylims_axs2))>100    % Only set y axis limits of haxs_2 if there is a long gap/unrealistic time difference
        set(haxs_2, 'YLim', [-2, 4])
    end
    %
    set(haxs_3, 'YLim', lims_plt_pressure)

    %
    set(haxs_4, 'XLim', lims_plt_datetime, 'YLim', lims_plt_pressure)

    %
    xlabel(haxs_3, 'Indices of each Spotter', 'Interpreter', 'Latex', 'FontSize', 16)
    xlabel(haxs_4, 'Time [PDT]', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    ylabel(haxs_1, 'Time [PDT]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_2, 'diff(time) [s]', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_3, 'pressure', 'Interpreter', 'Latex', 'FontSize', 16)
    ylabel(haxs_4, 'pressure', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_1, 'Timestamps', 'Interpreter', 'Latex', 'FontSize', 16)
    title(haxs_2, 'diff(time)', 'Interpreter', 'Latex', 'FontSize', 16)
    title(haxs_3, 'pressure', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_4, {['ROXSI 2022: ' moorID ' - SN ' SN];'pressure timeseries'}, ...
                   'Interpreter', 'Latex', 'FontSize', 24)

    %
    linkaxes([haxs_1, haxs_2, haxs_3], 'x')
