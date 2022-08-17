function fighandle = Aquadopp_pcolor_lvl_1(aquadoppL1, strplt)
%% fighandle = AQUADOPP_PCOLOR_LVL_1(aquadoppL1, strplt)
%
%   inputs
%       - aquadoppL1: Aquadopp level 1 ROXSI data structure.
%       - strplt (optional): string to indicate whether averaged
%                            or full velocity data is plotted
%                            (default is averaged fields).
%
%   outputs
%       - fighandle: figure handle.
%
%
% AQUADOPP_PCOLOR_LVL_1.m does some basic pcolor plots of the
% level 1 Aquadopp data.
%
%
% Olavo Badaro Marques, 16/Aug/2022.


%% Add cmocean in the Matlab path

%
addpath(genpath(fullfile(repo_dirpath(), 'cmocean')))


%% If second input not given, set default to plot averaged quantities

if ~exist('strplt', 'var')
    %
    strplt = 'averaged';
end


%% Define dynamic fields that allow for structures within structures

%
if strcmp(strplt, 'averaged')

    %
    time_dynfld = {'averaged', 'dtime'};
    %
    v1_dynfld = {'averaged', 'Ue'};
    v2_dynfld = {'averaged', 'Vn'};
    v3_dynfld = {'averaged', 'Wup'};
    %
    a1_dynfld = {'averaged', 'a1'};
    a2_dynfld = {'averaged', 'a2'};
    a3_dynfld = {'averaged', 'a3'};

%
elseif strcmp(strplt, 'full')

    %
    time_dynfld = {'dtime'};
    %
    v1_dynfld = {'Ue'};
    v2_dynfld = {'Vn'};
    v3_dynfld = {'Wup'};
    %
    a1_dynfld = {'a1'};
    a2_dynfld = {'a2'};
    a3_dynfld = {'a3'};

%
else
    %
    error('Invalid input strplt.')
end


% % % These dynamic fields (though not recommended
% % % by Mathworks), allow for
% % 
% % figure,
% %     pcolor(getfield(aquadopplvl1, v1_dynfld{:}))


%% Get time in separate variable for convenience

dtime = getfield(aquadoppL1, time_dynfld{:});


%% Make the figure

%
fighandle = figure;

    %
    haxs_v1 = axes('Position', [0.1, 0.65, 0.35, 0.25]);
    haxs_v2 = axes('Position', [0.1, 0.375, 0.35, 0.25]);
    haxs_v3 = axes('Position', [0.1, 0.1, 0.35, 0.25]);
    %
    haxs_a1 = axes('Position', [0.55, 0.65, 0.35, 0.25]);
    haxs_a2 = axes('Position', [0.55, 0.375, 0.35, 0.25]);
    haxs_a3 = axes('Position', [0.55, 0.1, 0.35, 0.25]);
    %
    haxs_all = [haxs_v1, haxs_a1, haxs_v2, haxs_a2, haxs_v3, haxs_a3];

    %
    hold(haxs_all, 'on')


        % ----------------------------------------------------
        %
        pcolor(haxs_v1, dtime, aquadoppL1.zhab, getfield(aquadoppL1, v1_dynfld{:}))
        pcolor(haxs_v2, dtime, aquadoppL1.zhab, getfield(aquadoppL1, v2_dynfld{:}))
        pcolor(haxs_v3, dtime, aquadoppL1.zhab, getfield(aquadoppL1, v3_dynfld{:}))
        %
        pcolor(haxs_a1, dtime, aquadoppL1.zhab, getfield(aquadoppL1, a1_dynfld{:}))
        pcolor(haxs_a2, dtime, aquadoppL1.zhab, getfield(aquadoppL1, a2_dynfld{:}))
        pcolor(haxs_a3, dtime, aquadoppL1.zhab, getfield(aquadoppL1, a3_dynfld{:}))


    %
    for i2 = 1:length(haxs_all)
        shading(haxs_all(i2), 'flat')
    end

    % Plot pressure (~depth)
    for i = 1:length(haxs_all)

        % Where the surface is
        plot(haxs_all(i), aquadoppL1.dtime, aquadoppL1.pressure, '-k', 'LineWidth', 2)

        % Where I should trim the data to remove sidelobe contamination
        plot(haxs_all(i), aquadoppL1.dtime, aquadoppL1.pressure - (2*aquadoppL1.binsize), '--k', 'LineWidth', 1)
    end


    % ----------------------------------------------------
    %
    % Colorbars
    for i2 = 1:length(haxs_all)
        hcb_aux(i2) = colorbar(haxs_all(i2));
        hcb_aux(i2).Position(3) = 0.02;
    end
    %
    for i2 = 1:2:length(haxs_all)
        hcb_aux(i2).Position(1) = 0.46;
    end
    for i2 = 2:2:length(haxs_all)
        hcb_aux(i2).Position(1) = 0.91;
    end

    %
    set(haxs_all, 'FontSize', 16, 'Box', 'on', ...
                  'XGrid', 'on', 'YGrid', 'on', ...
                  'YLim', [0, ceil(aquadoppL1.zhab(end))], ...
                  'Color', 0.7.*[1, 1, 1])

    % Simple fixed color limits
% %     set([haxs_v1, haxs_v2], 'CLim', 0.2.*[-1, 1])
    % Something that might look better
    lbelowmeanp = aquadoppL1.zhab < mean(aquadoppL1.pressure, 'omitnan');
    %
    U_ref_Clim = max([mean(std(aquadoppL1.averaged.Ue(lbelowmeanp, :), 0, 2, 'omitnan'), 'omitnan'), ...
                      mean(std(aquadoppL1.averaged.Vn(lbelowmeanp, :), 0, 2, 'omitnan'), 'omitnan')]);
    set([haxs_v1, haxs_v2], 'CLim', 1.25*U_ref_Clim .*[-1, 1])

    %
    set(haxs_v3, 'CLim', 0.05.*[-1, 1])
    

    %
    set([haxs_a1, haxs_a2, haxs_a3], 'CLim', [0, 200])

    %
    set(gcf, 'Units', 'normalized', 'Position', [0.3, 0.3, 0.4, 0.35]);
    %
    set(haxs_all(1:(end-2)), 'XTickLabel', [])
   
    %
    ylabel(haxs_all(3), 'Height above the bottom [m]', 'Interpreter', 'Latex', 'FontSize', 16)
    %
    title(haxs_all(1), ['Aquadopp ' char(aquadoppL1.mooringID) ' - SN ' char(aquadoppL1.SN) ': Ue, Vn, and Wup'], 'Interpreter', 'Latex', 'FontSize', 16)
    title(haxs_all(2), ['Aquadopp ' char(aquadoppL1.mooringID) ' - SN ' char(aquadoppL1.SN) ': a1, a2, and a3'], 'Interpreter', 'Latex', 'FontSize', 16)

    % ******* Requires cmocean toolbox *******
    set(haxs_all(1:2:end), 'Colormap', cmocean('balance'))
    set(haxs_all(2:2:end), 'Colormap', cmocean('amp'))

