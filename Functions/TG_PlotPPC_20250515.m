function [meanPPC, semPPC]  = TG_PlotPPC_20250515(ppc, colorData, varargin)
%
% TG_PlotPPC creates a bar plot for PPC values with a broken y-axis.
%
% Input:
%   ppc                 --> m * n matrix with PPC values
%                           with m PPC values and n conditions
%   labels              --> 1 * n cell array with bar labels
%   colorData           --> m * 3 array with rgb color data for each bar
%
% Tim Guth, 2025

%% check for parent input

% parent axis
parentAx = [];
for v = 1:numel(varargin)
    if isa(varargin{v}, 'matlab.graphics.axis.Axes') || isa(varargin{v}, 'matlab.ui.container.Panel')
        parentAx = varargin{v};
    end
end

%% compute mean and standard error of the mean

% mean
meanPPC         = mean(ppc, 1, 'omitnan');

% standard error of the mean
semPPC          = std(ppc, 1, 'omitnan') ./ sqrt(size(ppc, 1) - sum(isnan(ppc)));

%% define axes

% axes creation
if ~isempty(parentAx)

    % plot inside given parent using relative axes
    parentAxPos     = get(parentAx, 'Position');

    % first and second axis
    if any(strcmp(varargin, 'breakNegative'))

        % first axis
        ax1Pos          = parentAxPos;
        ax1Pos(2)       = ax1Pos(2) + (ax1Pos(4) * 0.3);
        ax1Pos(4)       = ax1Pos(4) * 0.4;
        ax1             = axes('units', 'normalized', 'Position', ax1Pos);

        % second axis
        ax2Pos          = parentAxPos;
        ax2Pos(2)       = ax2Pos(2) + (ax2Pos(4) * 0.7);
        ax2Pos(4)       = ax2Pos(4) * 0.2;
        ax2             = axes('units', 'normalized', 'Position', ax2Pos);
    else

        % first axis
        ax1Pos          = parentAxPos;
        ax1Pos(2)       = ax1Pos(2) + (ax1Pos(4) * 0.1);
        ax1Pos(4)       = ax1Pos(4) * 0.6;
        ax1             = axes('units', 'normalized', 'Position', ax1Pos);

        % second axis
        ax2Pos          = parentAxPos;
        ax2Pos(2)       = ax2Pos(2) + (ax2Pos(4) * 0.7);
        ax2Pos(4)       = ax2Pos(4) * 0.2;
        ax2             = axes('units', 'normalized', 'Position', ax2Pos);
    end
    hold(ax1, 'on');
    hold(ax2, 'on');

    % third axis
    if any(strcmp(varargin, 'breakNegative'))
        ax3Pos          = parentAxPos;
        ax3Pos(2)       = ax3Pos(2) + (ax3Pos(4) * 0.1);
        ax3Pos(4)       = ax3Pos(4) * 0.2;
        ax3             = axes('units', 'normalized', 'Position', ax3Pos);
        hold(ax3, 'on');
    end
else

    % first and second axis
    if any(strcmp(varargin, 'breakNegative'))
        ax1             = axes('units', 'normalized', 'Position', [0.2, 0.3, 0.7, 0.4]);
        ax2             = axes('units', 'normalized', 'Position', [0.2, 0.7, 0.7, 0.2]);
    else
        ax1             = axes('units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.6]);
        ax2             = axes('units', 'normalized', 'Position', [0.2, 0.7, 0.7, 0.2]);
    end
    hold(ax1, 'on');
    hold(ax2, 'on');

    % third axis
    if any(strcmp(varargin, 'breakNegative'))
        ax3             = axes('units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.2]);
        hold(ax3, 'on');
    end
end

%% plot lines
yline(ax1, 0.1, ':');

%% plots bars, distributions and error bars

% bar plot
nBars           = size(ppc, 2);
b               = bar(ax1, 1:nBars, meanPPC, 'FaceColor', 'flat');

if any(strcmp(varargin, 'breakNegative'))
    b.BaseLine.LineStyle    = ':';
end

% loop through bars
for iBar = 1:nBars

    % adjust color
    b.CData(iBar, :)            = colorData(iBar, :);

    % plot distributions on first axis
    distrPlot1                  = swarmchart(ax1, repmat(iBar, [size(ppc, 1), 1]), ...
        ppc(:, iBar), ...
        'k', 'filled', ...
        'MarkerFaceColor', [1, 1, 1], ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', colorData(iBar, :), ...
        'MarkerEdgeAlpha', 0.5);
    rng(1);
    distrPlot1.SizeData         = 5;
    distrPlot1.XJitter          = 'rand';
    distrPlot1.XJitterWidth     = 0.5;

    % plot distributions on second axis
    distrPlot2                  = swarmchart(ax2, repmat(iBar, [size(ppc, 1), 1]), ...
        ppc(:, iBar), ...
        'k', 'filled', ...
        'MarkerFaceColor', [1, 1, 1], ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', colorData(iBar, :), ...
        'MarkerEdgeAlpha', 0.5);
    rng(1);
    distrPlot2.SizeData         = 5;
    distrPlot2.XJitter          = 'rand';
    distrPlot2.XJitterWidth     = 0.5;

    % plot distributions on third axis
    if any(strcmp(varargin, 'breakNegative'))
        distrPlot3                  = swarmchart(ax3, repmat(iBar, [size(ppc, 1), 1]), ...
            ppc(:, iBar), ...
            'k', 'filled', ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', colorData(iBar, :), ...
            'MarkerEdgeAlpha', 0.5);
        rng(1);
        distrPlot3.SizeData         = 5;
        distrPlot3.XJitter          = 'rand';
        distrPlot3.XJitterWidth     = 0.5;
    end
end

% % plot error bars
% errorbar(ax1, meanPPC, semPPC, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k');

% first axis settings
if any(strcmp(varargin, 'breakNegative'))
    limAx1          = [0, 0.1];
else
    limAx1          = [floor(min(ppc, [], 'all') * 100) / 100, 0.1];
end
set(ax1, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'XLim', [0.4, nBars + 0.6], 'YLim', limAx1, ...
    'YTick', [0, 0.1]);
ylabel(ax1, 'PPC');

% second axis settings
tickValuesAx2   = 0.1:0.1:1;
set(ax2, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [0.1, 1], ...
    'YTick', tickValuesAx2, ...
    'YTickLabel', cat(2, repmat({''}, 1, numel(tickValuesAx2) - 1), {'1'}));

% third axis settings
if any(strcmp(varargin, 'breakNegative'))
    tickValuesAx3   = -1:0.1:0;
    set(ax3, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [-1, 0], ...
    'YTick', tickValuesAx3, ...
    'YTickLabel', cat(2, {'-1'}, repmat({''}, 1, numel(tickValuesAx3) - 1)));
end

% link axes
if any(strcmp(varargin, 'breakNegative'))
    linkaxes([ax1, ax2, ax3], 'x');
else
    linkaxes([ax1, ax2], 'x');
end

end