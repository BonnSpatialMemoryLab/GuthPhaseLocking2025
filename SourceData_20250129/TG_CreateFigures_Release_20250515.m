%==========================================================================
% This script recreates the main figures of the manuscript.
%
% Tim Guth, 2025
%==========================================================================

% start
clc; close all; clear;

% add functions
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\GuthPhaseLocking2025\Functions'));
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% set path
paths       = struct();
paths.data  = 'C:\Sciebo\GitCode\NeuroGuth\GuthPhaseLocking2025\SourceData_20250129';
paths.save  = fullfile(paths.data, 'Figures');
if ~isfolder(paths.save)
   mkdir(paths.save);
end

%%=========================================================================
% Figure 1
%%=========================================================================

%% Figure 1c_left

% load data
thisData    = load(fullfile(paths.data, 'Fig_1c_left.mat'));

% histogram of sesssion-wise memory performance
f           = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.7, 3.8, 4]);
objH        = histogram(thisData.avgObjRecallPerf, 0:0.05:1, 'FaceColor', [0.7, 0.7, 0.7]);
tmpAx       = get(gca);
set(gca, 'xlim', [-0.05, 1.05], 'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
xl          = xlabel('Memory performance');
yl          = ylabel('# sessions', 'units', 'normalized', 'position', [-0.05, 0.5, 0]);
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.5);
saveas(f, fullfile(paths.save, 'Fig_1c_left.png'));

%% Figure 1c_right

% load data
thisData    = load(fullfile(paths.data, 'Fig_1c_right.mat'));

% memory performance first versus second half
f           = figure('units', 'centimeters', 'position', [5, 5, 7, 7]);
axes('units', 'centimeters', 'position', [2, 2, 4, 4]);
hold on;
for iSess = 1:size(thisData.objFirstSecond, 1)
    plot([1, 2], thisData.objFirstSecond(iSess, :), '-', 'Color', [0.5, 0.5, 0.5]);
end
plot([1, 2], mean(thisData.objFirstSecond), '-', 'Color', 'blue', 'LineWidth', 3); % average across sessions
hold off;
set(gca, ...
    'xlim', [0.9, 2.1], 'xtick', [1, 2], 'xticklabel', {'First half', 'Second half'}, ...
    'ylim', [0, 1], 'ytick', [0, 1], 'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
xlabel('Trial', ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
ylabel('Memory performance', ...
    'units', 'normalized', 'position', [-0.15, 0.5, 0], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);

% t-test first vs second half
[~, objP, ~, objT] = ttest(thisData.objFirstSecond(:, 1), thisData.objFirstSecond(:, 2));
title(strcat('{P = }', num2str(objP, 3)));
saveas(f, fullfile(paths.save, 'Fig_1c_right.png'));

%% Figure 1e_left

% load data
thisData    = load(fullfile(paths.data, 'Fig_1e_left.mat'));

% histogram of trial-wise memory performance
f           = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.7, 3.8, 4]);
hold on;
locH = histogram(thisData.trialLocRecallPerf, 0:0.05:1, ...
    'FaceColor', [0.7, 0.7, 0.7]);
plot([0.5, 0.5], [0, max(locH.Values)], ':', 'Color', [1 0 0], 'LineWidth', 2);
hold off;
tmpAx       = get(gca);
set(gca, 'xlim', [-0.05, 1.05], 'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
xl          = xlabel('Memory performance');
yl          = ylabel('# recalls', 'units', 'normalized', 'position', [-0.05, 0.5, 0]);
set(gca, 'ylim', [0, max(locH.Values)], 'ytick', [0, max(locH.Values)]);
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.5);
saveas(f, fullfile(paths.save, 'Fig_1e_left.png'));

%% Figure 1e_right

% load data
thisData    = load(fullfile(paths.data, 'Fig_1e_right.mat'));

% memory performance first versus second half
f           = figure('units', 'centimeters', 'position', [5, 5, 7, 7]);
axes('units', 'centimeters', 'position', [2, 2, 4, 4]);
hold on;
for iSess = 1:size(thisData.locFirstSecond, 1)
    plot([1, 2], thisData.locFirstSecond(iSess, :), '-', 'Color', [0.5, 0.5, 0.5]);
end
plot([1, 2], mean(thisData.locFirstSecond), '-', 'Color', 'blue', 'LineWidth', 3); % average across sessions
hold off;
set(gca, ...
    'xlim', [0.9, 2.1], 'xtick', [1, 2], 'xticklabel', {'First half', 'Second half'}, ...
    'ylim', [0, 1], 'ytick', [0, 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
xlabel('Trial', ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
ylabel('Memory performance', ...
    'units', 'normalized', 'position', [-0.15, 0.5, 0], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);

% t-test first vs second half
[~, locP, ~, locT] = ttest(thisData.locFirstSecond(:, 1), thisData.locFirstSecond(:, 2));
title(strcat('{P = }', num2str(locP, 3)));
saveas(f, fullfile(paths.save, 'Fig_1e_right.png'));

%%=========================================================================
% Figure 2
%%=========================================================================

%% Figure 2a_upperRight

% load data
thisData    = load(fullfile(paths.data, 'Fig_2a_upperRight.mat'));

% filtered LFP signal
f           = figure('units', 'centimeters', 'position', [5, 5, 20, 7]);
timeInSecs  = thisData.timeInMilliseconds / 1000;
secs        = 30;
bIncluded   = timeInSecs < secs;
plot(timeInSecs(bIncluded), thisData.signalInMicrovolts(bIncluded), 'Color', 'k');
yline(thisData.thresholdInMicrovolts, 'Color', 'r');
xlabel('Time (s)');
ylabel('Voltage (µV)');
saveas(f, fullfile(paths.save, 'Fig_2a_upperRight.png'));

%% Figure 2a_lowerRight

% load data
thisData                = load(fullfile(paths.data, 'Fig_2a_lowerRight.mat'));

% convert integers to doubles
thisSpikesInMicrovolts  = double(thisData.spikesInMicrovolts);

% filtered LFP signal
f                       = figure('units', 'centimeters', 'position', [5, 5, 20, 5]);
timeAxis                = (0:size(thisSpikesInMicrovolts, 2) - 1) / thisData.samplingRate * 1000; % convert to milliseconds

% first example unit
ax1                     = subplot(1, 2, 1);
thisSpikeData           = thisSpikesInMicrovolts(thisData.clusterClass == 2, :);
plot(timeAxis, thisSpikeData, 'Color', [0.3, 0.3, 0.3]);
hold on;
plot(timeAxis, mean(thisSpikeData), 'Color', 'k');
plot(timeAxis, mean(thisSpikeData) + std(thisSpikeData), 'Color', [0.5, 0.5, 0.5]);
plot(timeAxis, mean(thisSpikeData) - std(thisSpikeData), 'Color', [0.5, 0.5, 0.5]);
xlabel('Time (ms)');
ylabel('Voltage (µV)');

% second example unit
ax2 = subplot(1, 2, 2);
thisSpikeData           = thisSpikesInMicrovolts(thisData.clusterClass == 3, :);
plot(timeAxis, thisSpikeData, 'Color', [0.3, 0.3, 0.3]);
hold on;
plot(timeAxis, mean(thisSpikeData), 'Color', 'k');
plot(timeAxis, mean(thisSpikeData) + std(thisSpikeData), 'Color', [0.5, 0.5, 0.5]);
plot(timeAxis, mean(thisSpikeData) - std(thisSpikeData), 'Color', [0.5, 0.5, 0.5]);
xlabel('Time (ms)');
ylabel('Voltage (µV)');

% link axes
linkaxes([ax1, ax2], 'y');
ylim([-200, 100]);
saveas(f, fullfile(paths.save, 'Fig_2a_lowerRight.png'));

%% Figure 2b

% load data
thisData        = load(fullfile(paths.data, 'Fig_2b.mat'));

% create figure
f               = figure('units', 'centimeters', 'position', [5, 5, 22, 8]);

% raw data
xLimit          = [-1.5, 3];
yLimit          = [-500, 200];
rectangle('Position', [0, yLimit(1, 1) + 1, 1.5, abs(yLimit(1, 1)) + yLimit(1, 2)], 'EdgeColor', 'none', 'FaceColor', [0.9, 0.9, 0.9]);
hold on;
plot(thisData.timeInSeconds, thisData.signalInMicrovolts, 'Color', 'k', 'LineWidth', 1);
xlim(xLimit);
ylim(yLimit);

% phase
cl              = cline(thisData.timeInSeconds, thisData.filteredSignalInMicrovolts, [], thisData.phaseInRadians);
set(cl, 'LineWidth', 3.5);

% colormap
hmap(1:256, 1)  = linspace(0, 1, 256);
hmap(:, [2, 3]) = 0.8; % brightness
huemap          = hsv2rgb(hmap);
colormap(huemap);

% spiketrain
spiketrain      = line([thisData.spikeTimesInSeconds'; thisData.spikeTimesInSeconds'], [-450; -400], 'Color', 'k', 'LineWidth', 1);

% labels
xlabel('Time (s)');
ylabel('Voltage (µV)');
saveas(f, fullfile(paths.save, 'Fig_2b.png'));

%% Figure 2c

% load data
thisData                = load(fullfile(paths.data, 'Fig_2c.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.025]);
saveas(f, fullfile(paths.save, 'Fig_2c.png'));

%%=========================================================================
% Figure 3
%%=========================================================================

%% Figure 3a_top

% load data
thisData                = load(fullfile(paths.data, 'Fig_3a_top.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.025]);
saveas(f, fullfile(paths.save, 'Fig_3a_top.png'));

%% Figure 3a_right

% load data
thisData                = load(fullfile(paths.data, 'Fig_3a_bottom.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.025]);
saveas(f, fullfile(paths.save, 'Fig_3a_bottom.png'));

%% Figure 3b_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_3b_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.baselineHighPowerPPC, thisData.baselineLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_3b_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineHighPowerPPC, thisData.baselineLowPowerPPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_3b_left_ppc.png'));

%% Figure 3b_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_3b_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.encodingHighPowerPPC, thisData.encodingLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_3b_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingHighPowerPPC, thisData.encodingLowPowerPPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_3b_middle_ppc.png'));

%% Figure 3b_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_3b_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.recallHighPowerPPC, thisData.recallLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_3b_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallHighPowerPPC, thisData.recallLowPowerPPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_3b_right_ppc.png'));

%%=========================================================================
% Figure 4
%%=========================================================================

%% Figure 4a

% load data
thisData        = load(fullfile(paths.data, 'Fig_4a.mat'));

% colors for baseline, encoding and recall
colorData       = [0.6, 0.4, 0.2; 0, 0.447, 0.741; 0.85, 0.325, 0.098];

% figure
percFig         = figure('units', 'centimeters', 'Position', [10, 10, 14, 10]);
b               = bar(thisData.unitsWithPhaseLocking.Row, thisData.unitsWithPhaseLocking.Variables, 'FaceColor', 'flat');
hold on;

% loop through conditions
[nCond, nReg]   = size(thisData.unitsWithPhaseLocking.Variables');
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x coordinate of the bars
    xBar(iCond, :) = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData = colorData(iCond, :);
end

% settings
set(gca, 'tickDir', 'out', 'box', 'off');
ylabel('Units with phase locking (%)'); % adjust y-axis
ylim([0, 100]);
yticks([0, 50, 100]);
saveas(percFig, fullfile(paths.save, 'Fig_4a.png'));

%% Figure 4b

% load data
thisData        = load(fullfile(paths.data, 'Fig_4bcd.mat'));

% colors for baseline, encoding and recall
colorData       = [0.6, 0.4, 0.2; 0, 0.447, 0.741; 0.85, 0.325, 0.098];

% mean PPC
uniqueRegions   = unique(thisData.data4LME.regions);
uniquePeriods   = unique(thisData.data4LME.periods);
meanPPC         = nan(size(uniqueRegions, 1), size(uniquePeriods, 1));

% loop through regions
for iBar = 1:size(uniqueRegions, 1)

    % loop through periods
    for iPer = 1:size(uniquePeriods, 1)

        % index
        bThisReg            = strcmp(cellstr(thisData.data4LME.regions), char(uniqueRegions(iBar))) & ...
            strcmp(cellstr(thisData.data4LME.periods), char(uniquePeriods(iPer)));

        % mean
        meanPPC(iBar, iPer) = mean(thisData.data4LME.ppc(bThisReg));
    end
end

% figure
ppcFig          = figure('units', 'centimeters', 'Position', [10, 10, 14, 10]);

% first, second, and third axis
ax1             = axes('units', 'normalized', 'Position', [0.2, 0.3, 0.7, 0.4]);
ax2             = axes('units', 'normalized', 'Position', [0.2, 0.7, 0.7, 0.2]);
ax3             = axes('units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.2]);

% hold on
hold(ax1, 'on');
hold(ax2, 'on');
hold(ax3, 'on');

% bar plots
b                       = bar(ax1, 1:size(uniqueRegions, 1), meanPPC, 'FaceColor', 'flat');

% loop through conditions
nCond           = size(unique(thisData.data4LME.periods), 1);
nReg            = size(unique(thisData.data4LME.regions), 1);
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x coordinate of the bars
    xBar(iCond, :)                  = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData                  = colorData(iCond, :);
    
    % adjust baseline
    b(iCond).BaseLine.LineStyle     = ':';
end

% plot distributions
for iBar = 1:size(uniqueRegions, 1)

    % loop through periods
    for iPer = 1:size(uniquePeriods, 1)

        % index
        bThisReg    = strcmp(cellstr(thisData.data4LME.regions), char(uniqueRegions(iBar))) & ...
            strcmp(cellstr(thisData.data4LME.periods), char(uniquePeriods(iPer)));

        % plot distributions on first axis
        distrPlot1                  = swarmchart(ax1, repmat(xBar(iPer, iBar), size(thisData.data4LME.ppc(bThisReg))), ...
            thisData.data4LME.ppc(bThisReg), ...
            'k', 'filled', ...
        'MarkerFaceColor', [1, 1, 1], ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', colorData(iPer, :), ...
        'MarkerEdgeAlpha', 0.5);
        distrPlot1.SizeData         = 5;
        distrPlot1.XJitter          = 'rand';
        distrPlot1.XJitterWidth     = 0.08;

        % plot distributions on second axis
        distrPlot2                  = swarmchart(ax2, repmat(xBar(iPer, iBar), size(thisData.data4LME.ppc(bThisReg))), ...
            thisData.data4LME.ppc(bThisReg), ...
            'k', 'filled', ...
        'MarkerFaceColor', [1, 1, 1], ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', colorData(iPer, :), ...
        'MarkerEdgeAlpha', 0.5);
        distrPlot2.SizeData         = 5;
        distrPlot2.XJitter          = 'rand';
        distrPlot2.XJitterWidth     = 0.08;

        % plot distributions on third axis
        distrPlot3                  = swarmchart(ax3, repmat(xBar(iPer, iBar), size(thisData.data4LME.ppc(bThisReg))), ...
            thisData.data4LME.ppc(bThisReg), ...
            'k', 'filled', ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', colorData(iPer, :), ...
            'MarkerEdgeAlpha', 0.5);
        distrPlot3.SizeData         = 5;
        distrPlot3.XJitter          = 'rand';
        distrPlot3.XJitterWidth     = 0.08;
    end
end

% first axis settings
set(ax1, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [0, 0.1], ...
    'YTick', [0, 0.1]);
ylabel(ax1, 'PPC');
yline(ax1, 0.1, ':');

% second axis settings
tickValuesAx2   = 0.2:0.1:1;
set(ax2, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [0.1, 1], ...
    'YTick', tickValuesAx2, ...
    'YTickLabel', cat(2, repmat({''}, 1, numel(tickValuesAx2) - 1), {'1'}));

% third axis settings
tickValuesAx3   = -1:0.1:0;
set(ax3, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [-1, 0], ...
    'YTick', tickValuesAx3, ...
    'YTickLabel', cat(2, {'-1'}, repmat({''}, 1, numel(tickValuesAx3) - 1)));

% link axes
linkaxes([ax1, ax2, ax3], 'x');

% save
set(ppcFig, 'renderer', 'painters');
saveas(ppcFig, fullfile(paths.save, 'Fig_4b.png'));

%% Figure 4c

% load data
thisData                = load(fullfile(paths.data, 'Fig_4bcd.mat'));

% extract regions and periods
regions                 = cellstr(thisData.data4LME.regions);
periods                 = cellstr(thisData.data4LME.periods);

% grouping variables
regionNames             = unique(regions, 'stable');

% prepare data for boxplot
groupLabels             = strcat(regions, "_", periods);
[groupLabels, sortIdx]  = sort(groupLabels);

% figure
f                       = figure('units', 'centimeters', 'Position', [10, 10, 18, 8]);
boxplot(thisData.data4LME.spikenumber(sortIdx), groupLabels, 'Colors', 'k', 'Symbol', '.');

% axes
set(gca, 'tickDir', 'out', 'box', 'off', 'YScale', 'log');
ylabel('Number of spikes');
saveas(f, fullfile(paths.save, 'Fig_4c.png'));

%% Figure 4d

% load data
thisData        = load(fullfile(paths.data, 'Fig_4bcd.mat'));

% colors for baseline, encoding and recall
colorData       = [0.6, 0.4, 0.2; 0, 0.447, 0.741; 0.85, 0.325, 0.098];

% mean PPC
uniqueRegions   = unique(thisData.data4LME.regions);
uniquePeriods   = unique(thisData.data4LME.periods);
meanPPC         = nan(size(uniqueRegions, 1), size(uniquePeriods, 1));

% loop through regions
for iBar = 1:size(uniqueRegions, 1)

    % loop through periods
    for iPer = 1:size(uniquePeriods, 1)

        % index
        bThisReg            = strcmp(cellstr(thisData.data4LME.regions), char(uniqueRegions(iBar))) & ...
            strcmp(cellstr(thisData.data4LME.periods), char(uniquePeriods(iPer)));

        % mean
        meanPPC(iBar, iPer) = mean(thisData.data4LME.spikepower(bThisReg));
    end
end

% bar plots
f               = figure('units', 'centimeters', 'Position', [10, 10, 14, 10]);
b               = bar(uniqueRegions, meanPPC, 'FaceColor', 'flat');
hold on;

% loop through conditions
nCond           = size(unique(thisData.data4LME.periods), 1);
nReg            = size(unique(thisData.data4LME.regions), 1);
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x coordinate of the bars
    xBar(iCond, :) = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData = colorData(iCond, :);
end

% plot distributions
for iBar = 1:size(uniqueRegions, 1)

    % loop through periods
    for iPer = 1:size(uniquePeriods, 1)

        % index
        bThisReg            = strcmp(cellstr(thisData.data4LME.regions), char(uniqueRegions(iBar))) & ...
            strcmp(cellstr(thisData.data4LME.periods), char(uniquePeriods(iPer)));

        % plot distributions on second axis
        distrPlot                   = swarmchart(repmat(xBar(iPer, iBar), size(thisData.data4LME.spikepower(bThisReg))), ...
            thisData.data4LME.spikepower(bThisReg), ...
            'k', 'filled', ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', colorData(iPer, :), ...
            'MarkerEdgeAlpha', 0.5);
        distrPlot.SizeData          = 5;
        distrPlot.XJitter           = 'rand';
        distrPlot.XJitterWidth      = 0.08;
    end
end

% settings
set(gca, 'tickDir', 'out', 'box', 'off');
ylabel('Mean log-power');
ylim([0, 5]);
set(f, 'renderer', 'painters');
saveas(f, fullfile(paths.save, 'Fig_4d.png'));

%%=========================================================================
% Figure 5
%%=========================================================================

%% Figure 5a

% load data
thisData        = load(fullfile(paths.data, 'Fig_5ab.mat'));

% figure
f               = figure('units', 'centimeters', 'position', [5, 5, 30, 12]);
rectangle('Position', [0, -400, 1.5, 600], 'EdgeColor', 'none', 'FaceColor', [0.9, 0.9, 0.9]);
hold on;
plot(thisData.timeInSeconds, thisData.signalInMicrovolts, 'Color', [0.588, 0, 0]);
hold on;
plot(thisData.timeInSeconds(~thisData.isOscillation), thisData.signalInMicrovolts(~thisData.isOscillation), 'Color', 'k');

% spiketrain
spiketrain      = line([thisData.spikeTimesInSeconds'; thisData.spikeTimesInSeconds'], [-350; -300], 'Color', 'k', 'LineWidth', 1);

% box settings
scaleFactor         = 100; % rescaling log-power results to make them fit to axis
boxBottom           = -10 * scaleFactor; % position of box bottom
boxHeight           = 5 * scaleFactor;
boxWidth            = 0.5;

% figure
hold on;

% plot box and slope level
for iSprint = 1:size(thisData.sprintSegTime, 2)

    % plot box
    rectangle('Position', [thisData.sprintSegTime(1, iSprint) - (boxWidth * 0.5), boxBottom, boxWidth, boxHeight], 'FaceColor', [1, 1, 1], 'EdgeColor', [0.5, 0.5, 0.5]);
end

% plot sprint results
xAxisReference      = thisData.sprintSegTime + (log10(thisData.sprintFreqs) - (log10(thisData.sprintFreqs(end)) / 2)) * (mean(diff(thisData.sprintSegTime)) / log10(thisData.sprintFreqs(end))) * 0.8;
powPlot             = plot(xAxisReference, log10(thisData.sprintSegPow) * scaleFactor + boxBottom + scaleFactor, 'Color', 'k', 'LineWidth', 2);
fooofPlot           = plot(xAxisReference, log10(thisData.sprintSegFooof) * scaleFactor + boxBottom + scaleFactor, 'Color', [0.588, 0, 0], 'LineStyle', ':', 'LineWidth', 2);
aperiodicPlot       = plot(xAxisReference, log10(thisData.sprintSegAperiodic) * scaleFactor + boxBottom + scaleFactor, 'Color', [0, 0, 0.588], 'LineWidth', 2);

% labels
set(gca, 'tickDir', 'out', 'box', 'off');
xlabel('Time (s)');
ylabel('Voltage (µV)');
saveas(f, fullfile(paths.save, 'Fig_5a.png'));

%% Figure 5b

% load data
thisData            = load(fullfile(paths.data, 'Fig_5ab.mat'));

% plot one example
f                   = figure('units', 'centimeters', 'position', [5, 5, 6, 10]);
rectangle('Position', [1, 0.1, 9, 9999.9], 'EdgeColor', 'none', 'FaceColor', [0.9, 0.9, 0.9]);
hold on;
powPlot1            = plot(thisData.sprintFreqs, thisData.sprintSegPow(:, 1), 'Color', 'k', 'LineWidth', 2);
fooofPlot1          = plot(thisData.sprintFreqs, thisData.sprintSegFooof(:, 1), 'Color', [0.588, 0, 0], 'LineStyle', ':', 'LineWidth', 2);
aperiodicPlot1      = plot(thisData.sprintFreqs, thisData.sprintSegAperiodic(:, 1), 'Color', [0, 0, 0.588], 'LineWidth', 2);
apSlope             = mean(diff(log10(thisData.sprintSegAperiodic(:, 1))) ./ diff(log10(thisData.sprintFreqs)));
sprintSegShift      = thisData.sprintSegFooof(:, 1) .* (thisData.sprintSegAperiodic(:, 1) ./ thisData.sprintSegAperiodic(end, 1));
ylim([0.1, 10000]);
set(gca, 'XScale', 'log', 'YScale', 'log', 'TickDir', 'out', 'visible', 'on');
xlabel('Frequency (Hz)');
ylabel('Power');
saveas(f, fullfile(paths.save, 'Fig_5b.png'));

%%=========================================================================
% Figure 6
%%=========================================================================

%% Figure 6b_left

% load data
thisData            = load(fullfile(paths.data, 'Fig_6b_left.mat'));

% slope
medianSlopeFr   = median([thisData.highSlopeFiringRate, thisData.lowSlopeFiringRate], 'omitnan');
f    	        = figure('units', 'centimeters', 'Position', [10, 10, 6, 12]);
parallelcoords([thisData.highSlopeFiringRate, thisData.lowSlopeFiringRate], 'Color', [0.4, 0.4, 0.4, 0.6], ...
    'LineStyle', '-', 'LineWidth', 0.1);
hold on;
parallelcoords(medianSlopeFr, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 1);
set(gca, 'yscale', 'log');
set(gca, 'XTickLabel',{'High', 'Low'});
xlim([0.75, 2.25]);
ylim([0.1, 30]);
set(gca, 'tickDir', 'out');
set(f, 'renderer', 'painters');
saveas(f, fullfile(paths.save, 'Fig_6b_left_FR.png'));

% parameters
params                  = [];
params.nSur             = 10001; % 10001
randSeed                = 444;
rng(randSeed, 'twister');
randSeedNum             = randi(100000, params.nSur, 1);

% permutation test
[slopeRank, slopeT, surSlopeT]  = TG_PermTest1D_2PS_20230727(thisData.highSlopeFiringRate, thisData.lowSlopeFiringRate, params.nSur, randSeedNum);
slopePval                       = min(slopeRank, 1 - slopeRank) * 2;

% surrogate distribution
surSlopeFigure     = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
surSlopeTCounts    = histcounts(surSlopeT, linspace(-10, 10, 201));
surSlopeHist       = histogram('BinCounts', surSlopeTCounts, 'BinEdges', linspace(-10, 10, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(slopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-10, 10]);
xticks([-10, 0, 10]);
ylim([0, 500]);
yticks([0, 500]);
set(gca,'TickDir', 'out');
box off;
saveas(surSlopeFigure, fullfile(paths.save, 'Fig_6b_left_sur.png'));

%% Figure 6b_right

% load data
thisData            = load(fullfile(paths.data, 'Fig_6b_right.mat'));

% slope
medianSlopeFr   = median([thisData.oscillationDetectedFiringRate, thisData.noOscillationDetectedFiringRate], 'omitnan');
f    	        = figure('units', 'centimeters', 'Position', [10, 10, 6, 12]);
parallelcoords([thisData.oscillationDetectedFiringRate, thisData.noOscillationDetectedFiringRate], 'Color', [0.4, 0.4, 0.4, 0.6], ...
    'LineStyle', '-', 'LineWidth', 0.1);
hold on;
parallelcoords(medianSlopeFr, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 1);
set(gca, 'yscale', 'log');
set(gca, 'XTickLabel',{'Yes', 'No'});
xlim([0.75, 2.25]);
ylim([0.1, 30]);
set(gca, 'tickDir', 'out');
set(f, 'renderer', 'painters');
saveas(f, fullfile(paths.save, 'Fig_6b_right_FR.png'));

% parameters
params                  = [];
params.nSur             = 10001;
randSeed                = 444;
rng(randSeed, 'twister');
randSeedNum             = randi(100000, params.nSur, 1);

% permutation test
[boutRank, boutT, surBoutT]     = TG_PermTest1D_2PS_20230727(thisData.oscillationDetectedFiringRate, thisData.noOscillationDetectedFiringRate, params.nSur, randSeedNum);
boutPval                        = min(boutRank, 1 - boutRank) * 2;

% surrogate distribution
surBoutFigure     = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
surBoutTCounts    = histcounts(surBoutT, linspace(-10, 10, 201));
surBoutHist       = histogram('BinCounts', surBoutTCounts, 'BinEdges', linspace(-10, 10, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(boutT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-10, 10]);
xticks([-10, 0, 10]);
ylim([0, 500]);
yticks([0, 500]);
set(gca,'TickDir', 'out');
box off;
saveas(surBoutFigure, fullfile(paths.save, 'Fig_6b_right_sur.png'));

%% Figure 6c_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_6c_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.baselineHighSlopePPC, thisData.baselineLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_6c_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineHighSlopePPC, thisData.baselineLowSlopePPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_6c_left_ppc.png'));

%% Figure 6c_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_6c_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.encodingHighSlopePPC, thisData.encodingLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_6c_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingHighSlopePPC, thisData.encodingLowSlopePPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_6c_middle_ppc.png'));

%% Figure 6c_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_6c_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.recallHighSlopePPC, thisData.recallLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_6c_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallHighSlopePPC, thisData.recallLowSlopePPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_6c_right_ppc.png'));

%% Figure 6d

% load data
thisData            = load(fullfile(paths.data, 'Fig_6d.mat'));

% parameters
params              = struct();
params.percEdges    = linspace(0, 1, 11);

% print histogram of burst percentages
prctBurstFig        = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
percCounts          = histcounts(thisData.oscillationDetected, params.percEdges);
histogram('BinCounts', percCounts, 'BinEdges', params.percEdges * 100, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1);
box off;
set(gca, 'TickDir', 'out');
xlabel('Oscillation detected (%)');
ylabel('Number of wires');
ylim([0, 400]);
title(strcat(num2str(size(thisData.oscillationDetected, 1)), 32, 'wires'));
saveas(prctBurstFig, fullfile(paths.save, 'Fig_6d.png'));

%% Figure 6e

% load data
thisData            = load(fullfile(paths.data, 'Fig_6e.mat'));

% parameters
params              = struct();
params.wireEdges    = linspace(1, 10, 10);

% print histogram of theta frequencies for all wires
wireBurstFig        = figure('units', 'centimeters', 'position', [15, 15, 6, 10]);
wireBurstCounts     = histcounts(thisData.modeFrequency, params.wireEdges);
histogram('BinCounts', wireBurstCounts, 'BinEdges', params.wireEdges, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1);
box off;
set(gca, 'TickDir', 'out');
xlabel('Frequency (Hz)');
ylabel('Number of wires');
title(strcat(num2str(size(thisData.modeFrequency, 1)), 32, 'wires'));
saveas(wireBurstFig, fullfile(paths.save, 'Fig_6e.png'));

%% Figure 6f

% load data
thisData            = load(fullfile(paths.data, 'Fig_6f.mat'));

% parameters
params              = struct();
params.binEdges     = linspace(1, 10, 19);

% print histogram of theta frequencies for all bursts
binBurstFig         = figure('units', 'centimeters', 'position', [15, 15, 6, 10]);
TG_ShadeSEM_20210714(thisData.binCenters, thisData.oscillationDetected, 'k', 0.5);
box off;
set(gca, 'TickDir', 'out');
xlabel('Frequency (Hz)');
ylabel('Probability (%)');
xticks(0:2:10);
ylim([0, 25]);
grid on;
saveas(binBurstFig, fullfile(paths.save, 'Fig_6f.png'));

%% Figure 6g_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_6g_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_6g_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_6g_left_ppc.png'));

%% Figure 6g_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_6g_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_6g_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_6g_middle_ppc.png'));

%% Figure 6g_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_6g_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_6g_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_6g_right_ppc.png'));

%%=========================================================================
% Figure 7
%%=========================================================================

%% Figure 7a

% load data
thisData                            = load(fullfile(paths.data, 'Fig_7a.mat'));
thisDataFieldnames                  = fieldnames(thisData);

% parameters
params                              = struct();
params.polarHistEdges               = linspace(-pi, pi, 21);
params.polarHistLabels              = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% loop through units
for iUnit = 1:length(thisDataFieldnames)
    
    % get data of this unit
    thisUnit    = thisData.(thisDataFieldnames{iUnit});

    % density plot
    f           = figure;
    pcolor(thisUnit.spikeTimeInSeconds, thisUnit.ybins, thisUnit.density');
    shading interp;
    xlabel('ms');
    ylabel('µV');
    title(thisDataFieldnames{iUnit});
    saveas(f, fullfile(paths.save, strcat('Fig_7a_', thisDataFieldnames{iUnit}, '_densityPlot', '.png')));

    % polar histograms
    f           = figure('units', 'centimeters', 'Position', [10, 10, 10, 4]);

    % encoding
    subplot(1, 2, 1);
    encHistc    = histcounts(thisUnit.encodingSpikePhases, params.polarHistEdges);
    encHistg    = polarhistogram('BinEdges', params.polarHistEdges, 'BinCounts', encHistc, 'FaceColor', 'k');
    hold on;
    polarplot([circ_mean(thisUnit.encodingSpikePhases), circ_mean(thisUnit.encodingSpikePhases)], [0, max(encHistc)], 'r', 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rlim([0, max(encHistc)]);
    title('Encoding');

    % retrieval
    subplot(1, 2, 2);
    recHistc    = histcounts(thisUnit.retrievalSpikePhases, params.polarHistEdges);
    recHistg    = polarhistogram('BinEdges', params.polarHistEdges, 'BinCounts', recHistc, 'FaceColor', 'k');
    hold on;
    polarplot([circ_mean(thisUnit.retrievalSpikePhases), circ_mean(thisUnit.retrievalSpikePhases)], [0, max(recHistc)], 'r', 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rlim([0, max(recHistc)]);
    title('Retrieval');

    % save
    saveas(f, fullfile(paths.save, strcat('Fig_7a_', thisDataFieldnames{iUnit}, '_polarHistograms', '.png')));

    % spike-triggered averages
    f           = figure('units', 'centimeters', 'Position', [10, 10, 16, 5]);

    % encoding
    subplot(1, 2, 1);
    plot(thisUnit.spikeTriggeredAverageTimeInSeconds, thisUnit.encodingSpikeTriggeredAverage, 'Color', 'k');
    xline(0, ':', 'Color', 'k');
    xlim([-0.5, 0.5]);
    set(gca, 'tickDir', 'out', 'box', 'off');
    xlabel('Time (s)');
    ylabel('Voltage (µV)');

    % retrieval
    subplot(1, 2, 2);
    plot(thisUnit.spikeTriggeredAverageTimeInSeconds, thisUnit.retrievalSpikeTriggeredAverage, 'Color', 'k');
    xline(0, ':', 'Color', 'k');
    xlim([-0.5, 0.5]);
    set(gca, 'tickDir', 'out', 'box', 'off');
    xlabel('Time (s)');
    ylabel('Voltage (µV)');

    % save
    saveas(f, fullfile(paths.save, strcat('Fig_7a_', thisDataFieldnames{iUnit}, '_STA', '.png')));
end

%% Figure 7b

% load data
thisData    = load(fullfile(paths.data, 'Fig_7b.mat'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encSuccPPC, thisData.encFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_7b_ppc.png'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[encRank, encT, surEncT]            = TG_PermTest1D_2PS_20230727(thisData.encSuccPPC, thisData.encFailPPC, params.nSur, randSeedNum);

% surrogate test results
encPval                             = (1 - encRank) * 2; % Bonferroni correction
if encPval > 1
    encPval = 1;
end

% plot histogram of surrogates and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surEncT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encT, 'color', 'k', 'LineWidth', 1);
xlim([-10, 10]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_7b_sur.png'));

%% Figure 7c

% load data
thisData                = load(fullfile(paths.data, 'Fig_7c.mat'));

% parameters
params                  = struct();
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure 
f                       = figure('units', 'centimeters', 'Position', [30, 10, 16, 8]);

% plot all phases for successful encoding units
subplot(1, 2, 1);
polarplot([thisData.encSuccAngle, thisData.encSuccAngle]', [zeros(size(thisData.encSuccAngle, 1), 1), thisData.encSuccMrv]', 'Color', [0.1, 0.6, 0.1, 0.6], 'LineWidth', 1);
hold on;
polarplot([circ_mean(thisData.encSuccAngle), circ_mean(thisData.encSuccAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
rticks([0, 0.5, 1]);

% plot all phases for unsuccessful encoding units
subplot(1, 2, 2);
polarplot([thisData.encFailAngle, thisData.encFailAngle]', [zeros(size(thisData.encFailAngle, 1), 1), thisData.encFailMrv]', 'Color', [0.6, 0.1, 0.1, 0.6], 'LineWidth', 1);
hold on;
polarplot([circ_mean(thisData.encFailAngle), circ_mean(thisData.encFailAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
rticks([0, 0.5, 1]);

% save
saveas(f, fullfile(paths.save, 'Fig_7c.png'));

%% Figure 7d

% load data
thisData                = load(fullfile(paths.data, 'Fig_7d.mat'));

% figure 
f                       = figure('units', 'centimeters', 'Position', [30, 10, 16, 8]);
TG_ShadeSEM_20210714(thisData.frequencies, thisData.sfcSuccessfulEncoding, [0.1, 0.6, 0.1], 0.4);
hold on;
TG_ShadeSEM_20210714(thisData.frequencies, thisData.sfcUnsuccessfulEncoding, [0.6, 0.1, 0.1], 0.4);
xlabel('Frequency (Hz)');
xlim([0, 100]);
set(gca, 'XScale', 'log', 'tickDir', 'out');
box off;
ylabel('SFC (%)');
ylim([0, 4]);
title('Encoding');
saveas(f, fullfile(paths.save, 'Fig_7d.png'));

%% Figure 7e

% load data
thisData        = load(fullfile(paths.data, 'Fig_7e.mat'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recSuccPPC, thisData.recFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_7e_ppc.png'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[recRank, recT, surRecT]            = TG_PermTest1D_2PS_20230727(thisData.recSuccPPC, thisData.recFailPPC, params.nSur, randSeedNum);

% surrogate test results
recPval                             = (1 - recRank) * 2; % Bonferroni correction
if recPval > 1
    recPval = 1;
end

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surRecT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(recT, 'color', 'k', 'LineWidth', 1);
xlim([-10, 10]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_7e_sur.png'));

%% Figure 7f

% load data
thisData                = load(fullfile(paths.data, 'Fig_7f.mat'));

% parameters
params                  = struct();
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure 
f                       = figure('units', 'centimeters', 'Position', [30, 10, 16, 8]);

% plot all phases for successful recall units
subplot(1, 2, 1);
polarplot([thisData.recSuccAngle, thisData.recSuccAngle]', [zeros(size(thisData.recSuccAngle, 1), 1), thisData.recSuccMrv]', 'Color', [0.1, 0.6, 0.1, 0.6], 'LineWidth', 1);
hold on;
polarplot([circ_mean(thisData.recSuccAngle), circ_mean(thisData.recSuccAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
rticks([0, 0.5, 1]);

% plot all phases for unsuccessful recall units
subplot(1, 2, 2);
polarplot([thisData.recFailAngle, thisData.recFailAngle]', [zeros(size(thisData.recFailAngle, 1), 1), thisData.recFailMrv]', 'Color', [0.6, 0.1, 0.1, 0.6], 'LineWidth', 1);
hold on;
polarplot([circ_mean(thisData.recFailAngle), circ_mean(thisData.recFailAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
rticks([0, 0.5, 1]);

% save
saveas(f, fullfile(paths.save, 'Fig_7f.png'));

%% Figure 7g

% load data
thisData                = load(fullfile(paths.data, 'Fig_7g.mat'));

% figure 
f                       = figure('units', 'centimeters', 'Position', [30, 10, 16, 8]);
TG_ShadeSEM_20210714(thisData.frequencies, thisData.sfcSuccessfulRetrieval, [0.1, 0.6, 0.1], 0.4);
hold on;
TG_ShadeSEM_20210714(thisData.frequencies, thisData.sfcUnsuccessfulRetrieval, [0.6, 0.1, 0.1], 0.4);
xlabel('Frequency (Hz)');
xlim([0, 100]);
set(gca, 'XScale', 'log', 'tickDir', 'out');
box off;
ylabel('SFC (%)');
ylim([0, 3]);
title('Retrieval');
saveas(f, fullfile(paths.save, 'Fig_7g.png'));

%% Figure 8

% load data
encData     = load(fullfile(paths.data, 'Fig_8a.mat'));
recData     = load(fullfile(paths.data, 'Fig_8b.mat'));

% set random seed
randSeed    = 444;
rng(randSeed, 'twister');
randSeedNum = randi(100000, 100000, 1);

% number of surrogates
params      = struct();
params.nSur = 10001;

% order
plotNames   = ...
    {'highThetaPower', 'lowThetaPower', ...
    'highAperiodicSlope', 'lowAperiodicSlope', ...
    'oscillationDetected', 'noOscillationDetected', ...
    'objectResponsiveUnits', 'nonObjectResponsiveUnits', ...
    'phaseLockingUnits', 'nonPhaseLockingUnits', ...
    'singleUnits', 'multiUnits', ...
    'objectRecall', 'locationRecall'};

% figure
f           = figure('units', 'centimeters', 'Position', [2, 0, 37, 36]);

% loop through plots
allSum      = cell(size(plotNames));
allSize     = cell(size(plotNames));
for iPlot = 1:size(plotNames, 2)

    % this encoding data
    encCond                    = encData.(plotNames{iPlot});

    % permutation test
    [encRank, encT, surEncT]    = TG_PermTest1D_2PS_20230727(encCond.encSuccPPC, encCond.encFailPPC, params.nSur, randSeedNum);
    encPval                     = (1 - encRank) * 4;
    if encPval > 1
        encPval     = 1;
    else
        encPval                 = round(encPval * 1000) / 1000;
    end

    % plot PPC
    targetAx    = axes('units', 'normalized', 'Position', ...
        [0.04 + ((iPlot - 1) / (size(plotNames, 2) + 1)), 0.75, (1 / (size(plotNames, 2) + 4)), 0.18], ...
        'Visible', 'off');
    TG_PlotPPC_20250515([encCond.encSuccPPC, encCond.encFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative', targetAx);

    % plot histogram of surrogate and empirical t-value
    axes('units', 'normalized', 'Position', [0.04 + ((iPlot - 1) / (size(plotNames, 2) + 1)), 0.68, (1 / (size(plotNames, 2) + 4)), 0.05]);
    histogram(surEncT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(encT, 'color', 'k', 'LineWidth', 1);
    xlim([-15, 15]);
    ylim([0, 1000]);
    yticks([0, 1000]);
    set(gca,'TickDir', 'out');
    box off;
    title(['P = ', num2str(encPval), '; n = ', num2str(sum(~isnan(encCond.encSuccPPC) & ~isnan((encCond.encSuccPPC))))]);
    if iPlot == 1
        xlabel(gca, 't');
        ylabel('Surrogates');
    end

    % this retrieval data
    recCond                    = recData.(plotNames{iPlot});

    % permutation test
    [recRank, recT, surRecT]    = TG_PermTest1D_2PS_20230727(recCond.recSuccPPC, recCond.recFailPPC, params.nSur, randSeedNum);
    recPval                     = (1 - recRank) * 4;
    if recPval > 1
        recPval     = 1;
    else
        recPval                 = round(recPval * 1000) / 1000;
    end

    % plot PPC
    targetAx    = axes('units', 'normalized', 'Position', ...
        [0.04 + ((iPlot - 1) / (size(plotNames, 2) + 1)), 0.45, (1 / (size(plotNames, 2) + 4)), 0.18], ...
        'Visible', 'off');
    TG_PlotPPC_20250515([recCond.recSuccPPC, recCond.recFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative', targetAx);

    % plot histogram of surrogate and empirical t-value
    axes('units', 'normalized', 'Position', [0.04 + ((iPlot - 1) / (size(plotNames, 2) + 1)), 0.38, (1 / (size(plotNames, 2) + 4)), 0.05]);
    histogram(surRecT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(recT, 'color', 'k', 'LineWidth', 1);
    xlim([-15, 15]);
    ylim([0, 1000]);
    yticks([0, 1000]);
    set(gca,'TickDir', 'out');
    box off;
    title(['P = ', num2str(recPval), '; n = ', num2str(sum(~isnan(recCond.recSuccPPC) & ~isnan((recCond.recSuccPPC))))]);
    if iPlot == 1
        xlabel(gca, 't');
        ylabel('Surrogates');
    end
end

% save
set(f, 'renderer', 'painters');
saveas(f, fullfile(paths.save, 'Fig_8.png'));

%%=========================================================================
% Figure 9
%%=========================================================================

%% Figure 9b

% load data
thisData                            = load(fullfile(paths.data, 'Fig_9b.mat'));
thisDataFieldnames                  = fieldnames(thisData);

% parameters
params                              = struct();
params.nSur                         = 10001;
params.polarHistEdges               = linspace(-pi, pi, 21);
params.polarHistLabels              = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% loop through units
for iUnit = 1:length(thisDataFieldnames)

    % get data of this unit
    thisUnit    = thisData.(thisDataFieldnames{iUnit});

    % set random seed
    rng(444, 'twister');

    % density plot
    f           = figure;
    pcolor(thisUnit.spikeTimeInSeconds, thisUnit.ybins, thisUnit.density');
    shading interp;
    xlabel('ms');
    ylabel('µV');
    title(thisDataFieldnames{iUnit});
    saveas(f, fullfile(paths.save, strcat('Fig_9b_', thisDataFieldnames{iUnit}, '_densityPlot', '.png')));

    % polar histograms
    f           = figure('units', 'centimeters', 'Position', [10, 10, 10, 4]);

    % encoding
    subplot(1, 2, 1);
    encHistc    = histcounts(cat(1, thisUnit.encodingSpikePhases{:}), params.polarHistEdges);
    encHistg    = polarhistogram('BinEdges', params.polarHistEdges, 'BinCounts', encHistc, 'FaceColor', 'k');
    hold on;
    polarplot([circ_mean(cat(1, thisUnit.encodingSpikePhases{:})), circ_mean(cat(1, thisUnit.encodingSpikePhases{:}))], [0, max(encHistc)], 'r', 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rlim([0, max(encHistc)]);
    title('Encoding');

    % retrieval
    subplot(1, 2, 2);
    recHistc    = histcounts(cat(1, thisUnit.retrievalSpikePhases{:}), params.polarHistEdges);
    recHistg    = polarhistogram('BinEdges', params.polarHistEdges, 'BinCounts', recHistc, 'FaceColor', 'k');
    hold on;
    polarplot([circ_mean(cat(1, thisUnit.retrievalSpikePhases{:})), circ_mean(cat(1, thisUnit.retrievalSpikePhases{:}))], [0, max(recHistc)], 'r', 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rlim([0, max(recHistc)]);
    title('Retrieval');

    % save
    saveas(f, fullfile(paths.save, strcat('Fig_9b_', thisDataFieldnames{iUnit}, '_polarHistograms', '.png')));

    % spike-triggered averages
    f           = figure('units', 'centimeters', 'Position', [10, 10, 16, 5]);

    % encoding
    subplot(1, 2, 1);
    plot(thisUnit.spikeTriggeredAverageTimeInSeconds, thisUnit.encodingSpikeTriggeredAverage, 'Color', 'k');
    xline(0, ':', 'Color', 'k');
    xlim([-0.5, 0.5]);
    set(gca, 'tickDir', 'out', 'box', 'off');
    xlabel('Time (s)');
    ylabel('Voltage (µV)');

    % retrieval
    subplot(1, 2, 2);
    plot(thisUnit.spikeTriggeredAverageTimeInSeconds, thisUnit.retrievalSpikeTriggeredAverage, 'Color', 'k');
    xline(0, ':', 'Color', 'k');
    xlim([-0.5, 0.5]);
    set(gca, 'tickDir', 'out', 'box', 'off');
    xlabel('Time (s)');
    ylabel('Voltage (µV)');

    % save
    saveas(f, fullfile(paths.save, strcat('Fig_9b_', thisDataFieldnames{iUnit}, '_STA', '.png')));

    %% Watson-Williams test for phase shifts

    % concatenate phases of two groups
    encRecPhases                = cat(1, thisUnit.encodingSpikePhases, thisUnit.retrievalSpikePhases);

    % class variable
    classEncRec                 = [ones(size(thisUnit.encodingSpikePhases, 1), 1); ones(size(thisUnit.retrievalSpikePhases, 1), 1) * 2];

    % empirical Watson-Williams test
    [~, wwEncRecTable]          = circ_wwtest(cat(1, encRecPhases{classEncRec == 1}), cat(1, encRecPhases{classEncRec == 2}));

    % F-values
    wwEncRec                    = wwEncRecTable{2, 5};

    % surrogate tests
    surWwEncRec                 = nan(params.nSur, 1);
    for iSur = 1:params.nSur

        % shuffle class for surrogate dataset
        surClassEncRec              = datasample(classEncRec, numel(classEncRec), 'replace', false);

        % surrogate Watson-Williams test and surrogate F-values
        [~, surWwEncRecTable]       = circ_wwtest(cat(1, encRecPhases{surClassEncRec == 1}), cat(1, encRecPhases{surClassEncRec == 2}));
        surWwEncRec(iSur, 1)        = surWwEncRecTable{2, 5};
    end

    % rank of empirical WW-test F-value in surrogate dataset
    encRecRank                  = sum(wwEncRec > surWwEncRec) / params.nSur;

    %% plot distribution of phase shifts

    % p-value
    wwP         = 1 - (sum(wwEncRec > surWwEncRec) / params.nSur);

    % correct low p-values of permutation tests (see Phipson and Smyth, 2010)
    wwP(wwP < (1 / params.nSur)) = 1 / params.nSur;

    % plot
    permFig     = figure;
    permDistr   = histogram(surWwEncRec, 'FaceColor', [0.5, 0.5, 0.5]);
    hold on;
    empVal      = xline(wwEncRec, 'k');
    title({strcat('{Permutation test (P = }', num2str(wwP, 3), ')')});
    xlabel('f-value Watson-Williams test');
    ylabel('Number of surrogates');
    ylim([0, ceil(max(permDistr.Values) / 1000) * 1000]);
    xlim([0, ceil(permDistr.BinEdges(end) / 10) * 10]);
    set(gca, 'tickDir', 'out', 'box', 'off');
    saveas(permFig, fullfile(paths.save, strcat('Fig_9b_', thisDataFieldnames{iUnit}, '_sur', '.png')));
end

%% Figure 9c

% load data
thisData                = load(fullfile(paths.data, 'Fig_9c.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiff, params.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSig, params.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_9c.png'));

%% Figure 9d

% load data
thisData                = load(fullfile(paths.data, 'Fig_9d.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiffSucc, params.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSigSucc, params.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_9d.png'));

%% Figure 9e

% load data
thisData                = load(fullfile(paths.data, 'Fig_9e.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiffFail, params.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSigFail, params.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_9e.png'));

%% Figure 9f

% load data
thisData        = load(fullfile(paths.data, 'Fig_9f.mat'));

% figure
shiftFig        = figure('units', 'centimeters', 'position', [10, 10, 18, 12]);
ax1             = axes('units', 'normalized', 'position', [0.1, 0.2, 0.6, 0.7]);

% define categories for the x-axis
shiftBounds     = categorical({'All', 'PL', 'Succ', 'Succ PL', 'Fail', 'Fail PL'});

% reorder categories
desiredOrder    = {'All', 'PL', 'Succ', 'Succ PL', 'Fail', 'Fail PL'};
shiftBounds     = reordercats(shiftBounds, desiredOrder);

% define colors for each category
shiftColors     = [
    0.16, 0.16, 0.16; ...
    0.16, 0.16, 0.16; ...
    0.1, 0.6, 0.1; ...
    0.1, 0.6, 0.1; ...
    0.6, 0.1, 0.1; ...
    0.6, 0.1, 0.1
    ];

% define transparency levels for each category
shiftAlphas     = [0.2, 0.6, 0.2, 0.6, 0.2, 0.6];

% initialize bar objects
b               = gobjects(size(shiftBounds));

% clear current axes and begin plotting
hold on;

% loop through each category to create bars with corresponding properties
for iLine = 1:numel(shiftBounds)
    b(iLine) = bar(shiftBounds(iLine), thisData.shiftPercentage(iLine));
    set(b(iLine), 'FaceColor', shiftColors(iLine, :), 'FaceAlpha', shiftAlphas(iLine));
end

% adjustments
ylabel('Ratio');
ylim([0, 15]);
yticks(0:5:15);
set(gca, 'box', 'off', 'TickDir', 'out');

% surrogate distribution
ax2             = axes('units', 'normalized', 'position', [0.75, 0.2, 0.2, 0.7]);
histogram(thisData.surrogates, 'FaceColor', [0.2, 0.2, 0.2], 'Orientation', 'horizontal');
yline(mean(thisData.surrogates));
xlim([0, 140]);
xticks([0, 140]);
ylim([0, 15]);
set(gca, 'tickDir', 'out');
yticks([]);
box off;

% chance level
ax3             = axes('units', 'normalized', 'position', [0.1, 0.2, 0.85, 0.7]);
set(ax3, 'Visible', 'off');
ylim([0, 15]);
yline(mean(thisData.surrogates), 'LineStyle', ':');
saveas(shiftFig, fullfile(paths.save, 'Fig_9f.png'));

%%=========================================================================
% Figure S1
%%=========================================================================

%% Figure S1a

% load data
thisData        = load(fullfile(paths.data, 'Fig_S1a.mat'));

% plot
unitsFig        = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
unitsPerWire    = histogram(thisData.unitsPerWire, 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('Units per wire');
ylabel('Number of wires');
set(gca, 'TickDir', 'out', 'box', 'off');
saveas(unitsFig, fullfile(paths.save, 'Fig_S1a.png'));

%% Figure S1b

% load data
thisData        = load(fullfile(paths.data, 'Fig_S1b.mat'));

% plot
isiFig          = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
lowIsiPerUnit   = histogram(thisData.percISIbelowThreeMs, 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('% ISI <3 ms');
ylabel('Number of units');
set(gca, 'TickDir', 'out', 'box', 'off');
saveas(isiFig, fullfile(paths.save, 'Fig_S1b.png'));

%% Figure S1c

% load data
thisData        = load(fullfile(paths.data, 'Fig_S1c.mat'));

% plot
frFig           = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
meanFRperUnit   = histogram(thisData.meanFR, 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('Mean FR (Hz)');
ylabel('Number of units');
set(gca, 'TickDir', 'out', 'box', 'off');
saveas(frFig, fullfile(paths.save, 'Fig_S1c.png'));

%% Figure S1d

% load data
thisData        = load(fullfile(paths.data, 'Fig_S1d.mat'));

% plot
snrFig          = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
SNRperUnit      = histogram(thisData.peakSNR, 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('Peak SNR');
ylabel('Number of units');
set(gca, 'TickDir', 'out', 'box', 'off');
saveas(snrFig, fullfile(paths.save, 'Fig_S1d.png'));

%%=========================================================================
% Figure S2
%%=========================================================================

%% Figure S2

% load data
thisData        = load(fullfile(paths.data, 'Fig_S2.mat'));

% create figure
f               = figure('units', 'centimeters', 'position', [1, 1, 18, 22]);
numAxes         = size(thisData.complexSignal, 1);

% loop through examples
for iEx = 1:numAxes
    
    % create axis
    axes('units', 'normalized', 'Position', [0.1, 0.85 - ((iEx - 1) / (numAxes + 1)), 0.65, (1 / (numAxes + 4))]);
 
    % raw data
    xLimit          = [-1.5, 3];
    yLimit          = [round(min(thisData.signalInMicrovolts{iEx, 1}) / 10) * 10 - 100, ...
        round(max(thisData.signalInMicrovolts{iEx, 1}) / 10) * 10 + 20];
    rectangle('Position', [0, yLimit(1, 1) + 1, 1.5, abs(yLimit(1, 1)) + yLimit(1, 2)], 'EdgeColor', 'none', 'FaceColor', [0.9, 0.9, 0.9]);
    hold on;
    plot(thisData.timeInSeconds{iEx, 1}, thisData.signalInMicrovolts{iEx, 1}, 'Color', 'k', 'LineWidth', 1);
    xlim(xLimit);
    ylim(yLimit);

    % phase
    cl              = cline(thisData.timeInSeconds{iEx, 1}, thisData.filteredSignalInMicrovolts{iEx, 1}, [], angle(thisData.complexSignal{iEx, 1}));
    set(cl, 'LineWidth', 3.5);

    % colormap
    hmap(1:256, 1)  = linspace(0, 1, 256);
    hmap(:, [2, 3]) = 0.8; % brightness
    huemap          = hsv2rgb(hmap);
    colormap(huemap);

    % spiketrain
    spiketrain      = line([thisData.spikeTimesInSeconds{iEx, 1}'; thisData.spikeTimesInSeconds{iEx, 1}'], [yLimit(1, 1); yLimit(1, 1) + (range(yLimit) / 8)], 'Color', 'k', 'LineWidth', 1);

    % labels
    xlabel('Time (s)');
    ylabel('Voltage (µV)');

    % create axis
    axes('units', 'normalized', 'Position', [0.85, 0.85 - ((iEx - 1) / (numAxes + 1)), 0.1, (1 / (numAxes + 4))]);

    % plot
    c2              = cline(real(thisData.complexSignal{iEx, 1}), imag(thisData.complexSignal{iEx, 1}), [], angle(thisData.complexSignal{iEx, 1}));
    set(c2, 'LineWidth', 0.5);

    % colormap
    hmap(1:256, 1)  = linspace(0, 1, 256);
    hmap(:, [2, 3]) = 0.8; % brightness
    huemap          = hsv2rgb(hmap);
    colormap(huemap);

    %  adjust figure
    axis equal;
    xline(0);
    yline(0);
    figLim          = ceil(max(abs([real(thisData.complexSignal{iEx, 1}), imag(thisData.complexSignal{iEx, 1})])) / 100) * 100;
    xlim([-figLim, figLim]);
    ylim([-figLim, figLim]);
    xticks([-figLim, figLim]);
    yticks([-figLim, figLim]);
    xlabel('Real (µV)');
    ylabel('Imaginary (µV)');
    set(gca, 'tickDir', 'out', 'box', 'off');
end

% save
saveas(f, fullfile(paths.save, 'Fig_S2.png'));

%%=========================================================================
% Figure S3
%%=========================================================================

%% Figure S3a

% load data
thisData            = load(fullfile(paths.data, 'Fig_S3abcdef.mat'));

% plot raw signal
rawSignalFig        = figure;
plot(thisData.timeInSeconds, thisData.signalInMicrovolts, 'Color', 'k', 'LineWidth', 2);thisData.signalInMicrovolts
hold on;
plot(thisData.timeInSeconds(thisData.isInterpolated), thisData.signalInMicrovolts(thisData.isInterpolated), 'Color', 'b', 'LineWidth', 2);
set(gca, 'tickDir', 'out');
ylim([-400, 400]);
yticks([-400, 0, 400]);
box off;
xlabel('Time (s)');
ylabel('Amplitude (µV)');
set(rawSignalFig, 'Renderer', 'painter');
saveas(rawSignalFig, fullfile(paths.save, 'Fig_S3a.png'));

%% Figure S3b

% load data
thisData            = load(fullfile(paths.data, 'Fig_S3abcdef.mat'));

% plot filtered signal
filteredSignalFig   = figure;
plot(thisData.timeInSeconds, thisData.filteredSignalInMicrovolts, 'Color', 'k', 'LineWidth', 2);
hold on;
plot(thisData.timeInSeconds(thisData.isInterpolated), thisData.filteredSignalInMicrovolts(thisData.isInterpolated), 'Color', 'b', 'LineWidth', 2);
set(gca, 'tickDir', 'out');
ylim([-400, 400]);
yticks([-400, 0, 400]);
box off;
xlabel('Time (s)');
ylabel('Amplitude (µV)');
set(filteredSignalFig, 'Renderer', 'painter');
saveas(filteredSignalFig, fullfile(paths.save, 'Fig_S3b.png'));

%% Figure S3c

% load data
thisData            = load(fullfile(paths.data, 'Fig_S3abcdef.mat'));

% plot hilbert signal
hilbertFig          = figure;
hilbertPlotTime     = thisData.timeInSeconds;
hilbertPlotReal     = real(thisData.hilbertComplexSignal);
hilbertPlotImag     = imag(thisData.hilbertComplexSignal);
hilbertPlotTime(thisData.isInterpolated)    = NaN;
hilbertPlotReal(thisData.isInterpolated)    = NaN;
hilbertPlotImag(thisData.isInterpolated)    = NaN;
plot3(hilbertPlotTime, hilbertPlotReal, hilbertPlotImag, 'Color', 'k', 'LineWidth', 2);
hold on;
plot3(thisData.timeInSeconds(thisData.isInterpolated), real(thisData.hilbertComplexSignal(thisData.isInterpolated)), imag(thisData.hilbertComplexSignal(thisData.isInterpolated)), 'Color', 'b', 'LineWidth', 2);
set(gca, 'tickDir', 'out');
xlabel('Time (s)');
ylabel('Real part (µV)');
zlabel('Imaginary part (µV)');
ylim([-400, 400]);
yticks([-400, 0, 400]);
zlim([-400, 400]);
zticks([-400, 0, 400]);
view(45, 45);
set(hilbertFig, 'Renderer', 'painter');
saveas(hilbertFig, fullfile(paths.save, 'Fig_S3c.png'));

%% Figure S3d

% load data
thisData            = load(fullfile(paths.data, 'Fig_S3abcdef.mat'));

% plot in polar histogram
polarPlotFig        = figure;
polarplot(angle(thisData.hilbertComplexSignal), abs(thisData.hilbertComplexSignal), 'Color', 'k', 'LineWidth', 2);
hold on;
polarplot(angle(thisData.hilbertComplexSignal(thisData.isInterpolated)), abs(thisData.hilbertComplexSignal(thisData.isInterpolated)), 'Color', 'b', 'LineWidth', 2);
set(polarPlotFig, 'Renderer', 'painter');
saveas(polarPlotFig, fullfile(paths.save, 'Fig_S3d.png'));

%% Figure S3e

% load data
thisData            = load(fullfile(paths.data, 'Fig_S3abcdef.mat'));

% plot hilbert phase estimate
hilbertPhaseFig   = figure;
plot(thisData.timeInSeconds, thisData.filteredSignalInMicrovolts, 'Color', 'k', 'LineWidth', 2);
hold on;
plot(thisData.timeInSeconds(thisData.isInterpolated), thisData.filteredSignalInMicrovolts(thisData.isInterpolated), 'Color', 'b', 'LineWidth', 2);
ylim([-400, 400]);
yticks([-400, 0, 400]);
yyaxis right;
plot(thisData.timeInSeconds, rad2deg(angle(thisData.hilbertComplexSignal)), 'Color', [0.848, 0.324, 0.098], 'LineWidth', 2);
hold on;
plot(thisData.timeInSeconds(thisData.isInterpolated), rad2deg(angle(thisData.hilbertComplexSignal(thisData.isInterpolated))), 'Color', [0.117, 0.508, 0.508], 'LineWidth', 2, 'LineStyle', '-');
set(gca, 'tickDir', 'out');
ax                  = gca;
ax.YAxis(2).Color   = [0.848, 0.324, 0.098];
ylim([-350, 350]);
yticks([-180, 0, 180]);
box off;
xlabel('Time (s)');
ylabel('Angle (degree)');
set(hilbertPhaseFig, 'Renderer', 'painter');
saveas(hilbertPhaseFig, fullfile(paths.save, 'Fig_S3e.png'));

%% Figure S3f

% load data
thisData            = load(fullfile(paths.data, 'Fig_S3abcdef.mat'));

% plot generalized phase estimate
generalizedPhaseFig   = figure;
plot(thisData.timeInSeconds, thisData.filteredSignalInMicrovolts, 'Color', 'k', 'LineWidth', 2);
hold on;
plot(thisData.timeInSeconds(thisData.isInterpolated), thisData.filteredSignalInMicrovolts(thisData.isInterpolated), 'Color', 'b', 'LineWidth', 2);
ylim([-400, 400]);
yticks([-400, 0, 400]);
yyaxis right;
plot(thisData.timeInSeconds, rad2deg(angle(thisData.generalizedPhaseComplexSignal)), 'Color', [0.848, 0.324, 0.098], 'LineWidth', 2);
hold on;
plot(thisData.timeInSeconds(thisData.isInterpolated), rad2deg(angle(thisData.generalizedPhaseComplexSignal(thisData.isInterpolated))), 'Color', [0.117, 0.508, 0.508], 'LineWidth', 2, 'LineStyle', '-');
set(gca, 'tickDir', 'out');
ax                  = gca;
ax.YAxis(2).Color   = [0.848, 0.324, 0.098];
ylim([-350, 350]);
yticks([-180, 0, 180]);
box off;
xlabel('Time (s)');
ylabel('Angle (degree)');
set(generalizedPhaseFig, 'Renderer', 'painter');
saveas(generalizedPhaseFig, fullfile(paths.save, 'Fig_S3f.png'));

%% Figure S3g

% load data
thisData            = load(fullfile(paths.data, 'Fig_S3g.mat'));

% plot mean powerspectrum
powerSpctrmFig      = figure;
TG_ShadeSEM_20210714(thisData.frequencies, thisData.powerspectra, 'k', 0.5);
xlabel('Frequency (Hz)');
ylabel('Power');
xticks(0:10);
xlim([0, 11]);
ylim([0, 0.5]);
set(gca, 'tickDir', 'out', 'box', 'off');
saveas(powerSpctrmFig, fullfile(paths.save, 'Fig_S3g.png'));

%% Figure S3h

% load data
thisData    = load(fullfile(paths.data, 'Fig_S3h.mat'));

% plot
anglesFig   = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, ...
    'FaceAlpha', thisData.angleHist.FaceAlpha);
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.02]);
saveas(anglesFig, fullfile(paths.save, 'Fig_S3h.png'));

%% Figure S3i

% load data
thisData                                    = load(fullfile(paths.data, 'Fig_S3i.mat'));

% logarithmize heatmap
logHeatMapCounts                            = log10(thisData.heatmapCounts);
logHeatMapCounts(isinf(logHeatMapCounts))   = NaN;

% plot the 2D histogram
complexFig  = figure();
pcolor(thisData.heatmapCenters, thisData.heatmapCenters, logHeatMapCounts);
shading interp;
axis equal;
box off;
set(gca, 'TickDir', 'out');
xlim([-500, 500]);
ylim([-500, 500]);
xline(0);
yline(0);
xlabel('Real (µV)');
ylabel('Imaginary (µV)');

% colorbar
c1              = colorbar;
ylabel(c1, 'log(sample count)');

% save
saveas(complexFig, fullfile(paths.save, 'Fig_S3i.png'));

%% Figure S3j

% load data
thisData                = load(fullfile(paths.data, 'Fig_S3j.mat'));

% plot
anglesCorrFig           = figure;
polarhistogram('BinEdges', thisData.angleCorrectionHist.BinEdges, ...
    'BinCounts', thisData.angleCorrectionHist.Values, ...
    'FaceColor', thisData.angleCorrectionHist.FaceColor, ...
    'FaceAlpha', thisData.angleCorrectionHist.FaceAlpha);
set(gca, 'ThetaTickLabel', thisData.angleCorrectionHist.polarHistLabels);
rlim([0, 0.4]);
saveas(anglesCorrFig, fullfile(paths.save, 'Fig_S3j.png'));

%%=========================================================================
% Figure S4
%%=========================================================================

%% Figure S4a

% load data
thisData            = load(fullfile(paths.data, 'Fig_S4a.mat'));

% figure
distrFig            = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
bD                  = bar([1, 2], [mean(thisData.successfulTaskFrequencies, 'omitnan'), mean(thisData.unsuccessfulTaskFrequencies, 'omitnan')], 'FaceColor', [0.5, 0.5, 0.5]);
hold on;
s1                  = swarmchart(ones(size(thisData.successfulTaskFrequencies)), thisData.successfulTaskFrequencies, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
    'MarkerEdgeAlpha', 0.5);
s1.SizeData          = 20;
s1.XJitter           = 'rand';
s1.XJitterWidth      = 0.5;
s2                  = swarmchart(ones(size(thisData.unsuccessfulTaskFrequencies)) * 2, thisData.unsuccessfulTaskFrequencies, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
    'MarkerEdgeAlpha', 0.5);
s2.SizeData          = 20;
s2.XJitter           = 'rand';
s2.XJitterWidth      = 0.5;
set(gca, 'tickDir', 'out', 'box', 'off', 'XTicklabel', {'Successful', 'Unsuccessful'});
ylim([0, 4]);
ylabel('Frequency (Hz)');
saveas(distrFig, fullfile(paths.save, 'Fig_S4a.png'));

%% Figure S4b

% load data
thisData            = load(fullfile(paths.data, 'Fig_S4b.mat'));

% plot
rhoSlopeCycleFig    = figure;
h1                  = histogram(thisData.correlationSlopeAndCycleFrequency, -0.3:0.01:-0.1, 'FaceColor', [0.5, 0.5, 0.5]);
set(gca, 'tickDir', 'out', 'box', 'off');
xticks([-0.3, -0.2, -0.1]);
xlabel('Spearman rho');
ylabel('Number of sessions');
saveas(rhoSlopeCycleFig, fullfile(paths.save, 'Fig_S4b.png'));

%% Figure S4c

% load data
thisData            = load(fullfile(paths.data, 'Fig_S4c.mat'));

% plot
fitSlopeCycleFig    = figure;
h2                  = histogram(thisData.slopeOfLinearFitBetweenSlopeAndCycleFrequency, -2:0.1:0, 'FaceColor', [0.5, 0.5, 0.5]);
set(gca, 'tickDir', 'out', 'box', 'off');
xticks([-2, -1, 0]);
xlabel('Slope of linear fit');
ylabel('Number of sessions');
saveas(fitSlopeCycleFig, fullfile(paths.save, 'Fig_S4c.png'));

%%=========================================================================
% Figure S6
%%=========================================================================

%% Figure S6a

% load data
thisData            = load(fullfile(paths.data, 'Fig_S6a.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.025]);
saveas(f, fullfile(paths.save, 'Fig_S6a.png'));

%% Figure S6b

% load data
thisData        = load(fullfile(paths.data, 'Fig_S6b.mat'));

% colors for baseline, encoding and recall
colorData       = [0.6, 0.4, 0.2; 0, 0.447, 0.741; 0.85, 0.325, 0.098];

% figure
percFig         = figure('units', 'centimeters', 'Position', [10, 10, 14, 10]);
b               = bar(thisData.unitsWithPhaseLocking.Row, thisData.unitsWithPhaseLocking.Variables, 'FaceColor', 'flat');
hold on;

% loop through conditions
[nCond, nReg]   = size(thisData.unitsWithPhaseLocking.Variables');
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x coordinate of the bars
    xBar(iCond, :) = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData = colorData(iCond, :);
end

% settings
set(gca, 'tickDir', 'out', 'box', 'off');
ylabel('Units with phase locking (%)'); % adjust y-axis
ylim([0, 100]);
yticks([0, 50, 100]);
saveas(percFig, fullfile(paths.save, 'Fig_S6b.png'));

%% Figure S6c

% load data
thisData    = load(fullfile(paths.data, 'Fig_S6c.mat'));

% mean PPC
uniqueFreq  = unique(thisData.ppcData.freqBand);
uniqueOsc   = unique(thisData.ppcData.oscYesNo, 'stable');
meanPPC     = nan(size(uniqueFreq, 1), size(uniqueOsc, 1));

% loop through bars
for iBar = 1:size(uniqueFreq, 1)

    % loop through oscillation index
    for iOsc = 1:size(uniqueOsc, 1)

        % index
        bIdx                = strcmp(cellstr(thisData.ppcData.freqBand), char(uniqueFreq(iBar))) & ...
            strcmp(cellstr(thisData.ppcData.oscYesNo), char(uniqueOsc(iOsc)));

        % mean
        meanPPC(iBar, iOsc) = mean(thisData.ppcData.PPC(bIdx));
    end
end

% figure with first, second, and third axis
ppcFig          = figure;
ax1             = axes('units', 'normalized', 'Position', [0.2, 0.3, 0.7, 0.4]);
ax2             = axes('units', 'normalized', 'Position', [0.2, 0.7, 0.7, 0.2]);
ax3             = axes('units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.2]);

% hold on
hold(ax1, 'on');
hold(ax2, 'on');
hold(ax3, 'on');

% bar plots
b               = bar(ax1, [1, 2], meanPPC, 'FaceColor', 'flat');

% loop through conditions
nCond           = size(unique(thisData.ppcData.oscYesNo), 1);
nReg            = size(unique(thisData.ppcData.freqBand), 1);
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x-coordinate of the bars
    xBar(iCond, :) = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData = [0.5, 0.5, 0.5];
end

% plot distributions
for iBar = 1:size(uniqueFreq, 1)

    % loop through oscillation index
    for iOsc = 1:size(uniqueOsc, 1)

        % index
        bIdx            = strcmp(cellstr(thisData.ppcData.freqBand), char(uniqueFreq(iBar))) & ...
            strcmp(cellstr(thisData.ppcData.oscYesNo), char(uniqueOsc(iOsc)));

        % swarmchart
        s1              = swarmchart(ax1, repmat(xBar(iOsc, iBar), size(thisData.ppcData.PPC(bIdx))), thisData.ppcData.PPC(bIdx), ...
            'k', 'filled', ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
            'MarkerEdgeAlpha', 0.5);
        s1.SizeData          = 5;
        s1.XJitter           = 'rand';
        s1.XJitterWidth      = 0.15;

        % swarmchart
        s2              = swarmchart(ax2, repmat(xBar(iOsc, iBar), size(thisData.ppcData.PPC(bIdx))), thisData.ppcData.PPC(bIdx), ...
            'k', 'filled', ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
            'MarkerEdgeAlpha', 0.5);
        s2.SizeData          = 5;
        s2.XJitter           = 'rand';
        s2.XJitterWidth      = 0.15;

        % swarmchart
        s3              = swarmchart(ax3, repmat(xBar(iOsc, iBar), size(thisData.ppcData.PPC(bIdx))), thisData.ppcData.PPC(bIdx), ...
            'k', 'filled', ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
            'MarkerEdgeAlpha', 0.5);
        s3.SizeData          = 5;
        s3.XJitter           = 'rand';
        s3.XJitterWidth      = 0.15;
    end
end

% first axis settings
set(ax1, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [0, 0.1], ...
    'YTick', [0, 0.1]);
ylabel(ax1, 'PPC');
yline(ax1, 0.1, ':');

% second axis settings
tickValuesAx2   = 0.2:0.1:1;
set(ax2, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [0.1, 1], ...
    'YTick', tickValuesAx2, ...
    'YTickLabel', cat(2, repmat({''}, 1, numel(tickValuesAx2) - 1), {'1'}));

% third axis settings
tickValuesAx3   = -1:0.1:0;
set(ax3, 'TickDir', 'out', 'Box', 'off', ...
    'XColor', 'none', 'YLim', [-1, 0], ...
    'YTick', tickValuesAx3, ...
    'YTickLabel', cat(2, {'-1'}, repmat({''}, 1, numel(tickValuesAx3) - 1)))

% link axes
linkaxes([ax1, ax2, ax3], 'x');

% save figure
set(ppcFig, 'renderer', 'painters');
saveas(ppcFig, fullfile(paths.save, 'Fig_S6c.png'));

%% Figure S6d_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S6d_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S6d_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S6d_left_ppc.png'));

%% Figure S6d_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S6d_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S6d_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S6d_middle_ppc.png'));

%% Figure S6d_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S6d_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S6d_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S6d_right_ppc.png'));

%%=========================================================================
% Figure S7
%%=========================================================================

%% Figure S7a

% load data
thisData                = load(fullfile(paths.data, 'Fig_S7a.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.1]);
saveas(f, fullfile(paths.save, 'Fig_S7a.png'));

%% Figure S7b

% load data
thisData                = load(fullfile(paths.data, 'Fig_S7b.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.025]);
saveas(f, fullfile(paths.save, 'Fig_S7b.png'));

%% Figure S7c

% load data
thisData        = load(fullfile(paths.data, 'Fig_S7c.mat'));

% colors for baseline, encoding and recall
colorData       = [0.6, 0.4, 0.2; 0, 0.447, 0.741; 0.85, 0.325, 0.098];

% figure
percFig         = figure('units', 'centimeters', 'Position', [10, 10, 14, 10]);
b               = bar(thisData.unitsWithPhaseLocking.Row, thisData.unitsWithPhaseLocking.Variables, 'FaceColor', 'flat');
hold on;

% loop through conditions
[nCond, nReg]   = size(thisData.unitsWithPhaseLocking.Variables');
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x coordinate of the bars
    xBar(iCond, :) = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData = colorData(iCond, :);
end

% settings
set(gca, 'tickDir', 'out', 'box', 'off');
ylabel('Units with phase locking (%)'); % adjust y-axis
ylim([0, 100]);
yticks([0, 50, 100]);
saveas(percFig, fullfile(paths.save, 'Fig_S7c.png'));

%% Figure S7d_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7d_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.baselineHighPowerPPC, thisData.baselineLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7d_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineHighPowerPPC, thisData.baselineLowPowerPPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7d_left_ppc.png'));

%% Figure S7d_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7d_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.encodingHighPowerPPC, thisData.encodingLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7d_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingHighPowerPPC, thisData.encodingLowPowerPPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7d_middle_ppc.png'));

%% Figure S7d_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7d_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.recallHighPowerPPC, thisData.recallLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7d_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallHighPowerPPC, thisData.recallLowPowerPPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7d_right_ppc.png'));

%% Figure S7e_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7e_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.baselineHighSlopePPC, thisData.baselineLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7e_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineHighSlopePPC, thisData.baselineLowSlopePPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7e_left_ppc.png'));

%% Figure S7e_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7e_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.encodingHighSlopePPC, thisData.encodingLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7e_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingHighSlopePPC, thisData.encodingLowSlopePPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7e_middle_ppc.png'));

%% Figure S7e_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7e_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.recallHighSlopePPC, thisData.recallLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7e_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallHighSlopePPC, thisData.recallLowSlopePPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7e_right_ppc.png'));

%% Figure S7f_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7f_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7f_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7f_left_ppc.png'));

%% Figure S7f_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7f_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7f_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7f_middle_ppc.png'));

%% Figure S7f_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S7f_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7f_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7f_right_ppc.png'));

%% Figure S7g

% load data
thisData    = load(fullfile(paths.data, 'Fig_S7g.mat'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encSuccPPC, thisData.encFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7g_ppc.png'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[encRank, encT, surEncT]            = TG_PermTest1D_2PS_20230727(thisData.encSuccPPC, thisData.encFailPPC, params.nSur, randSeedNum);

% surrogate test results
encPval                             = (1 - encRank) * 2; % Bonferroni correction
if encPval > 1
    encPval = 1;
end

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surEncT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encT, 'color', 'k', 'LineWidth', 1);
xlim([-10, 10]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7g_sur.png'));

%% Figure S7h

% load data
thisData        = load(fullfile(paths.data, 'Fig_S7h.mat'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recSuccPPC, thisData.recFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S7h_ppc.png'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[recRank, recT, surRecT]            = TG_PermTest1D_2PS_20230727(thisData.recSuccPPC, thisData.recFailPPC, params.nSur, randSeedNum);

% surrogate test results
recPval                             = (1 - recRank) * 2; % Bonferroni correction
if recPval > 1
    recPval = 1;
end

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surRecT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(recT, 'color', 'k', 'LineWidth', 1);
xlim([-10, 10]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S7h_sur.png'));

%% Figure S7i

% load data
thisData                = load(fullfile(paths.data, 'Fig_S7i.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiff, params.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSig, params.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_S7i.png'));

%% Figure S7j

% load data
thisData                = load(fullfile(paths.data, 'Fig_S7j.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiffSucc, params.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSigSucc, params.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_S7j.png'));

%% Figure S7k

% load data
thisData                = load(fullfile(paths.data, 'Fig_S7k.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiffFail, params.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSigFail, params.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_S7k.png'));

%% Figure S7l

% load data
thisData        = load(fullfile(paths.data, 'Fig_S7l.mat'));

% figure
shiftFig        = figure('units', 'centimeters', 'position', [10, 10, 18, 12]);

% define categories for the x-axis
shiftBounds     = categorical({'All', 'PL', 'Succ', 'Succ PL', 'Fail', 'Fail PL'});

% reorder categories
desiredOrder    = {'All', 'PL', 'Succ', 'Succ PL', 'Fail', 'Fail PL'};
shiftBounds     = reordercats(shiftBounds, desiredOrder);

% define colors for each category
shiftColors     = [
    0.16, 0.16, 0.16; ...
    0.16, 0.16, 0.16; ...
    0.1, 0.6, 0.1; ...
    0.1, 0.6, 0.1; ...
    0.6, 0.1, 0.1; ...
    0.6, 0.1, 0.1
    ];

% define transparency levels for each category
shiftAlphas     = [0.2, 0.6, 0.2, 0.6, 0.2, 0.6];

% initialize bar objects
b               = gobjects(size(shiftBounds));

% clear current axes and begin plotting
hold on;

% loop through each category to create bars with corresponding properties
for iLine = 1:numel(shiftBounds)
    b(iLine) = bar(shiftBounds(iLine), thisData.shiftPercentage(iLine));
    set(b(iLine), 'FaceColor', shiftColors(iLine, :), 'FaceAlpha', shiftAlphas(iLine));
end

% adjustments
ylabel('Ratio');
ylim([0, 15]);
yticks(0:5:15);
set(gca, 'box', 'off', 'TickDir', 'out');

% save
saveas(shiftFig, fullfile(paths.save, 'Fig_S7l.png'));

%%=========================================================================
% Figure S8
%%=========================================================================

%% Figure S8a

% load data
thisData                = load(fullfile(paths.data, 'Fig_S8a.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.1]);
saveas(f, fullfile(paths.save, 'Fig_S8a.png'));

%% Figure S8b

% load data
thisData                = load(fullfile(paths.data, 'Fig_S8b.mat'));

% polar histogram
f                       = figure;
polarhistogram('BinEdges', thisData.angleHist.BinEdges, ...
    'BinCounts', thisData.angleHist.Values, ...
    'FaceColor', thisData.angleHist.FaceColor, 'FaceAlpha', thisData.angleHist.FaceAlpha, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', thisData.angleHist.polarHistLabels);
rlim([0, 0.025]);
saveas(f, fullfile(paths.save, 'Fig_S8b.png'));

%% Figure S8c

% load data
thisData        = load(fullfile(paths.data, 'Fig_S8c.mat'));

% colors for baseline, encoding and recall
colorData       = [0.6, 0.4, 0.2; 0, 0.447, 0.741; 0.85, 0.325, 0.098];

% figure
percFig         = figure('units', 'centimeters', 'Position', [10, 10, 14, 10]);
b               = bar(thisData.unitsWithPhaseLocking.Row, thisData.unitsWithPhaseLocking.Variables, 'FaceColor', 'flat');
hold on;

% loop through conditions
[nCond, nReg]   = size(thisData.unitsWithPhaseLocking.Variables');
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x coordinate of the bars
    xBar(iCond, :) = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData = colorData(iCond, :);
end

% settings
set(gca, 'tickDir', 'out', 'box', 'off');
ylabel('Units with phase locking (%)'); % adjust y-axis
ylim([0, 100]);
yticks([0, 50, 100]);
saveas(percFig, fullfile(paths.save, 'Fig_S8c.png'));

%% Figure S8d_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8d_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.baselineHighPowerPPC, thisData.baselineLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8d_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineHighPowerPPC, thisData.baselineLowPowerPPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8d_left_ppc.png'));

%% Figure S8d_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8d_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.encodingHighPowerPPC, thisData.encodingLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8d_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingHighPowerPPC, thisData.encodingLowPowerPPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8d_middle_ppc.png'));

%% Figure S8d_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8d_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.recallHighPowerPPC, thisData.recallLowPowerPPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
xlabel('t');
ylabel('Count');
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8d_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallHighPowerPPC, thisData.recallLowPowerPPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8d_right_ppc.png'));

%% Figure S8e_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8e_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.baselineHighSlopePPC, thisData.baselineLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8e_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineHighSlopePPC, thisData.baselineLowSlopePPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8e_left_ppc.png'));

%% Figure S8e_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8e_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.encodingHighSlopePPC, thisData.encodingLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8e_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingHighSlopePPC, thisData.encodingLowSlopePPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8e_middle_ppc.png'));

%% Figure S8e_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8e_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(thisData.recallHighSlopePPC, thisData.recallLowSlopePPC, params.nSur, randSeedNum);

% surrogate test results
topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(topLowT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8e_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallHighSlopePPC, thisData.recallLowSlopePPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8e_right_ppc.png'));

%% Figure S8f_left

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8f_left.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.6, 0.4, 0.2], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8f_left_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.baselineOscillationDetectedPPC, thisData.baselineNoOscillationDetectedPPC], [0.6, 0.4, 0.2; 0.6, 0.4, 0.2], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8f_left_ppc.png'));

%% Figure S8f_middle

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8f_middle.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0, 0.447, 0.741], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8f_middle_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encodingOscillationDetectedPPC, thisData.encodingNoOscillationDetectedPPC], [0, 0.447, 0.741; 0, 0.447, 0.741], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8f_middle_ppc.png'));

%% Figure S8f_right

% load data
thisData                            = load(fullfile(paths.data, 'Fig_S8f_right.mat'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[oscNoRank, oscNoT, surOscNoT]      = TG_PermTest1D_2PS_20230727(thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC, params.nSur, randSeedNum);

% surrogate test results
oscNoPval                           = (1 - oscNoRank) * 3; % Bonferroni correction

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surOscNoT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(oscNoT, 'color', [0.85, 0.325, 0.098], 'LineWidth', 1);
xlim([-15, 15]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8f_right_sur.png'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recallOscillationDetectedPPC, thisData.recallNoOscillationDetectedPPC], [0.85, 0.325, 0.098; 0.85, 0.325, 0.098], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8f_right_ppc.png'));

%% Figure S8g

% load data
thisData    = load(fullfile(paths.data, 'Fig_S8g.mat'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.encSuccPPC, thisData.encFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8g_ppc.png'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[encRank, encT, surEncT]            = TG_PermTest1D_2PS_20230727(thisData.encSuccPPC, thisData.encFailPPC, params.nSur, randSeedNum);

% surrogate test results
encPval                             = (1 - encRank) * 2; % Bonferroni correction
if encPval > 1
    encPval = 1;
end

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surEncT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encT, 'color', 'k', 'LineWidth', 1);
xlim([-10, 10]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8g_sur.png'));

%% Figure S8h

% load data
thisData        = load(fullfile(paths.data, 'Fig_S8h.mat'));

% plot PPC
ppcFig      = figure('units', 'centimeters', 'Position', [10, 10, 5, 8]);
TG_PlotPPC_20250515([thisData.recSuccPPC, thisData.recFailPPC], [0.1, 0.6, 0.1; 0.6, 0.1, 0.1], 'breakNegative');
saveas(ppcFig, fullfile(paths.save, 'Fig_S8h_ppc.png'));

% set random seed
randSeed                            = 444;
rng(randSeed, 'twister');
randSeedNum                         = randi(100000, 100000, 1);

% number of surrogates
params.nSur                         = 10001;

% permutation test
[recRank, recT, surRecT]            = TG_PermTest1D_2PS_20230727(thisData.recSuccPPC, thisData.recFailPPC, params.nSur, randSeedNum);

% surrogate test results
recPval                             = (1 - recRank) * 2; % Bonferroni correction
if recPval > 1
    recPval = 1;
end

% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
histogram(surRecT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(recT, 'color', 'k', 'LineWidth', 1);
xlim([-10, 10]);
ylim([0, 1000]);
set(gca,'TickDir', 'out');
box off;
saveas(surFigure, fullfile(paths.save, 'Fig_S8h_sur.png'));

%% Figure S8i

% load data
thisData                = load(fullfile(paths.data, 'Fig_S8i.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiff, params.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSig, params.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_S8i.png'));

%% Figure S8j

% load data
thisData                = load(fullfile(paths.data, 'Fig_S8j.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiffSucc, params.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSigSucc, params.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_S8j.png'));

%% Figure S8k

% load data
thisData                = load(fullfile(paths.data, 'Fig_S8k.mat'));

% parameters
params                  = struct();
params.polarHistEdges   = linspace(-pi, pi, 24);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

% figure
circDistSig             = figure;
circDiffHist            = polarhistogram(thisData.circDiffFail, params.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1], 'FaceAlpha', 0.2);
hold on;
circDiffSigHist         = polarhistogram(thisData.circDiffSigFail, params.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1]);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
saveas(circDistSig, fullfile(paths.save, 'Fig_S8k.png'));

%% Figure S8l

% load data
thisData        = load(fullfile(paths.data, 'Fig_S8l.mat'));

% figure
shiftFig        = figure('units', 'centimeters', 'position', [10, 10, 18, 12]);

% define categories for the x-axis
shiftBounds     = categorical({'All', 'PL', 'Succ', 'Succ PL', 'Fail', 'Fail PL'});

% reorder categories
desiredOrder    = {'All', 'PL', 'Succ', 'Succ PL', 'Fail', 'Fail PL'};
shiftBounds     = reordercats(shiftBounds, desiredOrder);

% define colors for each category
shiftColors     = [
    0.16, 0.16, 0.16; ...
    0.16, 0.16, 0.16; ...
    0.1, 0.6, 0.1; ...
    0.1, 0.6, 0.1; ...
    0.6, 0.1, 0.1; ...
    0.6, 0.1, 0.1
    ];

% define transparency levels for each category
shiftAlphas     = [0.2, 0.6, 0.2, 0.6, 0.2, 0.6];

% initialize bar objects
b               = gobjects(size(shiftBounds));

% clear current axes and begin plotting
hold on;

% loop through each category to create bars with corresponding properties
for iLine = 1:numel(shiftBounds)
    b(iLine) = bar(shiftBounds(iLine), thisData.shiftPercentage(iLine));
    set(b(iLine), 'FaceColor', shiftColors(iLine, :), 'FaceAlpha', shiftAlphas(iLine));
end

% adjustments
ylabel('Ratio');
ylim([0, 15]);
yticks(0:5:15);
set(gca, 'box', 'off', 'TickDir', 'out');

% save
saveas(shiftFig, fullfile(paths.save, 'Fig_S8l.png'));

%%=========================================================================
% Figure S9
%%=========================================================================

%% Figure S9a

% load data
thisData    = load(fullfile(paths.data, 'Fig_S9a.mat'));

% create figure
swFig       = figure('units', 'centimeters', 'position', [2, 2, 20, 6]);

for iGroup = 1:2
    if iGroup == 1
        histogram(thisData.allSpikeWidth(thisData.isInterneuron), ...
            0:0.05:1.5, 'FaceColor', thisData.swGroupColors{iGroup}, 'FaceAlpha', 0.5);
        hold on;
    else
        histogram(thisData.allSpikeWidth(~thisData.isInterneuron), ...
            0:0.05:1.5, 'FaceColor', thisData.swGroupColors{iGroup}, 'FaceAlpha', 0.5);
    end
end
xlabel('Spike width (ms)');
ylabel('Number of single units');
ylim([0, 100]);
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painter');

% save figure
saveas(swFig, fullfile(paths.save, 'Fig_S9a.png'));

%% Figure S9b

% load data
thisData    = load(fullfile(paths.data, 'Fig_S9b.mat'));

% plot
corrSwFig       = figure('units', 'centimeters', 'Position', [2, 2, 8, 6]);
plot(thisData.allSpikeWidth(thisData.isInterneuron, :), thisData.PPC(thisData.isInterneuron, :), 'x', 'Color', [thisData.swGroupColors{1}, 0.5]);
hold on;
plot(thisData.allSpikeWidth(~thisData.isInterneuron, :), thisData.PPC(~thisData.isInterneuron, :), 'x', 'Color', [thisData.swGroupColors{2}, 0.5]);
box off;
set(gca, 'TickDir', 'out', 'YScale', 'log');
xlabel('Spike width (ms)');
ylabel('PPC');
xlim([0, 2]);
ylim([1e-8, 1]);
title('Correlation spike width and phase locking');
set(gcf, 'Renderer', 'painters');
saveas(corrSwFig, fullfile(paths.save, 'Fig_S9b.png'));

%% Figure S9c

% load data
thisData        = load(fullfile(paths.data, 'Fig_S9c.mat'));

% plot percentage of significantly phase-locking units
percFig         = figure;
b1              = bar([thisData.percSigInt, thisData.percSigPyr], 'FaceColor', 'flat');
b1.CData(1, :)  = thisData.swGroupColors{1, 1};
b1.CData(2, :)  = thisData.swGroupColors{1, 2};
b1.FaceAlpha    = 0.5;
set(gca, 'XTickLabel', {'Narrow', 'Wide'});
set(gca,'TickDir', 'out');
box off;
ylim([0, 100]);
yticks([0, 50, 100]);
set(gcf, 'Renderer', 'painter');
saveas(percFig, fullfile(paths.save, 'Fig_S9c.png'));

%% Figure S9d

% load data
thisData            = load(fullfile(paths.data, 'Fig_S9d.mat'));

% parameters
params              = struct(); 
params.conditions   = {'HighWidth', 'LowWidth'};

% plot firing rates as a function of slope
for iCond = 1:2

    % params.conditions
    if strcmp(params.conditions{iCond}, 'LowWidth')
        thisColor       = thisData.swGroupColors{1};
        thisGroupIdx    = thisData.isInterneuron;
    elseif strcmp(params.conditions{iCond}, 'HighWidth')
        thisColor       = thisData.swGroupColors{2};
        thisGroupIdx    = ~thisData.isInterneuron;
    end

    % surrogate tests
    [slopeRank, slopeT, surSlopeT]  = TG_PermTest1D_2PS_20230727(thisData.topSlopeFr(thisGroupIdx), ...
        thisData.lowSlopeFr(thisGroupIdx), thisData.nSur, thisData.randSeedNum);

    % p-value
    slopeLowFrPval                  = min(slopeRank, 1 - slopeRank) * 2 * size(params.conditions, 2); % Bonferroni corrected

    % firing rates
    medianSlopeFr                   = median([thisData.topSlopeFr(thisGroupIdx), thisData.lowSlopeFr(thisGroupIdx)], 'omitnan');
    slopeFrFig                      = figure('units', 'centimeters', 'Position', [10, 10, 6, 12]);
    parallelcoords([thisData.topSlopeFr(thisGroupIdx), thisData.lowSlopeFr(thisGroupIdx)], 'Color', [thisColor, 0.5], ...
        'LineStyle', '-', 'LineWidth', 0.5);
    hold on;
    parallelcoords(medianSlopeFr, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3);
    set(gca, 'yscale', 'log');
    set(gca, 'XTickLabel',{'high', 'low'});
    xlim([0.75, 2.25]);
    ylim([0.1, 30]);
    set(gca, 'tickDir', 'out');
    set(gcf, 'renderer', 'painters');
    title(['P-value = ', num2str(slopeLowFrPval)]);
    saveas(slopeFrFig, fullfile(paths.save, strcat('Fig_S9d_', params.conditions{iCond}, '_medianFR.png')));

    % surrogate distributions
    surSlopeFigure                  = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    surSlopeTCounts                 = histcounts(surSlopeT, -10:0.1:10);
    surSlopeHist                    = histogram('BinCounts', surSlopeTCounts, 'BinEdges', -10:0.1:10, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(slopeT, 'Color', [thisColor, 0.5], 'LineWidth', 1.5);
    ylim([0, 500]);
    set(gca,'TickDir','out');
    box off;
    saveas(surSlopeFigure, fullfile(paths.save, strcat('Fig_S9d_', params.conditions{iCond}, '_surDistr.png')));
end

%%=========================================================================
% Figure S14
%%=========================================================================

%% Figure S14

% load data
thisData            = load(fullfile(paths.data, 'Fig_S14.mat'));

% figure
oscSlopeFig         = figure('units', 'centimeters', 'position', [2, 2, 15, 6]);
bD                  = bar(1:4, ...
    mean([thisData.oscillationDetectedHighSlope, ...
    thisData.oscillationDetectedLowSlope, ...
    thisData.noOscillationDetectedHighSlope, ...
    thisData.noOscillationDetectedLowSlope, ...
    ], 'omitnan'), ...
    'FaceColor', [0.5, 0.5, 0.5]);
hold on;
s1                  = swarmchart(ones(size(thisData.oscillationDetectedHighSlope)), thisData.oscillationDetectedHighSlope, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
    'MarkerEdgeAlpha', 0.5);
s1.SizeData          = 20;
s1.XJitter           = 'rand';
s1.XJitterWidth      = 0.5;
s2                  = swarmchart(ones(size(thisData.oscillationDetectedLowSlope)) * 2, thisData.oscillationDetectedLowSlope, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
    'MarkerEdgeAlpha', 0.5);
s2.SizeData          = 20;
s2.XJitter           = 'rand';
s2.XJitterWidth      = 0.5;
s3                  = swarmchart(ones(size(thisData.noOscillationDetectedHighSlope)) * 3, thisData.noOscillationDetectedHighSlope, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
    'MarkerEdgeAlpha', 0.5);
s3.SizeData          = 20;
s3.XJitter           = 'rand';
s3.XJitterWidth      = 0.5;
s4                  = swarmchart(ones(size(thisData.noOscillationDetectedLowSlope)) * 4, thisData.noOscillationDetectedLowSlope, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
    'MarkerEdgeAlpha', 0.5);
s4.SizeData          = 20;
s4.XJitter           = 'rand';
s4.XJitterWidth      = 0.5;
set(gca, 'tickDir', 'out', 'box', 'off');
saveas(oscSlopeFig, fullfile(paths.save, 'Fig_S14.png'));

%%=========================================================================
% Figure S15
%%=========================================================================

%% Figure S15a

% load data
thisData            = load(fullfile(paths.data, 'Fig_S15ab.mat'));

% slope
[encRecSlopeRank, encRecSlopeT, encRecSlopeSurT]     = TG_PermTest1D_2PS_20230727(thisData.thisSessEncSlope, thisData.thisSessRecSlope, thisData.nSur, thisData.randSeedNum);
[encBlnSlopeRank, encBlnSlopeT, encBlnSlopeSurT]     = TG_PermTest1D_2PS_20230727(thisData.thisSessEncSlope, thisData.thisSessBlnSlope, thisData.nSur, thisData.randSeedNum);
[recBlnSlopeRank, recBlnSlopeT, recBlnSlopeSurT]     = TG_PermTest1D_2PS_20230727(thisData.thisSessRecSlope, thisData.thisSessBlnSlope, thisData.nSur, thisData.randSeedNum);

% p-values
encRecSlopePval                     = min(encRecSlopeRank, 1 - encRecSlopeRank) * 2 * 3; % Bonferroni corrected
if encRecSlopePval > 1
    encRecSlopePval = 1;
end
encBlnSlopePval                     = min(encBlnSlopeRank, 1 - encBlnSlopeRank) * 2 * 3; % Bonferroni corrected
if encBlnSlopePval > 1
    encBlnSlopePval = 1;
end
recBlnSlopePval                     = min(recBlnSlopeRank, 1 - recBlnSlopeRank) * 2 * 3; % Bonferroni corrected
if recBlnSlopePval > 1
    recBlnSlopePval = 1;
end

% figure
slopeFig            = figure('units', 'centimeters', 'position', [2, 2, 8, 10]);
bD                  = bar(1:3, ...
    mean([thisData.thisSessEncSlope, ...
    thisData.thisSessRecSlope, ...
    thisData.thisSessBlnSlope], 'omitnan'), ...
    'FaceColor', 'flat');
bD.CData            = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.6, 0.4, 0.2];
hold on;
s1                  = swarmchart(ones(size(thisData.thisSessEncSlope)), thisData.thisSessEncSlope, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0, 0.447, 0.741], ...
    'MarkerEdgeAlpha', 0.5);
s1.SizeData          = 20;
s1.XJitter           = 'rand';
s1.XJitterWidth      = 0.5;
s2                  = swarmchart(ones(size(thisData.thisSessRecSlope)) * 2, thisData.thisSessRecSlope, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.85, 0.325, 0.098], ...
    'MarkerEdgeAlpha', 0.5);
s2.SizeData          = 20;
s2.XJitter           = 'rand';
s2.XJitterWidth      = 0.5;
s3                  = swarmchart(ones(size(thisData.thisSessBlnSlope)) * 3, thisData.thisSessBlnSlope, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.6, 0.4, 0.2], ...
    'MarkerEdgeAlpha', 0.5);
s3.SizeData          = 20;
s3.XJitter           = 'rand';
s3.XJitterWidth      = 0.5;
set(gca, 'XTickLabel', {'Encoding', 'Recall', 'Baseline'});
ylim([0, 2]);
yticks([0, 1, 2]);
ylabel('Aperiodic slope');
set(gca, 'tickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(slopeFig, fullfile(paths.save, 'Fig_S15a_left.png'));

% encoding versus recall slope
encRecSlopeFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encRecSlopeTCounts  = histcounts(encRecSlopeSurT, linspace(-4, 4, 201));
encRecSlopeHist     = histogram('BinCounts', encRecSlopeTCounts, 'BinEdges', linspace(-4, 4, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encRecSlopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encRecSlopeFigure, fullfile(paths.save, 'Fig_S15a_encRecSurDistr.png'));

% encoding versus baseline slope
encBlnSlopeFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encBlnSlopeTCounts  = histcounts(encBlnSlopeSurT, linspace(-4, 4, 201));
encBlnSlopeHist     = histogram('BinCounts', encBlnSlopeTCounts, 'BinEdges', linspace(-4, 4, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encBlnSlopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encBlnSlopeFigure, fullfile(paths.save, 'Fig_S15a_encBlnSurDistr.png'));

% recall versus baseline slope
recBlnSlopeFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
recBlnSlopeTCounts  = histcounts(recBlnSlopeSurT, linspace(-4, 4, 201));
recBlnSlopeHist     = histogram('BinCounts', recBlnSlopeTCounts, 'BinEdges', linspace(-4, 4, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(recBlnSlopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(recBlnSlopeFigure, fullfile(paths.save, 'Fig_S15a_recBlnSurDistr.png'));

%% Figure S15b

% load data
thisData            = load(fullfile(paths.data, 'Fig_S15ab.mat'));

% theta
[encRecThetaRank, encRecThetaT, encRecThetaSurT]     = TG_PermTest1D_2PS_20230727(thisData.thisSessEncTheta, thisData.thisSessRecTheta, thisData.nSur, thisData.randSeedNum);
[encBlnThetaRank, encBlnThetaT, encBlnThetaSurT]     = TG_PermTest1D_2PS_20230727(thisData.thisSessEncTheta, thisData.thisSessBlnTheta, thisData.nSur, thisData.randSeedNum);
[recBlnThetaRank, recBlnThetaT, recBlnThetaSurT]     = TG_PermTest1D_2PS_20230727(thisData.thisSessRecTheta, thisData.thisSessBlnTheta, thisData.nSur, thisData.randSeedNum);

% p-values
encRecThetaPval                     = min(encRecThetaRank, 1 - encRecThetaRank) * 2 * 3; % Bonferroni corrected
if encRecThetaPval > 1
    encRecThetaPval = 1;
end
encBlnThetaPval                     = min(encBlnThetaRank, 1 - encBlnThetaRank) * 2 * 3; % Bonferroni corrected
if encBlnThetaPval > 1
    encBlnThetaPval = 1;
end
recBlnThetaPval                     = min(recBlnThetaRank, 1 - recBlnThetaRank) * 2 * 3; % Bonferroni corrected
if recBlnThetaPval > 1
    recBlnThetaPval = 1;
end

% figure
thetaFig            = figure('units', 'centimeters', 'position', [2, 2, 8, 10]);
bD                  = bar(1:3, ...
    mean([thisData.thisSessEncTheta, ...
    thisData.thisSessRecTheta, ...
    thisData.thisSessBlnTheta], 'omitnan'), ...
    'FaceColor', 'flat');
bD.CData            = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.6, 0.4, 0.2];
hold on;
s1                  = swarmchart(ones(size(thisData.thisSessEncTheta)), thisData.thisSessEncTheta, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0, 0.447, 0.741], ...
    'MarkerEdgeAlpha', 0.5);
s1.SizeData          = 20;
s1.XJitter           = 'rand';
s1.XJitterWidth      = 0.5;
s2                  = swarmchart(ones(size(thisData.thisSessRecTheta)) * 2, thisData.thisSessRecTheta, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.85, 0.325, 0.098], ...
    'MarkerEdgeAlpha', 0.5);
s2.SizeData          = 20;
s2.XJitter           = 'rand';
s2.XJitterWidth      = 0.5;
s3                  = swarmchart(ones(size(thisData.thisSessBlnTheta)) * 3, thisData.thisSessBlnTheta, ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', [0.6, 0.4, 0.2], ...
    'MarkerEdgeAlpha', 0.5);
s3.SizeData          = 20;
s3.XJitter           = 'rand';
s3.XJitterWidth      = 0.5;
set(gca, 'XTickLabel', {'Encoding', 'Recall', 'Baseline'});
ylim([0, 80]);
ylabel('Oscillation detected (%)');
set(gca, 'tickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(thetaFig, fullfile(paths.save, 'Fig_S15b_left.png'));

% encoding versus recall theta
encRecThetaFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encRecThetaTCounts  = histcounts(encRecThetaSurT, linspace(-4, 4, 201));
encRecThetaHist     = histogram('BinCounts', encRecThetaTCounts, 'BinEdges', linspace(-4, 4, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encRecThetaT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encRecThetaFigure, fullfile(paths.save, 'Fig_S15b_encRecSurDistr.png'));

% encoding versus baseline theta
encBlnThetaFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encBlnThetaTCounts  = histcounts(encBlnThetaSurT, linspace(-4, 4, 201));
encBlnThetaHist     = histogram('BinCounts', encBlnThetaTCounts, 'BinEdges', linspace(-4, 4, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encBlnThetaT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encBlnThetaFigure, fullfile(paths.save, 'Fig_S15b_encBlnSurDistr.png'));

% recall versus baseline theta
recBlnThetaFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
recBlnThetaTCounts  = histcounts(recBlnThetaSurT, linspace(-4, 4, 201));
recBlnThetaHist     = histogram('BinCounts', recBlnThetaTCounts, 'BinEdges', linspace(-4, 4, 201), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(recBlnThetaT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(recBlnThetaFigure, fullfile(paths.save, 'Fig_S15b_recBlnSurDistr.png'));

%%=========================================================================
% Figure S16
%%=========================================================================

%% Figure S16ab

% load data
thisData            = load(fullfile(paths.data, 'Fig_S16ab.mat'));

% loop through brain regions
for iReg = 1:size(thisData.brainRegForLoop, 1)

    %  selection index for specific regions
    if strcmp(thisData.brainRegForLoop{iReg}, 'ALL')
        bSel = true(size(thisData.brainRegWireIdx));
    else
        bSel = strcmp(thisData.brainRegForLoop{iReg}, thisData.brainRegWireIdx);
    end

    % print histogram of theta frequencies for all peaks (normalized)
    binBurstFig = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    TG_ShadeSEM_20210714(thisData.binCenters, thisData.burstCounts{iReg, 1} * 100, 'k', 0.5);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Frequency (Hz)');
    ylabel('Probability (%)');
    xticks(0:2:10);
    ylim([0, 25]);
    grid on;
    title(thisData.brainRegForLoop{iReg});
    saveas(binBurstFig, fullfile(paths.save, strcat('Fig_S16a_', thisData.brainRegForLoop{iReg}, '_Oscillations.png')));

     % print histogram of theta frequencies for all cycles (normalized)
    binBurstFig = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    TG_ShadeSEM_20210714(thisData.binCenters, thisData.cycleCounts{iReg, 1} * 100, 'k', 0.5);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Frequency (Hz)');
    ylabel('Probability (%)');
    xticks(0:2:10);
    ylim([0, 25]);
    grid on;
    title(thisData.brainRegForLoop{iReg});
    saveas(binBurstFig, fullfile(paths.save, strcat('Fig_S16b_', thisData.brainRegForLoop{iReg}, '_Cycles.png')));
end

%%=========================================================================
% Figure S18
%%=========================================================================

%% Figure S18

% load data
thisData            = load(fullfile(paths.data, 'Fig_S18ab.mat'));

% loop through brain regions
for iReg = 1:size(thisData.brainRegion, 1)

    % get index for cells from this region
    regIdx                      = strcmp(thisData.allBrainRegIdxSplit(:, 1), thisData.brainRegion{iReg});
    regUnitNum                  = sum(regIdx);

    % encoding
    regEncSucc                  = thisData.encSucc(regIdx, :);
    regEncFail                  = thisData.encFail(regIdx, :);

    % retrieval
    regRecSucc                  = thisData.recSucc(regIdx, :);
    regRecFail                  = thisData.recFail(regIdx, :);

    % plot - encoding
    encRegFigSFC = figure;
    TG_ShadeSEM_20210714(thisData.frequencies, regEncSucc * 100, [0.1, 0.6, 0.1], 0.4);
    hold on;
    TG_ShadeSEM_20210714(thisData.frequencies, regEncFail * 100, [0.6, 0.1, 0.1], 0.4);
    xlabel('Frequency (Hz)');
    xlim([0, 100]);
    set(gca, 'XScale', 'log', 'tickDir', 'out');
    box off;
    ylabel('SFC (%)');
    ylim([0, 5]);
    title(['Encoding: ', thisData.brainRegion{iReg}]);

    % print figure
    saveas(encRegFigSFC, fullfile(paths.save, strcat('Fig_S18a_', thisData.brainRegion{iReg}, '.png')));

    % plot - retrieval
    recRegFigSFC = figure;
    TG_ShadeSEM_20210714(thisData.frequencies, regRecSucc * 100, [0.1, 0.6, 0.1], 0.4);
    hold on;
    TG_ShadeSEM_20210714(thisData.frequencies, regRecFail * 100, [0.6, 0.1, 0.1], 0.4);
    xlabel('Frequency (Hz)');
    xlim([0, 100]);
    set(gca, 'XScale', 'log', 'tickDir', 'out');
    box off;
    ylabel('SFC (%)');
    ylim([0, 5]);
    title(['Retrieval: ', thisData.brainRegion{iReg}]);

    % print figure
    saveas(recRegFigSFC, fullfile(paths.save, strcat('Fig_S18b_', thisData.brainRegion{iReg}, '.png')));
end

%%=========================================================================
% Figure S19
%%=========================================================================

%% Figure S19

% load data
thisData        = load(fullfile(paths.data, 'Fig_S19.mat'));

% loop through cells
for iCell = 1:size(thisData.spikeTimesInMilliseconds, 2)

    % unit name
    unitName        = ['unit', num2str(iCell)];

    % create rasterplot
    rasFig          = figure('units', 'centimeters', 'position', [10, 10, 12, 12]);
    subplot(2, 1, 2);
    hold on;

    % loop through segments
    for iSeg = 1:size(thisData.spikeTimesInMilliseconds, 1)

        % plot horizontal lines for each element
        yVal    = iSeg * ones(size(thisData.spikeTimesInMilliseconds{iSeg, iCell}));
        for iSpike = 1:size(thisData.spikeTimesInMilliseconds{iSeg, iCell}, 1)
            line([thisData.spikeTimesInMilliseconds{iSeg, iCell}(iSpike), ...
                thisData.spikeTimesInMilliseconds{iSeg, iCell}(iSpike)], ...
                [yVal(iSpike) - 0.4, yVal(iSpike) + 0.4], 'Color', 'k', 'LineWidth', 1);
        end
    end

    % adjust x-axis
    xlim([-1500, 3000]);
    xlabel('Time (ms)');
    xline([0, 1500], 'LineWidth', 1);
    ylim([0, size(thisData.spikeTimesInMilliseconds, 1) + 1]);
    set(gca, 'YDir', 'reverse', 'tickDir', 'out');
    yticks([1, 50, 100]);
    ylabel('Segment number');
    set(gca, 'FontSize', 14);

    %% add mean firing rate

    % initialize firing rate vector
    spikeHistAll = zeros(size(thisData.spikeTimesInMilliseconds, 1), size(thisData.timeVector, 2) - 1);

    % compute firing rate for each segment
    for iSeg = 1:size(thisData.spikeTimesInMilliseconds, 1)

        % spike times
        spikeTimes              = thisData.spikeTimesInMilliseconds{iSeg, iCell};

        % compute histogram for the current cell
        spikeHist               = histcounts(spikeTimes, thisData.timeVector);

        % subtract baseline firing rate
        blRate                  = mean(spikeHist(thisData.baselineIndex));
        spikeHist               = spikeHist - blRate;

        % accumulate spike histogram
        spikeHistAll(iSeg, :)   = spikeHist;
    end

    % compute mean firing rate and SEM
    meanFiringRate      = mean(spikeHistAll, 1);
    semFiringRate       = std(spikeHistAll, 1) / sqrt(size(spikeHistAll, 1));

    % time axis
    timeAxis            = (thisData.timeVector(1 : end - 1) + thisData.timeVector(2 : end)) / 2;

    % plot firing rate with shaded area for SEM
    subplot(2, 1, 1);
    hold on;
    TG_ShadeSEM_20210714(timeAxis, spikeHistAll, [0.4, 0.4, 0.4], 0.6);

    % firing rate
    xlim([-1500, 3000]);
    xline([0, 1500], 'LineWidth', 1);
    ylabel('Firing Rate (Hz)');
    ylim([-0.5, 1]);
    set(gca, 'tickDir', 'out');
    set(gca, 'FontSize', 14);

    % save
    saveas(rasFig, fullfile(paths.save, strcat('Fig_S19_', unitName, '_spiking.png')));

    %% spike density plot

    % spike density plot
    densityFig      = figure('units', 'centimeters', 'position', [10, 10, 4, 3]);
    spikeTime       = thisData.spikeTimeInSeconds .* 1000; % [msec]
    pcolor(spikeTime, thisData.ybins{iCell, 1}, thisData.density{iCell, 1}');
    shading interp;
    box off;
    colormap(densityFig, 'parula');
    saveas(densityFig, fullfile(paths.save, strcat('Fig_S19_', unitName, '_densityPlot.png')));
end

%%=========================================================================
% Figure S20
%%=========================================================================

%% Figure S20a

% load data
thisData                = load(fullfile(paths.data, 'Fig_S20a.mat'));

% create TFR figure - difference
TFR_Diff                = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgSuccVsFail);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Encoding - difference');
xline([0, 1.5], ':', 'LineWidth', 2);
xlim([-1, 2.5]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Diff, fullfile(paths.save, 'Fig_S20a_difference.png'));

% create TFR figure - successful
TFR_Succ     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgSucc);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Encoding - successful');
xline([0, 1.5], ':', 'LineWidth', 2);
xlim([-1, 2.5]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Succ, fullfile(paths.save, 'Fig_S20a_successful.png'));

% create TFR figure - fail
TFR_Fail     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgFail);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Encoding - unsuccessful');
xline([0, 1.5], ':', 'LineWidth', 2);
xlim([-1, 2.5]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Fail, fullfile(paths.save, 'Fig_S20a_unsuccessful.png'));

%% Figure S20b

% load data
thisData                = load(fullfile(paths.data, 'Fig_S20b.mat'));

% create TFR figure - difference
TFR_Diff                = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgSuccVsFail);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Object recall - difference');
xline([0, 4], ':', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Diff, fullfile(paths.save, 'Fig_S20b_difference.png'));

% create TFR figure - successful
TFR_Succ     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgSucc);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Object recall - successful');
xline([0, 4], ':', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Succ, fullfile(paths.save, 'Fig_S20b_successful.png'));

% create TFR figure - fail
TFR_Fail     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgFail);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Object recall - unsuccessful');
xline([0, 4], ':', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Fail, fullfile(paths.save, 'Fig_S20b_unsuccessful.png'));

%% Figure S20c

% load data
thisData                = load(fullfile(paths.data, 'Fig_S20c.mat'));

% create TFR figure - difference
TFR_Diff                = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgSuccVsFail);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Location recall - difference');
xline([0, 4], ':', 'LineWidth', 2);
xTicks                  = get(gca, 'XTick');
xticklabels             = get(gca, 'XTickLabel');
xticklabels_numeric     = str2double(xticklabels);
new_xticklabels_numeric = xticklabels_numeric - 4;
new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
set(gca, 'XTickLabel', new_xticklabels);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Diff, fullfile(paths.save, 'Fig_S20c_difference.png'));

% create TFR figure - successful
TFR_Succ     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgSucc);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Location recall - successful');
xline([0, 4], ':', 'LineWidth', 2);
xTicks                  = get(gca, 'XTick');
xticklabels             = get(gca, 'XTickLabel');
xticklabels_numeric     = str2double(xticklabels);
new_xticklabels_numeric = xticklabels_numeric - 4;
new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
set(gca, 'XTickLabel', new_xticklabels);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Succ, fullfile(paths.save, 'Fig_S20c_successful.png'));

% create TFR figure - fail
TFR_Fail     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
ft_singleplotTFR(thisData.Plotcfg, thisData.grAvgFail);
yticks(thisData.Plotcfg.yticks); % change positions of y ticks
set(gca, 'YTickLabel', thisData.Plotcfg.newYLabels); % change values of y ticks
title('Location recall - unsuccessful');
xline([0, 4], ':', 'LineWidth', 2);
xTicks                  = get(gca, 'XTick');
xticklabels             = get(gca, 'XTickLabel');
xticklabels_numeric     = str2double(xticklabels);
new_xticklabels_numeric = xticklabels_numeric - 4;
new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
set(gca, 'XTickLabel', new_xticklabels);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painters');
saveas(TFR_Fail, fullfile(paths.save, 'Fig_S20c_unsuccessful.png'));

%%=========================================================================
% Figure S21
%%=========================================================================

%% Figure S21a

% load data
thisData    = load(fullfile(paths.data, 'Fig_S21a.mat'));

% plot
f           = figure;
TG_PolarplotShifts_20241015(thisData.encAngleShift, thisData.recAngleShift, [0, 0.447, 0.741], [0.85, 0.325, 0.098], [0.16, 0.16, 0.16, 0.2], thisData.param);
shiftLines  = findall(gca, 'Type', 'line');
for iLine = 1:length(shiftLines)
    if thisData.isPhaseLocking(iLine)
        shiftLines(iLine).Color = [0.16, 0.16, 0.16, 1];
    end
end
saveas(f, fullfile(paths.save, 'Fig_S21a.png'));

%% Figure S21b

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21b.mat'));

% plot
f               = figure;
TG_PolarplotShifts_20241015(thisData.encSuccAngleShift, thisData.recSuccAngleShift, [0, 0.447, 0.741], [0.85, 0.325, 0.098], [0.1, 0.6, 0.1, 0.2], thisData.param);
shiftLinesSucc  = findall(gca, 'Type', 'line');
for iLine = 1:length(shiftLinesSucc)
    if thisData.isPhaseLockingSucc(iLine)
        shiftLinesSucc(iLine).Color = [0.1, 0.6, 0.1, 1];
    end
end
saveas(f, fullfile(paths.save, 'Fig_S21b.png'));

%% Figure S21c

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21c.mat'));

% plot
f               = figure;

TG_PolarplotShifts_20241015(thisData.encFailAngleShift, thisData.recFailAngleShift, [0, 0.447, 0.741], [0.85, 0.325, 0.098], [0.6, 0.1, 0.1, 0.2], thisData.param);
shiftLinesFail  = findall(gca, 'Type', 'line');
for iLine = 1:length(shiftLinesFail)
    if thisData.isPhaseLockingFail(iLine)
        shiftLinesFail(iLine).Color = [0.6, 0.1, 0.1, 1];
    end
end
saveas(f, fullfile(paths.save, 'Fig_S21c.png'));

%% Figure S21d_top

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21d_top.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21d_top.png'));

%% Figure S21d_bottom

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21d_bottom.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21d_bottom.png'));

%% Figure S21e_top

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21e_top.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21e_top.png'));

%% Figure S21e_bottom

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21e_bottom.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21e_bottom.png'));

%% Figure S21f_top

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21f_top.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21f_top.png'));

%% Figure S21f_bottom

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21f_bottom.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21f_bottom.png'));

%% Figure S21g_top

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21g_top.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21g_top.png'));

%% Figure S21g_bottom

% load data
thisData        = load(fullfile(paths.data, 'Fig_S21g_bottom.mat'));

% plot
f               = figure;
polarhistogram(thisData.circDiff, thisData.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
rlim([0, 10]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
title(['n = ', num2str(size(thisData.circDiff, 1))]);
saveas(f, fullfile(paths.save, 'Fig_S21g_bottom.png'));