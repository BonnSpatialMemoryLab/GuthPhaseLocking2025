%==========================================================================
% This script distinguishes between putative interneurons and putative
% pyramidal cells.
%
% References:
% - Gast et al., Clinical Neurophysiology, 2016
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes              = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results'; % results folder
paths.artifact              = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.phaseData             = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.SUData                = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.slopeRes              = 'D:\TreasureHunt\SPRiNTAnalysis_20231212'; % firing rate results for high and low slopes
paths.save                  = 'D:\TreasureHunt\CellTypeAnalysis_20241025'; % save folder

% parameters
params                      = struct();
params.spikeWidthThreshold  = 0.65; % in ms
params.nSur                 = 10001;
params.conditions           = {'HighWidth', 'LowWidth'};

% set random seed
randSeed    = 444;
rng(randSeed, 'twister');
randSeedNum = randi(100000, 100000, 1);

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

%% subjects
subjects    = {...
    'TH01'; ...
    'TH02'; ...
    'TH03'; ...
    'TH04'; ...
    'TH05'; ...
    'TH06'; ...
    'TH07'; ...
    'TH08'; ...
    'TH09'; ...
    'TH10'; ...
    'TH11'; ...
    'TH12'; ...
    'TH13'; ...
    'TH14'; ...
    'TH15'; ...
    'TH16'; ...
    'TH17'; ...
    'TH18'; ...
    };

%% load phase results data
results             = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));
% allResults          = load(fullfile(paths.phaseRes, 'allResults.mat'));

% get data from all results
allIdx              = cat(1, results.phaseRes.idx); % unit index

% exclude multi units
bSU                 = cat(1, results.phaseRes.bSingleUnit);
allIdx              = allIdx(bSU, :);

% preallocate
allRes              = cell(size(allIdx, 1), 1);

%% loop through units
parfor iCell = 1:size(allIdx, 1)
    
    % display progress
    disp(iCell);

    % set random seed
    rng(randSeedNum(iCell));
    
    %% get cell information

    % get index of this unit
    cellIdx                 = allIdx(iCell, :);
    
    % get subject and session name
    cellSub                 = subjects{cellIdx(1, 1)};
    cellSess                = strcat('session_', num2str(cellIdx(1, 2)));
    
    % get wire name
    mwChanDir               = TG_GetChanDir_20210812(paths.SUData, cellSub, cellSess);  % unit data
    cellWire                = mwChanDir(cellIdx(1, 3)).name;

    % get unit number
    cellUnit                = cellIdx(1, 4);
    
    %% load data

    % load sampling rate and sampleinfo
    phaseData               = load(fullfile(paths.phaseData, cellSub, cellSess, cellWire, 'datacutFt2000HzBP_generalized_1_10.mat'), 'fsample', 'sampleinfo');

    % load unit data
    t                       = load(fullfile(paths.SUData, cellSub, cellSess, cellWire, 'times_datacut.mat'));
    
    %% unit data

    % data for this cluster
    thisCluster             = t.cluster_class(t.cluster_class(:, 1) == cellUnit, 2) / 1000; % convert to seconds
    thisSpike               = t.spikes(t.cluster_class(:, 1) == cellUnit, :);

    %% firing rate

    % mean firing rate
    meanFR                  = size(thisCluster, 1) / (phaseData.sampleinfo(2) / phaseData.fsample);

    %% spike width

    % original wave form
    waveForm                = mean(thisSpike, 1);
    waveTime                = (0:(numel(waveForm) - 1)) / t.par.sr * 1000; % in milliseconds

    % global minimum
    [~, minIdx]             = min(waveForm); % global minimum
    tMin                    = waveTime(minIdx);

    % maximum after global minimum
    tmpWaveForm             = waveForm;
    tmpWaveForm(1:minIdx)  	= nan;
    [~, maxIdx]         	= max(tmpWaveForm);
    tMax                  	= waveTime(maxIdx);

    % spike width
    spikeWidth              = tMax - tMin;

    %% collect results

    % average waveform
    unitRes                 = [];
    unitRes.waveForm        = waveForm;
    unitRes.waveTime  	    = waveTime;
    unitRes.tMin            = tMin;
    unitRes.tMax            = tMax;

    % metrics for classification
    unitRes.meanFR   	    = meanFR; % Hz
    unitRes.spikeWidth	    = spikeWidth; % duration (ms)

    % collect results across units
    allRes{iCell, 1}        = unitRes;
end

%% concatenate results
allRes          = cat(1, allRes{:});

%% save results
save(fullfile(paths.save, 'results'), 'allRes', 'params');

%% load results

% results
r               = load(fullfile(paths.save, 'results.mat'));

% average firing rates and spike width
allMeanFR    	= cell2mat({r.allRes.meanFR}'); % non-normal distribution
allSpikeWidth   = cell2mat({r.allRes.spikeWidth}');

%% clustering based on spike width

% 1 = narrow spikes (putative interneurons); 2 = wider spikes (putative pyramidal cells)
swGroupIdx          = ones(size(allSpikeWidth));
swGroupIdx(allSpikeWidth > params.spikeWidthThreshold) = 2;

%% report

% report number of cells
fprintf('\nNumber of putative interneurons: %d. Number of putative pyramidal cells: %d.\n', ...
    sum(swGroupIdx == 1), sum(swGroupIdx == 2));

% report differences in firing rate
[~, p, ~, stats] = ttest2(allMeanFR(swGroupIdx == 1), allMeanFR(swGroupIdx == 2));
fprintf('Mean firing rate (Hz): %.3f for putative interneurons; %.3f for putative pyramidal cells; t(%d) = %.3f, p = %.3f.\n', ...
    mean(allMeanFR(swGroupIdx == 1)), mean(allMeanFR(swGroupIdx == 2)), ...
    stats.df, stats.tstat, p);

% report differences in spike widths
[~, p, ~, stats] = ttest2(allSpikeWidth(swGroupIdx == 1), allSpikeWidth(swGroupIdx == 2));
fprintf('Spike width (ms): %.3f for putative interneurons; %.3f for putative pyramidal cells; t(%d) = %.3f, p = %.3f.\n', ...
    mean(allSpikeWidth(swGroupIdx == 1)), mean(allSpikeWidth(swGroupIdx == 2)), ...
    stats.df, stats.tstat, p);

%% figure for spike width clustering

% group colors
swGroupColors       = {[0.450, 0.223, 0.745], [0.301, 0.756, 0.831]};

% create figure
swFig = figure('units', 'centimeters', 'position', [2, 2, 20, 6]);

for iGroup = min(swGroupIdx):max(swGroupIdx)
    histogram(allSpikeWidth(swGroupIdx == iGroup), 0:0.05:1.5, 'FaceColor', swGroupColors{iGroup}, 'FaceAlpha', 0.5);
    hold on;
end
xlabel('Spike width (ms)');
ylabel('Number of single units');
ylim([0, 100]);
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painter');

% explanation of spike width
axes('units', 'centimeters', 'position', [2.2, 3.5, 2, 1.3]);
hold on;
plot(r.allRes(1).waveTime, r.allRes(1).waveForm, 'Color', [0, 0, 0]);
tmpY = min(r.allRes(1).waveForm);
plot([r.allRes(1).tMin, r.allRes(1).tMin, r.allRes(1).tMin, r.allRes(1).tMax, r.allRes(1).tMax, r.allRes(1).tMax], ...
    [tmpY * 1.1, tmpY * 1.3, tmpY * 1.2, tmpY * 1.2, tmpY * 1.3, tmpY * 1.1], '-', 'Color', [0.5, 0.5, 0.5]);
axis off;
set(gcf, 'Renderer', 'painter');

% save figure
saveas(swFig, fullfile(paths.save, 'SpikeWidthClustering.svg'));

%% firing rate as a function of slope

% load slope results
slopeRes                                            = load(fullfile(paths.slopeRes, 'results.mat'));
topSlopeFr                                          = slopeRes.topSlopeFr(bSU);
lowSlopeFr                                          = slopeRes.lowSlopeFr(bSU);

%% plot firing rates as a function of slope

% loop through different conditions
for iCond = 1:size(params.conditions, 2)

    % params.conditions
    if strcmp(params.conditions{iCond}, 'LowWidth')
        thisColor       = swGroupColors{1};
        thisGroupIdx    = swGroupIdx == 1;
    elseif strcmp(params.conditions{iCond}, 'HighWidth')
        thisColor       = swGroupColors{2};
        thisGroupIdx    = swGroupIdx == 2;
    end

    % surrogate tests
    [slopeRank, slopeT, surSlopeT]  = TG_PermTest1D_2PS_20230727(topSlopeFr(thisGroupIdx), ...
        lowSlopeFr(thisGroupIdx), slopeRes.param.nSur, slopeRes.randSeedNum);

    % p-value
    slopeLowFrPval                  = min(slopeRank, 1 - slopeRank) * 2 * size(params.conditions, 2); % Bonferroni corrected

    % firing rates
    medianSlopeFr                   = median([topSlopeFr(thisGroupIdx), lowSlopeFr(thisGroupIdx)], 'omitnan');
    slopeFrFig                      = figure('units', 'centimeters', 'Position', [10, 10, 6, 12]);
    parallelcoords([topSlopeFr(thisGroupIdx), lowSlopeFr(thisGroupIdx)], 'Color', [thisColor, 0.5], ...
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
    saveas(slopeFrFig, fullfile(paths.save, strcat(params.conditions{iCond}, '_slopeMedianFiringRate.svg')));

    % surrogate distributions
    surSlopeFigure                  = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    surSlopeTCounts                 = histcounts(surSlopeT, slopeRes.param.surHistEdges);
    surSlopeHist                    = histogram('BinCounts', surSlopeTCounts, 'BinEdges', slopeRes.param.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(slopeT, 'Color', [thisColor, 0.5], 'LineWidth', 1.5);
    ylim([0, 500]);
    set(gca,'TickDir','out');
    box off;
    saveas(surSlopeFigure, fullfile(paths.save, strcat(params.conditions{iCond}, '_slopeSurDistrFiringRate.svg')));
end

%% theta-phase locking strength for putative interneurons versus pyramidal cells

% load all phase locking results
plRes               = load(fullfile(paths.phaseRes, 'allResults.mat'), 'allRes', 'allIdx');

% get all spike-associated angles
allComplex          = cat(1, {plRes.allRes.allComplex})';
allComplex          = allComplex(~results.bExclude);
allComplex          = allComplex(bSU);
allAngle            = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
allAngleCounts      = cellfun(@length, allAngle);
minAngleCount       = min(allAngleCounts);
allAnglePooled      = cat(1, allAngle{:});

% rank for significant phase locking
allRank             = cat(1, {plRes.allRes.allRank})';
allRank             = cell2mat(allRank(~results.bExclude));
allRank             = allRank(bSU);
bSig                = allRank > 0.95;

% percentage of significantly phase locking putative interneurons and pyramidal cells
percSigInt          = sum(bSig(swGroupIdx == 1)) / sum(swGroupIdx == 1);
percSigPyr          = sum(bSig(swGroupIdx == 2)) / sum(swGroupIdx == 2);

% test PPC of putative interneurons versus putative pyramidal cells
PPC                 = cell2mat(cellfun(@(x) TG_PPC_20241128(x), allAngle, 'Uni', 0));
[~, ~, ~, stats]    = ttest2(PPC(swGroupIdx == 1, :), PPC(swGroupIdx == 2, :));
tVal                = stats.tstat;

%% permutation test for different number of significant units between putative interneurons and pyramidal cells

% difference
percDiff            = percSigInt - percSigPyr;

% create surrogates
surPercDiff         = nan(params.nSur, 1);
for iSur = 1:params.nSur

    % randomly permute group index
    surSwGroupIdx           = swGroupIdx(randperm(size(swGroupIdx, 1)));

    % surrogate percentage of significantly phase locking putative interneurons and pyramidal cells
    surPercSigInt           = sum(bSig(surSwGroupIdx == 1)) / sum(surSwGroupIdx == 1);
    surPercSigPyr           = sum(bSig(surSwGroupIdx == 2)) / sum(surSwGroupIdx == 2);

    % difference in shifting units
    surPercDiff(iSur, 1)    = surPercSigInt -  surPercSigPyr;
end

% rank
permDiffRank        = sum(percDiff > surPercDiff) / params.nSur;

% p-value
permDiffPval        = min([permDiffRank, 1 - permDiffRank]) * 2;

%% surrogate test putative interneurons versus pyramidal cells
surTval     = nan(params.nSur, 1);
parfor iSur = 1:params.nSur

    % set random seed
    rng(randSeedNum(iSur), 'twister');

    % create surrogate index
    surSwGroupIdx                           = swGroupIdx(randperm(length(swGroupIdx)));

    % test putative interneurons versus pyramidal cells in surrogate data
    [~, ~, ~, surStats]                     = ttest2(PPC(surSwGroupIdx == 1, :), PPC(surSwGroupIdx == 2, :));
    surTval(iSur, 1)                        = surStats.tstat;
end

% percentage of pyramidal cells versus interneurons
intVsPyrRank            = sum(tVal > surTval) / params.nSur;
intVsPyrPval            = min([intVsPyrRank, 1 - intVsPyrRank]) * 2;

%% correlations

% Spearman's rho for firing rate and PPC
[frRho, frPval]         = corr(allMeanFR, PPC, 'type', 'Spearman');

% test Spearman's rho for spike width and PPC
[swRho, swPval]         = corr(allSpikeWidth, PPC, 'type', 'Spearman');

%% plot correlation between firing rate and phase locking
corrFrFig               = figure('units', 'centimeters', 'Position', [2, 2, 8, 6]);
plot(allMeanFR(swGroupIdx == 2, :), PPC(swGroupIdx == 2, :), 'x', 'Color', [swGroupColors{2}, 0.5]);
hold on;
plot(allMeanFR(swGroupIdx == 1, :), PPC(swGroupIdx == 1, :), 'x', 'Color', [swGroupColors{1}, 0.5]);
box off;
set(gca, 'TickDir', 'out', 'XScale', 'log', 'YScale', 'log');
xlabel('Mean firing rate (Hz)');
ylabel('PPC');
xlim([0.1, 100]);
ylim([1e-8, 1]);
title('Correlation firing rate and phase locking');
set(gcf, 'renderer', 'painters');
saveas(corrFrFig, fullfile(paths.save, 'CorrFiringRateAndPhaseLocking.svg'));

%% plot correlation between spike width and phase locking
corrSwFig       = figure('units', 'centimeters', 'Position', [2, 2, 8, 6]);
plot(allSpikeWidth(swGroupIdx == 1, :), PPC(swGroupIdx == 1, :), 'x', 'Color', [swGroupColors{1}, 0.5]);
hold on;
plot(allSpikeWidth(swGroupIdx == 2, :), PPC(swGroupIdx == 2, :), 'x', 'Color', [swGroupColors{2}, 0.5]);
box off;
set(gca, 'TickDir', 'out', 'YScale', 'log');
xlabel('Spike width (ms)');
ylabel('PPC');
xlim([0, 2]);
ylim([1e-8, 1]);
title('Correlation spike width and phase locking');
set(gcf, 'Renderer', 'painters');
saveas(corrSwFig, fullfile(paths.save, 'CorrSpikeWidthAndPhaseLocking.svg'));

%% plot percentage of significantly phase-locking units
percFig         = figure;
b1              = bar([percSigInt * 100, percSigPyr * 100], 'FaceColor', 'flat');
b1.CData(1, :)  = swGroupColors{1, 1};
b1.CData(2, :)  = swGroupColors{1, 2};
b1.FaceAlpha    = 0.5;
set(gca, 'XTickLabel', {'Narrow', 'Wide'});
set(gca,'TickDir', 'out');
box off;
ylim([0, 100]);
yticks([0, 50, 100]);
set(gcf, 'Renderer', 'painter');
saveas(percFig, fullfile(paths.save, 'IntVsPyrSigPerc.svg'));

%% plot PPC
ppcFig          = figure;
meanPPC         = [mean(PPC(swGroupIdx == 1)), mean(PPC(swGroupIdx == 2))];
semPPC          = [std(PPC(swGroupIdx == 1)) / sqrt(sum(swGroupIdx == 1)), std(PPC(swGroupIdx == 2)) / sqrt(sum(swGroupIdx == 2))];
b2              = bar(meanPPC, 'FaceColor', 'flat');
b2.CData(1, :)  = swGroupColors{1, 1};
b2.CData(2, :)  = swGroupColors{1, 2};
b2.FaceAlpha    = 0.5;
hold on;
ppcErr          = errorbar(meanPPC, semPPC, 'LineStyle', 'none', 'Color', 'k');
ylim([0, 0.03]);
set(gca, 'XTickLabel', {'Narrow', 'Wide'});
set(gca,'TickDir', 'out');
box off;
set(gcf, 'Renderer', 'painter');
saveas(ppcFig, fullfile(paths.save, 'IntVsPyrPPC.svg'));

%% plot histogram of surrogate and empirical t-value
surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
surHist            = histogram(surTval, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 0.5);
xlim([-10, 10]);
ylim([0, 1000]);
xline(mean(tVal), 'color', 'k', 'LineWidth', 1);
set(gca,'TickDir','out');
box off;
set(gcf, 'Renderer', 'painter');
saveas(surFigure, fullfile(paths.save, 'surDistrIntVsPyrPPC.svg'));
