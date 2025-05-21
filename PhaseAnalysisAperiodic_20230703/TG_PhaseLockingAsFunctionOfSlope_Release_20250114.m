%==========================================================================
% This script is for analyzing the relationship between fooof slope
% and neuronal phase locking.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisAperiodic_20230703\1_10_Hz_PPC_Results'; % save folder
if ~isfolder(paths.save)
    mkdir(paths.save);
end

% parameters
params                  = [];
params.polarHistEdges   = linspace(-pi, pi, 21);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
params.cond             = {'baseline'; 'enc'; 'rec'}; % {'baseline'; 'enc'; 'rec'}; % or {'enc'; 'rec'};
params.nSur             = 10001; % 10001

% set random seed
randSeed    = 444;
rng(randSeed, 'twister');
randSeedNum = randi(100000, 100000, 1);

% own functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% load phase results data
results         = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));
if sum(strcmp(params.cond, 'baseline')) > 0
    allResults  = load(fullfile(paths.phaseRes, 'allResults.mat'));
end

%% loop through encoding and recall
for iCond = 1:size(params.cond, 1)

    %% get phases and slope index for all spikes
    if strcmp(params.cond{iCond, 1}, 'baseline')

        % get data
        allComplex      = {allResults.allRes.allComplex}';
        allEncRecIdx    = {allResults.allRes.allEncRecIdx}';
        allSlope        = {allResults.allRes.allSpikeSlope}';
        allOscIdx       = {allResults.allRes.allThetaIdx}';

        % exclude cells
        allComplex      = allComplex(~results.bExclude);
        allEncRecIdx    = allEncRecIdx(~results.bExclude);
        allSlope        = allSlope(~results.bExclude);
        allOscIdx       = allOscIdx(~results.bExclude);

        % exclude encoding and recall spikes
        allComplex      = cellfun(@(x, y) x(y == 0), allComplex, allEncRecIdx, 'Uni', 0);
        allSlope        = cellfun(@(x, y) x(y == 0), allSlope, allEncRecIdx, 'Uni', 0);
        allOscIdx       = cellfun(@(x, y) x(y == 0), allOscIdx, allEncRecIdx, 'Uni', 0);

        % get phase and power data
        allPhases       = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
        allPower        = cellfun(@(x) abs(x) .^ 2, allComplex, 'Uni', 0);

        % define color of histograms
        colorData       = [0.6, 0.4, 0.2];

    elseif strcmp(params.cond{iCond, 1}, 'enc')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encPhase}', 'Uni', 0);
        allSlope        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikeSlope}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encThetaIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikePower}', 'Uni', 0);
        colorData       = [0, 0.447, 0.741];
    elseif strcmp(params.cond{iCond, 1}, 'rec')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recPhase}', 'Uni', 0);
        allSlope        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikeSlope}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recThetaIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikePower}', 'Uni', 0);
        colorData       = [0.85, 0.325, 0.098];
    end

    %% calculate slope index and power index
    allSlopeIdx     = cellfun(@(x) double(x > median(x, 'omitnan')), allSlope, 'Uni', 0);
    allPowerIdx     = cellfun(@(x) double(x > median(x, 'omitnan')), allPower, 'Uni', 0);

    %% set slope index to NaN were slope results had NaNs
    for iCell = 1:size(allPhases, 1)
        
        % get index of NaNs
        slopeNanIdx                         = isnan(allSlope{iCell, 1});

        % adjust slope index
        allSlopeIdx{iCell, 1}(slopeNanIdx)  = NaN;
    end

    %% correlation

    % spearman correlation between slope and power
    slopePowCorr    = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allSlope, allPower, 'Uni', 0));
    slopePowRho     = mean(slopePowCorr, 'omitnan');

    % spearman correlation between slope and theta index
    slopeThetaCorr  = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allSlope, allOscIdx, 'Uni', 0));
    slopeThetaRho   = mean(slopeThetaCorr, 'omitnan');

    % spearman correlation between theta index and power
    oscPowerCorr    = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allOscIdx, allPower, 'Uni', 0));
    oscPowerRho     = mean(oscPowerCorr, 'omitnan');

    %% calculate PPCs and differences between different slopes

    % preallocate
    topPhases       = cell(size(allPhases, 1), 1);
    lowPhases       = cell(size(allPhases, 1), 1);
    topPPC          = nan(size(allPhases, 1), 1);
    lowPPC          = nan(size(allPhases, 1), 1);
    topAngle        = nan(size(allPhases, 1), 1);
    lowAngle        = nan(size(allPhases, 1), 1);

    % loop through cells
    parfor iCell = 1:size(allPhases, 1)

        % set random seed
        rng(randSeedNum(iCell, 1), 'twister');

        % display progress
        disp(iCell);

        % get slope index
        thisCellSlopeIdx         = allSlopeIdx{iCell, 1};

        % select spikes by slope index
        topPhases{iCell, 1}      = allPhases{iCell, 1}(thisCellSlopeIdx == 1);
        lowPhases{iCell, 1}      = allPhases{iCell, 1}(thisCellSlopeIdx == 0);

        % get mean and mean resultant vector length for each unit
        [~, topAngle(iCell, 1)]  = circ_axialmean(topPhases{iCell, 1});
        [~, lowAngle(iCell, 1)]  = circ_axialmean(lowPhases{iCell, 1});

        % get PPC
        try
            topPPC(iCell, 1)        = TG_PPC_20241128(topPhases{iCell, 1});
            lowPPC(iCell, 1)        = TG_PPC_20241128(lowPhases{iCell, 1});
        catch
            topPPC(iCell, 1)        = NaN;
            lowPPC(iCell, 1)        = NaN;
            continue;
        end
    end

    %% t-tests for phase locking across units
    
    % permutation test
    [topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(topPPC, lowPPC, params.nSur, randSeedNum);

    % surrogate test results
    topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

    %% Watson-William test

    % inclusion index
    incIdx              = ~isnan(topAngle) & ~isnan(lowAngle);

    % mean difference
    meanDiff            = circ_mean(angdiff(topAngle(incIdx), lowAngle(incIdx)));

    % empirical Watson-William test
    [~, wwTable]        = circ_wwtest(topAngle(incIdx), lowAngle(incIdx));
    topLowWw            = wwTable{2, 5};

    %% surrogate t-values of PPC
    
    % plot histogram of surrogate and empirical t-value
    surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    surHist            = histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(topLowT, 'color', colorData, 'LineWidth', 1);
    set(gca,'TickDir','out');
    box off;
    xlim([-10, 10]);
    ylim([0, 1000]);
    saveas(surFigure, fullfile(paths.save, strcat(params.cond{iCond}, '_surDistrPPC.svg')));

    %% histogram to compare PPCs of different slopes

    % calculate mean and standard error of the mean
    meanPPC     = mean([topPPC, lowPPC], 'omitnan');
    sem         = std([topPPC, lowPPC], 'omitnan') / sqrt(length([topPPC, lowPPC]));

    % create figure
    ppcFig      = figure;
    ppcMean     = bar(meanPPC, 'FaceColor', colorData);
    hold on;
    ppcErr      = errorbar(meanPPC, sem, 'LineStyle', 'none', 'Color', 'k');
    set(gca, 'XTickLabel', {'high', 'low'});
    set(gca,'TickDir','out');
    box off;
    ylim([0, 0.05]);
    title(strcat(params.cond{iCond}, 32, 'P =', 32, num2str(topLowPval)));
    saveas(ppcFig, fullfile(paths.save, strcat(params.cond{iCond}, '_meanPPC.svg')));
end
