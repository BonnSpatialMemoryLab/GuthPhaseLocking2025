%==========================================================================
% This script analyzes the relationship between power and
% neuronal phase locking
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisPower_20230718\1_10_Hz_PPC_Results'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% parameters
params                  = [];
params.polarHistEdges   = linspace(-pi, pi, 21);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
params.cond             = {'baseline'; 'enc'; 'rec'};
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
    
    %% get phases and power index for all spikes
    if strcmp(params.cond{iCond, 1}, 'baseline')
        
        % get data
        allComplex      = {allResults.allRes.allComplex}';
        allEncRecIdx    = {allResults.allRes.allEncRecIdx}';
        
        % exclude cells
        allComplex      = allComplex(~results.bExclude);
        allEncRecIdx    = allEncRecIdx(~results.bExclude);

        % exclude encoding and recall spikes
        allComplex      = cellfun(@(x, y) x(y == 0), allComplex, allEncRecIdx, 'Uni', 0);
        
        % get phase and power data
        allPhases       = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
        allPower        = cellfun(@(x) abs(x) .^ 2, allComplex, 'Uni', 0);

        % define color of histograms
        colorData       = [0.6, 0.4, 0.2];

    elseif strcmp(params.cond{iCond, 1}, 'enc')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encPhase}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikePower}', 'Uni', 0);
        colorData       = [0, 0.447, 0.741];
    elseif strcmp(params.cond{iCond, 1}, 'rec')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recPhase}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikePower}', 'Uni', 0);
        colorData       = [0.85, 0.325, 0.098];
    end

    %% calculate power index
    allPowerIdx     = cellfun(@(x) x > median(x), allPower, 'Uni', 0);

    %% calculate PPC and differences for oscillations present or not

    % preallocate
    topPhases       = cell(size(allPhases, 1), 1);
    lowPhases       = cell(size(allPhases, 1), 1);
    topAngle        = nan(size(allPhases, 1), 1);
    lowAngle        = nan(size(allPhases, 1), 1);
    topPPC          = nan(size(allPhases, 1), 1);
    lowPPC          = nan(size(allPhases, 1), 1);

    % loop through cells
    parfor iCell = 1:size(allPhases, 1)
        
        % set random seed
        rng(randSeedNum(iCell, 1), 'twister');

        % display progress
        disp(iCell);

        % select spikes by power index
        topPhases{iCell, 1}     = allPhases{iCell, 1}(allPowerIdx{iCell, 1} == 1);
        lowPhases{iCell, 1}     = allPhases{iCell, 1}(allPowerIdx{iCell, 1} == 0);
        
        % get mean and mean resultant vector length for each unit
        [~, topAngle(iCell, 1)] = circ_axialmean(topPhases{iCell, 1});
        [~, lowAngle(iCell, 1)] = circ_axialmean(lowPhases{iCell, 1});

        % get PPC
        topPPC(iCell, 1)        = TG_PPC_20241128(topPhases{iCell, 1});
        lowPPC(iCell, 1)        = TG_PPC_20241128(lowPhases{iCell, 1});
    end
    
    %% t-tests for phase locking across units

    % permutation test
    [topLowRank, topLowT, surTopLowT]   = TG_PermTest1D_2PS_20230727(topPPC, lowPPC, params.nSur, randSeedNum);

    % surrogate test results
    topLowPval                          = (1 - topLowRank) * 3; % Bonferroni correction

    %% Watson-William test

    % inclusion index
    incIdx              = ~isnan(topAngle) & ~isnan(lowAngle);

    % mean angular difference
    meanDiff            = circ_mean(angdiff(topAngle(incIdx), lowAngle(incIdx)));

    % empirical Watson-William test
    [~, wwTable]        = circ_wwtest(topAngle(incIdx), lowAngle(incIdx));
    topLowWw            = wwTable{2, 5};

    %% plot histogram of surrogate and empirical t-value
    surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    surHist            = histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(topLowT, 'color', colorData, 'LineWidth', 1);
    % limsx              = get(gca, 'XLim');
    % set(gca,'Xlim', [param.surHistEdges(1), limsx(2)]);
    xlim([-15, 15]);
    ylim([0, 1000]);
    set(gca,'TickDir','out');
    box off;
    saveas(surFigure, fullfile(paths.save, strcat(params.cond{iCond}, '_surDistr.svg')));

    %% histogram to compare PPC of different powers

    % calculate mean and standard error of the mean
    meanPPC     = mean([topPPC, lowPPC], 'omitnan');
    semPPC      = std([topPPC, lowPPC], 'omitnan') / sqrt(length([topPPC, lowPPC]));

    % create figure
    ppcFig      = figure;
    ppcMean     = bar(meanPPC, 'FaceColor', colorData);
    hold on;
    ppcErr      = errorbar(meanPPC, semPPC, 'LineStyle', 'none', 'Color', 'k');
    set(gca, 'XTickLabel', {'high', 'low'});
    ylim([0, 0.1]);
    set(gca,'TickDir','out');
    box off;
    title(strcat(params.cond{iCond}, 32, 'P =', 32, num2str(topLowPval)));
    saveas(ppcFig, fullfile(paths.save, strcat(params.cond{iCond}, '_meanPPC.svg')));
end
