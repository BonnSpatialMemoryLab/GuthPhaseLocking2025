%==========================================================================
% This script is for analyzing the relationship between the occurrence of
% theta oscillations and neuronal phase locking
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisOscillations_20230718\1_10_Hz_PPC_Results'; % save folder
if ~isfolder(paths.save)
    mkdir(paths.save);
end

% parameters
params                  = [];
params.polarHistEdges   = linspace(-pi, pi, 21);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
params.cond             = {'baseline'; 'enc'; 'rec'}; % or {'enc'; 'rec'}
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

    %% get phases and oscillation index for all spikes
    if strcmp(params.cond{iCond, 1}, 'baseline')

        % get data
        allComplex      = {allResults.allRes.allComplex}';
        allEncRecIdx    = {allResults.allRes.allEncRecIdx}';
        allOscIdx       = {allResults.allRes.allBurstIdx}';

        % exclude cells
        allComplex      = allComplex(~results.bExclude);
        allEncRecIdx    = allEncRecIdx(~results.bExclude);
        allOscIdx       = allOscIdx(~results.bExclude);

        % exclude encoding and recall spikes
        allComplex      = cellfun(@(x, y) x(y == 0), allComplex, allEncRecIdx, 'Uni', 0);
        allOscIdx       = cellfun(@(x, y) x(y == 0), allOscIdx, allEncRecIdx, 'Uni', 0);

        % get phase and power data
        allPhases       = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
        allPower        = cellfun(@(x) abs(x) .^ 2, allComplex, 'Uni', 0);

        % define color of histograms
        colorData       = [0.6, 0.4, 0.2];

    elseif strcmp(params.cond{iCond, 1}, 'enc')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encPhase}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encBurstIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikePower}', 'Uni', 0);
        colorData       = [0, 0.447, 0.741];
    elseif strcmp(params.cond{iCond, 1 }, 'rec')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recPhase}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recBurstIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikePower}', 'Uni', 0);
        colorData       = [0.85, 0.325, 0.098];
    end

    %% correlation between oscillation and power

    % Spearman correlation between power and theta index
    oscPowerCorr  = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allOscIdx, allPower, 'Uni', 0));
    oscPowerRho   = mean(oscPowerCorr, 'omitnan');

    %% calculate PPC and differences for oscillations present or not

    % preallocate
    oscPhases       = cell(size(allPhases, 1), 1);
    nonPhases       = cell(size(allPhases, 1), 1);
    oscAngle        = nan(size(allPhases, 1), 1);
    nonAngle        = nan(size(allPhases, 1), 1);
    oscPPC          = nan(size(allPhases, 1), 1);
    nonPPC          = nan(size(allPhases, 1), 1);

    % loop through cells
    parfor iCell = 1:size(allPhases, 1)
        
        % set random seed
        rng(randSeedNum(iCell, 1), 'twister');
        
        % display progress
        disp(iCell);

        % select spikes by oscillation index
        oscPhases{iCell, 1}         = allPhases{iCell, 1}(allOscIdx{iCell, 1} == 1);
        nonPhases{iCell, 1}         = allPhases{iCell, 1}(allOscIdx{iCell, 1} == 0);

        % get mean and mean resultant vector length for each unit
        [~, oscAngle(iCell, 1)]     = circ_axialmean(oscPhases{iCell, 1});
        [~, nonAngle(iCell, 1)]     = circ_axialmean(nonPhases{iCell, 1});

        % get PPC
        try
            oscPPC(iCell, 1)       = TG_PPC_20241128(oscPhases{iCell, 1});
            nonPPC(iCell, 1)       = TG_PPC_20241128(nonPhases{iCell, 1});
        catch
            oscPPC(iCell, 1)       = NaN;
            nonPPC(iCell, 1)       = NaN;
        end
    end

    %% inclusion index
    incIdx              = ~isnan(oscPPC) & ~isnan(nonPPC);
    oscPPC              = oscPPC(incIdx);
    nonPPC              = nonPPC(incIdx);
    
    %% t-tests for phase locking across units

    % permutation test
    [oscNonRank, oscNonT, surOscNonT]   = TG_PermTest1D_2PS_20230727(oscPPC, nonPPC, params.nSur, randSeedNum);

    % surrogate test results
    oscNonPval                          = (1 - oscNonRank) * 3; % Bonferroni correction

    %% Watson-William test

    % mean difference
    meanDiff            = circ_mean(angdiff(oscAngle(incIdx), nonAngle(incIdx)));

    % empirical Watson-William test
    [~, wwTable]        = circ_wwtest(oscAngle(incIdx), nonAngle(incIdx));
    oscNonWw            = wwTable{2, 5};
    
    %% plot histogram of surrogate and empirical t-value
    surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    surHist            = histogram(surOscNonT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xlim([-10, 10]);
    xticks([-10, 0, 10]);
    ylim([0, 1000]);
    xline(oscNonT, 'color', colorData, 'LineWidth', 1);
    set(gca,'TickDir','out');
    box off;
    saveas(surFigure, fullfile(paths.save, strcat(params.cond{iCond}, '_surDistrPPC.svg')));
    
    %% histogram to compare PPC of theta and no theta
    
    % calculate mean and standard error of the mean
    meanPPC     = mean([oscPPC, nonPPC]);
    semPPC      = std([oscPPC, nonPPC]) / sqrt(length([oscPPC, nonPPC]));

    % create figure
    ppcFig      = figure;
    ppcMean     = bar(meanPPC, 'FaceColor', colorData);
    hold on;
    ppcErr      = errorbar(meanPPC, semPPC, 'LineStyle', 'none', 'Color', 'k');
    set(gca, 'XTickLabel', {'theta', 'no theta'});
    set(gca,'TickDir', 'out');
    box off;
    ylim([0, 0.05]);
    title(strcat(params.cond{iCond}, 32, 'P =', 32, num2str(oscNonPval)));
    saveas(ppcFig, fullfile(paths.save, strcat(params.cond{iCond}, '_meanPPCBouts.svg')));
end