%==========================================================================
% This script compares theta-phase locking for 1-10 Hz and 3-10 in the
% presence and absence of detected oscillations.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.fullPhaseRes      = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results'; % results folder 1-10 Hz
paths.highPhaseRes      = 'D:\TreasureHunt\PhaseAnalysis_20230921\3_10_Hz_PPC_Results'; % results folder 3-10 Hz
paths.save              = 'D:\TreasureHunt\PhaseAnalysisOscillations_20230718\1_10_Hz_3_10_Hz_PPC_Comparison'; % save folder
if ~isfolder(paths.save)
    mkdir(paths.save);
end

% parameters
params                  = [];
params.nSur             = 10001; % 10001

% set random seed
randSeed                = 444;
rng(randSeed, 'twister');
randSeedNum             = randi(100000, 100000, 1);

% own functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% load phase results data
additionalRes   = load(fullfile(paths.fullPhaseRes, 'additionalResultsSmall.mat'));
fullResults     = load(fullfile(paths.fullPhaseRes, 'allResults.mat'));
highResults     = load(fullfile(paths.highPhaseRes, 'allResults.mat'));

%% get data

% get data 1-10 Hz
fullComplex     = {fullResults.allRes.allComplex}';
fullOscIdx      = {fullResults.allRes.allBurstIdx}';
fullComplex     = fullComplex(~additionalRes.bExclude);
fullOscIdx      = fullOscIdx(~additionalRes.bExclude);

% get data 3-10 Hz
highComplex     = {highResults.allRes.allComplex}';
highOscIdx      = {highResults.allRes.allBurstIdx}';
highComplex     = highComplex(~additionalRes.bExclude);
highOscIdx      = highOscIdx(~additionalRes.bExclude);

% get phase 
fullPhases       = cellfun(@(x) angle(x), fullComplex, 'Uni', 0);
highPhases       = cellfun(@(x) angle(x), highComplex, 'Uni', 0);

%% calculate PPC and differences for oscillations present or not

% preallocate
fullOscPPC          = nan(size(fullPhases, 1), 1);
fullNonPPC          = nan(size(fullPhases, 1), 1);
highOscPPC          = nan(size(highPhases, 1), 1);
highNonPPC          = nan(size(highPhases, 1), 1);

% loop through cells
parfor iCell = 1:size(fullPhases, 1)

    % set random seed
    rng(randSeedNum(iCell, 1), 'twister');

    % display progress
    disp(iCell);

    %% get PPC for 1-10 Hz

    % select spikes by oscillation index
    fullOscPhases           = fullPhases{iCell, 1}(fullOscIdx{iCell, 1} == 1);
    fullNonPhases           = fullPhases{iCell, 1}(fullOscIdx{iCell, 1} == 0);
    
    % get PPC
    fullOscPPC(iCell, 1)    = TG_PPC_20241128(fullOscPhases);
    fullNonPPC(iCell, 1)    = TG_PPC_20241128(fullNonPhases);

    %% get PPC for 3-10 Hz

    % select spikes by oscillation index
    highOscPhases           = highPhases{iCell, 1}(highOscIdx{iCell, 1} == 1);
    highNonPhases           = highPhases{iCell, 1}(highOscIdx{iCell, 1} == 0);

    % get PPC
    highOscPPC(iCell, 1)    = TG_PPC_20241128(highOscPhases);
    highNonPPC(iCell, 1)    = TG_PPC_20241128(highNonPhases);
end

% prepare for 2x2 ANOVA
PPC         =  [fullOscPPC; fullNonPPC; highOscPPC; highNonPPC];

% predictor one: frequency band
freqBand    = categorical([repmat("1-10Hz", length(fullOscPPC) + length(fullNonPPC), 1); 
                 repmat("3-10Hz", length(highOscPPC) + length(highNonPPC), 1)]);

% predictor two: oscillation presence
oscYesNo    = categorical([repmat("Yes", length(fullOscPPC), 1);
               repmat("No", length(fullNonPPC), 1);
               repmat("Yes", length(highOscPPC), 1);
               repmat("No", length(highNonPPC), 1)]);

% session index
unitIdx             = cat(1, additionalRes.phaseRes.idx);
[~, ~, sessIdx]     = unique(unitIdx(:, 1:2), 'rows');
session             = categorical(repmat(sessIdx, 4, 1));

% fit linear mixed effects model
tbl         = table(freqBand, oscYesNo, session, PPC);
LMEPPC      = fitlme(tbl, 'PPC ~ 1 + freqBand * oscYesNo + (1 | session)');

% ANOVA to test significance
disp(anova(LMEPPC));

% t-test 1-10 Hz no oscillations vs. 3-10 Hz oscillations
[~, pOscNon, ~, statsOscNon]    = ttest(fullNonPPC, highOscPPC);

%% histogram to compare PPCs

% calculate mean and standard error of the mean
meanPPC     = mean([fullOscPPC, fullNonPPC, highOscPPC, highNonPPC]);
semPPC      = std([fullOscPPC, fullNonPPC, highOscPPC, highNonPPC]) / ...
    sqrt(length([fullOscPPC, fullNonPPC, highOscPPC, highNonPPC]));

% create figure
ppcFig      = figure;
ppcMean     = bar(meanPPC, 'FaceColor', [0.5, 0.5, 0.5]);
hold on;
ppcErr      = errorbar(meanPPC, semPPC, 'LineStyle', 'none', 'Color', 'k');
set(gca, 'XTickLabel', {'1-10 Hz Osc', '1-10 Hz No', '3-10 Hz Osc', '3-10 Hz No'});
set(gca,'TickDir', 'out');
box off;
ylim([0, 0.05]);
saveas(ppcFig, fullfile(paths.save, '1_10_Hz_3_10_Hz_meanPPC.svg'));
