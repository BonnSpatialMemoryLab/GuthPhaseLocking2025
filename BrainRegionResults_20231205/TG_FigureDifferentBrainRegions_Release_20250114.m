%=======================================================================
% This script creates a figure to show phase locking, firing rate,
% power, and artifact ratio in different brain regions.
%
% Tim Guth, 2025
%=======================================================================

%% settings
clc; close all; clear;

% set random seed
randSeedNum         = 444;
rng(randSeedNum, 'twister');
randSeedNum         = randi(100000, 100000, 1); % for randomseed

% paths
paths.load              = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results';
paths.save              = 'D:\TreasureHunt\BrainRegionResults_20231205\1_10_Hz_PPC_Results'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% own functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\BrainRegionResults_20231205'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% parameters
params                  = [];
params.numSubsamples    = 1000;
params.polarHistEdges   = linspace(-pi, pi, 21);
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};

%% load previously saved phase analysis results
results                 = load(fullfile(paths.load, 'additionalResultsSmall.mat'));
allResults              = load(fullfile(paths.load, 'allResults.mat'), 'allRes');
allResults              = allResults.allRes;
incResults              = allResults(~results.bExclude);
fprintf('Total number of cells: %d.\n', size(allResults, 1));

%% analysis

% subject index of cells
cellIdx                 = cat(1, incResults.idx);
sessIdx                 = cellIdx(:, 1:2);

% cells with significant phase locking
bBlSigCell              = cat(1, incResults.blRank) > 0.95;
bEncSigCell             = cat(1, results.phaseRes.encBothRank) > 0.95;
bRecSigCell             = cat(1, results.phaseRes.recBothRank) > 0.95;

% cell's brain regions
allBrainRegIdx          = {incResults.brainRegion}';
allBrainRegIdxSplit     = split(allBrainRegIdx, '_');
brainRegion             = unique(allBrainRegIdxSplit(:, 1));

% spike number
blSpikeNum              = cell2mat(cellfun(@(x) size(cat(1, x{:}), 1), {results.phaseRes.blSpikePower}', 'Uni', 0));
encSpikeNum             = cell2mat(cellfun(@(x) size(cat(1, x{:}), 1), {results.phaseRes.encSpikePower}', 'Uni', 0));
recSpikeNum             = cell2mat(cellfun(@(x) size(cat(1, x{:}), 1), {results.phaseRes.recSpikePower}', 'Uni', 0));

% pairwise phase consistency
blPPC                   = cat(1, incResults.blPPC);
encPPC                  = cat(1, incResults.encBothPPC);
recPPC                  = cat(1, incResults.recBothPPC);

% get power data
blPower                 = cellfun(@(x) cat(1, x{:}), {results.phaseRes.blSpikePower}', 'Uni', 0);
encPower                = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikePower}', 'Uni', 0);
recPower                = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikePower}', 'Uni', 0);

% mean log-power for each cell
meanLogPowerBl          = cell2mat(cellfun(@(x) mean(log10(x), 'omitnan'), blPower, 'Uni', 0));
meanLogPowerEnc         = cell2mat(cellfun(@(x) mean(log10(x), 'omitnan'), encPower, 'Uni', 0));
meanLogPowerRec         = cell2mat(cellfun(@(x) mean(log10(x), 'omitnan'), recPower, 'Uni', 0));

%% preallocate

% unit and session information
regUnitNum              = nan(size(brainRegion, 1), 1);
regSessNum              = nan(size(brainRegion, 1), 1);
regSessIdxCell          = cell(size(brainRegion, 1), 1);

% phase locking
regSigPLBl              = cell(size(brainRegion, 1), 1);
regSigPLEnc             = cell(size(brainRegion, 1), 1);
regSigPLRec             = cell(size(brainRegion, 1), 1);
regPL                   = nan(size(brainRegion, 1), 1);
regPLEnc                = nan(size(brainRegion, 1), 1);
regPLRec                = nan(size(brainRegion, 1), 1);

% power
regCellPowerBl          = cell(size(brainRegion, 1), 1);
regCellPowerEnc         = cell(size(brainRegion, 1), 1);
regCellPowerRec         = cell(size(brainRegion, 1), 1);
regPower                = nan(size(brainRegion, 1), 2);
regPowerEnc             = nan(size(brainRegion, 1), 2);
regPowerRec             = nan(size(brainRegion, 1), 2);

% spike number
regCellSpikeNumBl       = cell(size(brainRegion, 1), 1);
regCellSpikeNumEnc      = cell(size(brainRegion, 1), 1);
regCellSpikeNumRec      = cell(size(brainRegion, 1), 1);

% pairwise phase consistency
regCellPPCBl           = cell(size(brainRegion, 1), 1);
regCellPPCEnc          = cell(size(brainRegion, 1), 1);
regCellPPCRec          = cell(size(brainRegion, 1), 1);
regPPCBl               = nan(size(brainRegion, 1), 2);
regPPCEnc              = nan(size(brainRegion, 1), 2);
regPPCRec              = nan(size(brainRegion, 1), 2);

% loop through brain regions and count cells with significant phase locking
for iReg = 1:size(brainRegion, 1)

    % get index for cells from this region
    regIdx                      = strcmp(allBrainRegIdxSplit(:, 1), brainRegion{iReg});
    regUnitNum(iReg, 1)         = sum(regIdx);

    % get session information for those cells
    regSessIdx                  = sessIdx(regIdx, :);
    regSessIdxCell{iReg, 1}     = regSessIdx;
    regSessIdxUnique            = unique(regSessIdx, 'rows');

    % number of sessions with cells in this region
    regSessNum(iReg, 1)         = size(regSessIdxUnique, 1);
    
    %% phase locking
    
    % baseline
    regSigPLBl{iReg, 1}         = bBlSigCell(regIdx);
    regPL(iReg, 1)              = mean(regSigPLBl{iReg, 1}) * 100;

    % encoding
    regSigPLEnc{iReg, 1}        = bEncSigCell(regIdx);
    regPLEnc(iReg, 1)           = mean(regSigPLEnc{iReg, 1}) * 100;

    % recall
    regSigPLRec{iReg, 1}        = bRecSigCell(regIdx);
    regPLRec(iReg, 1)           = mean(regSigPLRec{iReg, 1}) * 100;

    %% power

    % baseline
    regCellPowerBl{iReg, 1}     = meanLogPowerBl(regIdx);
    regPower(iReg, 1)           = mean(regCellPowerBl{iReg, 1});
    regPower(iReg, 2)           = std(regCellPowerBl{iReg, 1}) / sqrt(length(regCellPowerBl{iReg, 1}(~isnan(regCellPowerBl{iReg, 1}))));

    % encoding
    regCellPowerEnc{iReg, 1}    = meanLogPowerEnc(regIdx);
    regPowerEnc(iReg, 1)        = mean(regCellPowerEnc{iReg, 1});
    regPowerEnc(iReg, 2)        = std(regCellPowerEnc{iReg, 1}) / sqrt(length(regCellPowerEnc{iReg, 1}(~isnan(regCellPowerEnc{iReg, 1}))));

    % recall
    regCellPowerRec{iReg, 1}    = meanLogPowerRec(regIdx);
    regPowerRec(iReg, 1)        = mean(regCellPowerRec{iReg, 1});
    regPowerRec(iReg, 2)        = std(regCellPowerRec{iReg, 1}) / sqrt(length(regCellPowerRec{iReg, 1}(~isnan(regCellPowerRec{iReg, 1}))));

    %% number of spikes

    % baseline
    regCellSpikeNumBl{iReg, 1}  = blSpikeNum(regIdx);

    % encoding
    regCellSpikeNumEnc{iReg, 1} = encSpikeNum(regIdx);

    % recall
    regCellSpikeNumRec{iReg, 1} = recSpikeNum(regIdx);

    %% Pairwise phase consistency

    % baseline
    regCellPPCBl{iReg, 1}       = blPPC(regIdx);
    regPPCBl(iReg, 1)           = mean(regCellPPCBl{iReg, 1});
    regPPCBl(iReg, 2)           = std(regCellPPCBl{iReg, 1}) / sqrt(length(regCellPPCBl{iReg, 1}(~isnan(regCellPPCBl{iReg, 1}))));

    % encoding
    regCellPPCEnc{iReg, 1}      = encPPC(regIdx);
    regPPCEnc(iReg, 1)          = mean(regCellPPCEnc{iReg, 1});
    regPPCEnc(iReg, 2)          = std(regCellPPCEnc{iReg, 1}) / sqrt(length(regCellPPCEnc{iReg, 1}(~isnan(regCellPPCEnc{iReg, 1}))));

    % recall
    regCellPPCRec{iReg, 1}      = recPPC(regIdx);
    regPPCRec(iReg, 1)          = mean(regCellPPCRec{iReg, 1});
    regPPCRec(iReg, 2)          = std(regCellPPCRec{iReg, 1}) / sqrt(length(regCellPPCRec{iReg, 1}(~isnan(regCellPPCRec{iReg, 1}))));
end

%% get results

% only include brain areas with data from at least 5 sessions
includeIdx                  = regSessNum >= 5;
incBrainRegCell             = brainRegion(includeIdx);
incBrainReg                 = categorical(incBrainRegCell);
incUnitNum                  = regUnitNum(includeIdx);
incSessIdxCell              = regSessIdxCell(includeIdx);

% colors for baseline, encoding and recall
colorData                   = [0.6, 0.4, 0.2; 0, 0.447, 0.741; 0.85, 0.325, 0.098];

%% create bar plot for percentage of phase-locking units
percFig                     = TG_CreateRegionFigure_20231205(incBrainReg, regPL, regPLEnc, regPLRec, includeIdx, colorData);
ylabel('Percentage of phase-locking units'); % adjust y-axis
ylim([0, 100]);
yticks([0, 50, 100]);
print(percFig, fullfile(paths.save, 'PhaseLockingInDifferentBrainRegions_20231017'), '-dsvg', '-r300');

%% create bar plot for average pairwise phase consistency
brainRegRayFig              = TG_CreateRegionFigure_20231205(incBrainReg, regPPCBl, regPPCEnc, regPPCRec, includeIdx, colorData);
ylabel('Pairwise phase consistency');
print(brainRegRayFig, fullfile(paths.save, 'PPCInDifferentBrainRegions_20231017'), '-dsvg', '-r300');

%% create box plot for spike numbers

% include only specific regions
incCellSpikeNum         = regCellSpikeNumBl(includeIdx);
incCellSpikeNumEnc      = regCellSpikeNumEnc(includeIdx);
incCellSpikeNumRec      = regCellSpikeNumRec(includeIdx);

% group data
numBoxplots             = reshape(1:sum(includeIdx) * 3, 3, sum(includeIdx));
allCatCellSpikeNum      = [];
allCatGroup             = [];
for iReg = 1:sum(includeIdx)

    % group data from baseline, encoding and recall
    catCellSpikeNum     = cat(1, incCellSpikeNum{iReg}, incCellSpikeNumEnc{iReg}, incCellSpikeNumRec{iReg});
    catGroup            = [ones(size(incCellSpikeNum{iReg})) * numBoxplots(1, iReg); ...
        ones(size(incCellSpikeNumEnc{iReg})) * numBoxplots(2, iReg); ...
        ones(size(incCellSpikeNumRec{iReg})) * numBoxplots(3, iReg)];

    % concatenate all data
    allCatCellSpikeNum  = cat(1, allCatCellSpikeNum, catCellSpikeNum);
    allCatGroup         = cat(1, allCatGroup, catGroup);
end

% boxplot figure
brainRegSpikeNumFig     = figure('unit', 'centimeters', 'position', [10, 10, 10, 10]);
bp                      = boxplot(allCatCellSpikeNum, allCatGroup);
h                       = findobj('LineStyle', '--');
set(h, 'LineStyle', '-');
h                       = findobj('Marker', '+');
set(h, 'Marker', '.', 'MarkerEdgeColor', 'k');

% adjust y-axis
set(gca, 'tickDir', 'out', 'yscale', 'log');
ylabel('Median spike number');
box off;

% print figure
print(brainRegSpikeNumFig, fullfile(paths.save, 'SpikeNumberInDifferentBrainRegions_20231017'), '-dsvg', '-r300');

%% create bar plot for average power
brainRegPowFig              = TG_CreateRegionFigure_20231205(incBrainReg, regPower, regPowerEnc, regPowerRec, includeIdx, colorData);

% adjust y-axis
ylabel('Average log-power');
ylim([0, 4]);

% print figure
print(brainRegPowFig, fullfile(paths.save, 'PowerInDifferentBrainRegions_20231017'), '-dsvg', '-r300');

%% Chi-Square test for percentages

% create contingency table - baseline
blTable                     = nan(2, size(incBrainReg, 1));
blTable(1, :)               = cell2mat(cellfun(@(x) sum(x), regSigPLBl(includeIdx), 'Uni', 0));
blTable(2, :)               = cell2mat(cellfun(@(x) size(x, 1), regSigPLBl(includeIdx), 'Uni', 0))' - blTable(1, :);

% perform chi-squared test
[pBl, qBl]                  = chi2test(blTable);

% create contingency table - encoding
encTable                    = nan(2, size(incBrainReg, 1));
encTable(1, :)              = cell2mat(cellfun(@(x) sum(x), regSigPLEnc(includeIdx), 'Uni', 0));
encTable(2, :)              = cell2mat(cellfun(@(x) size(x, 1), regSigPLEnc(includeIdx), 'Uni', 0))' - encTable(1, :);

% perform chi-squared test
[pEnc, qEnc]                = chi2test(encTable);

% create contingency table - recall
recTable                    = nan(2, size(incBrainReg, 1));
recTable(1, :)              = cell2mat(cellfun(@(x) sum(x), regSigPLRec(includeIdx), 'Uni', 0));
recTable(2, :)              = cell2mat(cellfun(@(x) size(x, 1), regSigPLRec(includeIdx), 'Uni', 0))' - recTable(1, :);

% perform chi-squared test
[pRec, qRec]                = chi2test(recTable);

% Bonferroni-correct p-values
pBl                         = pBl * 3;
pEnc                        = pEnc * 3;
pRec                        = pRec * 3;

%% LME for pairwise phase consistency and powers

% pairwise phase consistency data
ppcMatrix           = nan(sum(incUnitNum) * 3, 1);
ppcMatrix(:, 1)     = cat(1, ...
    cat(1, regCellPPCBl{includeIdx}), ...
    cat(1, regCellPPCEnc{includeIdx}), ...
    cat(1, regCellPPCRec{includeIdx}));

% organize power data
powerMatrix         = nan(sum(incUnitNum) * 3, 1);
powerMatrix(:, 1)   = cat(1, ...
    cat(1, regCellPowerBl{includeIdx}), ...
    cat(1, regCellPowerEnc{includeIdx}), ...
    cat(1, regCellPowerRec{includeIdx}));

% brain regions
regions             = categorical(repmat(repelem(incBrainReg, incUnitNum), 3, 1));

% experiment periods
periods             = categorical(repelem({'Baseline', 'Encoding', 'Recall'}, sum(incUnitNum))');

% sessions
incSessIdx          = cat(1, incSessIdxCell{:});
[~, ~, sessNum]     = unique(incSessIdx, 'rows', 'stable');
sessions            = categorical(repmat(sessNum, 3, 1));

% spike numbers
spikeNumbers        = cat(1, incCellSpikeNum{:}, incCellSpikeNumEnc{:}, incCellSpikeNumRec{:});

% data table for LME
data4LME            = table(ppcMatrix, regions, periods, sessions, spikeNumbers, powerMatrix, ...
    'VariableNames', {'ppc', 'regions', 'periods', 'sessions', 'spikenumber', 'spikepower'});

% model formula (fixed effect: regions, random effects: sessions and periods)
modelFormulaPPC             = 'ppc ~ 1 + regions + (1|sessions) + (1|periods)';
modelFormulaPower           = 'spikepower ~ 1 + regions + (1|sessions) + (1|periods)';

% ANOVA results pairwise phase consistency
LMEppc              = fitlme(data4LME, modelFormulaPPC);
statsLMEppc         = anova(LMEppc);
disp(statsLMEppc(2, :));

% ANOVA results power
LMEpower            = fitlme(data4LME, modelFormulaPower);
statsLMEpower       = anova(LMEpower);
disp(statsLMEpower(2, :));

% loop through different brain region references
refComparisons      = [];
refStatsPPC         = [];
refStatsPower       = [];
for iRef = 1:size(incBrainReg, 1)
    
    % reorder references
    regionOrder         = circshift(string(incBrainReg), iRef);
    data4LME.regions    = reordercats(data4LME.regions, regionOrder);

    % fit linear mixed effects model for ppc and power
    LMEppc              = fitlme(data4LME, modelFormulaPPC);
    LMEpower            = fitlme(data4LME, modelFormulaPower);

    % reference brain region
    referenceRegion     = regionOrder(1);

    % regions to compare to
    thisLMERegions      = cellfun(@(x) x(2), cellfun(@(x) strsplit(x, '_'), LMEppc.CoefficientNames(2:end), 'Uni', 0))';

    % collect region comparison names
    refComparisons      = cat(1, refComparisons, strcat(thisLMERegions(iRef:end), " vs. ", referenceRegion));

    % region comparison stats for ppc and power
    thisRefStatsPPC     = double(LMEppc.Coefficients(2:end, [4,5,6]));
    thisRefStatsPower   = double(LMEpower.Coefficients(2:end, [4,5,6]));

    % collect region comparison stats for ppc and power
    refStatsPPC         = cat(1, refStatsPPC, thisRefStatsPPC(iRef:end, :));
    refStatsPower       = cat(1, refStatsPower, thisRefStatsPower(iRef:end, :));
end

% sort data
[refComparisons, sortIdx]   = sort(refComparisons);
refStatsPPC                 = refStatsPPC(sortIdx, :);
refStatsPower               = refStatsPower(sortIdx, :);

% Bonferroni-correct and round P-values
refStatsPPC(:, 3)                           = round(refStatsPPC(:, 3) * size(refStatsPPC, 1), 3);
refStatsPower(:, 3)                         = round(refStatsPower(:, 3) * size(refStatsPower, 1), 3);
refStatsPPC(refStatsPPC(:, 3) > 1, 3)       = 1;
refStatsPower(refStatsPower(:, 3) > 1, 3)   = 1;

%% Kruskal-Wallis tests for spike numbers

% perform Kruskal-Wallis tests
[pBlSpikeNum, tblBlSpikeNum, statsBlSpikeNum]     = kruskalwallis(cat(1, incCellSpikeNum{:}), repelem(incBrainRegCell, incUnitNum));
[pEncSpikeNum, tblEncSpikeNum, statsEncSpikeNum]  = kruskalwallis(cat(1, incCellSpikeNumEnc{:}), repelem(incBrainRegCell, incUnitNum));
[pRecSpikeNum, tblRecSpikeNum, statsRecSpikeNum]  = kruskalwallis(cat(1, incCellSpikeNumRec{:}), repelem(incBrainRegCell, incUnitNum));

% Bonferroni-correct p-values
pBlSpikeNum                 = pBlSpikeNum * 3;
pEncSpikeNum                = pEncSpikeNum * 3;
pRecSpikeNum                = pRecSpikeNum * 3;
