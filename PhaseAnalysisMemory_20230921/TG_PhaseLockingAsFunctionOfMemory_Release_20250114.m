%==========================================================================
% This script applies surrogate tests for comparing phase locking between
% successful and unsuccessful encoding-recall pairs.
%
% Tim Guth, 2025
%==========================================================================

% start
clear; close all; clc;

% add paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results'; % phase results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisMemory_20230921\1_10_Hz_PPC_Results';
if ~isfolder(paths.save)
    mkdir(paths.save);
end

% add functions
addpath(genpath('D:\External\Functions\'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% load result
results                 = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

% set random seed
randSeed                = 444;
rng(randSeed, 'twister');
randSeedNum             = randi(100000, 100000, 1); % for randomseed

% settings
params.nSur             = 10001;
params.polarHistLabels  = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
params.locRecThreshold  = 'median'; % 'median', or 'absolute'
params.recallType       = 'all'; % 'all', 'object', or 'location'
params.phaseLocking     = 'all'; % 'all', 'sig', or 'no'
params.unitType         = 'all'; % 'all', 'single', or 'multi'
params.cellType         = 'all'; % 'all', 'object', or 'noobject'
params.powLvl           = 'all'; % 'all', 'high', or 'low'
params.slopeLvl         = 'all'; % 'all', 'high', or 'low'
params.thetaLvl         = 'all'; % 'all', 'theta' or 'notheta'
params.resultName       = strcat(params.recallType, 'Recall', '_', ...
    params.phaseLocking, 'PL', '_', ...
    params.unitType, 'Units', '_', ...
    params.cellType, 'Cells', '_', ...
    params.powLvl, 'Powers', '_', ...
    params.slopeLvl, 'Slopes', '_', ...
    params.thetaLvl, 'Osc', '_');

% include only specified units (all, single units or multi units)
if strcmp(params.unitType, 'all') % all cells
    bUnit = true(1, size([results.phaseRes], 1));
elseif strcmp(params.unitType, 'single')
    bUnit = [results.phaseRes.bSingleUnit];
elseif strcmp(params.unitType, 'multi')
    bUnit = ~[results.phaseRes.bSingleUnit];
end

% include only specified units (all, object cells, non-object cells)
if strcmp(params.cellType, 'all') % all cells
    bCell = true(1, size([results.phaseRes], 1));
elseif strcmp(params.cellType, 'object') % cells with increased firing rate during encoding
    bCell = [results.phaseRes.bObjectCell];
elseif strcmp(params.cellType, 'noobject')
    bCell = ~[results.phaseRes.bObjectCell];
end

% include only specified units (all, significantly phase locking, non-PL cells)
if strcmp(params.phaseLocking, 'all') % all cells
    bSig = true(1, size([results.phaseRes], 1));
elseif strcmp(params.phaseLocking, 'sig')
    bSig = [results.phaseRes.encBothRank] > 0.95 & [results.phaseRes.recBothRank] > 0.95;
elseif strcmp(params.phaseLocking, 'no')
    bSig = ~([results.phaseRes.encBothRank] > 0.95 & [results.phaseRes.recBothRank] > 0.95);
end

% extract phase results
bInc                    = bUnit & bCell & bSig;
phaseRes                = results.phaseRes(bInc);

% preallocate mean phases for succ. and unsucc. encoding and recall
encPPC                  = nan(size(phaseRes, 1), 1);
recPPC                  = nan(size(phaseRes, 1), 1);
encSuccPPC              = nan(size(phaseRes, 1), 1);
encFailPPC              = nan(size(phaseRes, 1), 1);
recSuccPPC              = nan(size(phaseRes, 1), 1);
recFailPPC              = nan(size(phaseRes, 1), 1);
encMrv                  = nan(size(phaseRes, 1), 1);
recMrv                  = nan(size(phaseRes, 1), 1);
encSuccMrv              = nan(size(phaseRes, 1), 1);
encFailMrv              = nan(size(phaseRes, 1), 1);
recSuccMrv              = nan(size(phaseRes, 1), 1);
recFailMrv              = nan(size(phaseRes, 1), 1);
encAngle                = nan(size(phaseRes, 1), 1);
recAngle                = nan(size(phaseRes, 1), 1);
encSuccAngle            = nan(size(phaseRes, 1), 1);
encFailAngle            = nan(size(phaseRes, 1), 1);
recSuccAngle            = nan(size(phaseRes, 1), 1);
recFailAngle            = nan(size(phaseRes, 1), 1);
allSizeEncSucc          = nan(size(phaseRes, 1), 1);
allSizeEncFail          = nan(size(phaseRes, 1), 1);
allSizeRecSucc          = nan(size(phaseRes, 1), 1);
allSizeRecFail          = nan(size(phaseRes, 1), 1);
surEncSuccAngle         = nan(size(phaseRes, 1), params.nSur);
surEncFailAngle         = nan(size(phaseRes, 1), params.nSur);
surRecSuccAngle         = nan(size(phaseRes, 1), params.nSur);
surRecFailAngle         = nan(size(phaseRes, 1), params.nSur);

% loop through units and unfold the trial-wise information
parfor iCell = 1:size(phaseRes, 1)

    % initialize variables
    bGoodMemEnc     = [];
    bGoodMemRec     = [];
    encRecallIdx    = [];
    recRecallIdx    = [];
    encPowerIdx     = [];
    recPowerIdx     = [];
    encSlopeIdx     = [];
    recSlopeIdx     = [];
    encThetaIdx     = [];
    recThetaIdx     = [];
    recNotEmptyIdx  = [];

    % set random seed
    rng(randSeedNum(iCell, 1), 'twister');

    % display progress
    disp(iCell);

    %% spike-wise memory performance

    % number of spikes
    encNumSpikes            = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).encSpikePower, 'UniformOutput', 0));
    recNumSpikes            = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).recSpikePower, 'UniformOutput', 0));

    % memory performance index per spike
    if strcmp(params.locRecThreshold, 'median')
        bGoodMemEnc         = repelem(phaseRes(iCell).bGoodMemEnc, encNumSpikes);
        bGoodMemRec         = repelem(phaseRes(iCell).bGoodMemRec, recNumSpikes);
    elseif strcmp(params.locRecThreshold, 'absolute')
        bGoodMemEnc         = repelem(phaseRes(iCell).bGoodMemEncAbs, encNumSpikes);
        bGoodMemRec         = repelem(phaseRes(iCell).bGoodMemRecAbs, recNumSpikes);
    end

    %% extract phases, power and slope

    % phases
    allEncPhase         = cat(1, phaseRes(iCell).encPhase{:});
    allRecPhase         = cat(1, phaseRes(iCell).recPhase{:});

    % power
    encPower            = cat(1, phaseRes(iCell).encSpikePower{:});
    recPower            = cat(1, phaseRes(iCell).recSpikePower{:});

    % theta
    encTheta            = cat(1, phaseRes(iCell).encBurstIdx{:});
    recTheta            = cat(1, phaseRes(iCell).recBurstIdx{:});

    % slope
    encSlope            = cat(1, phaseRes(iCell).encSpikeSlope{:});
    recSlope            = cat(1, phaseRes(iCell).recSpikeSlope{:});

    %% inclusion index for spikes according to settings

    % extract only spikes from specific recall type
    if strcmp(params.recallType, 'all')
        encRecallIdx    = ones(size(allEncPhase));
        recRecallIdx    = ones(size(allRecPhase));
    elseif strcmp(params.recallType, 'object')
        encRecallIdx    = repelem(phaseRes(iCell).encObjOrLoc, encNumSpikes) == 1;
        recRecallIdx    = repelem(phaseRes(iCell).recObjOrLoc, recNumSpikes) == 1;
    elseif strcmp(params.recallType, 'location')
        encRecallIdx    = repelem(phaseRes(iCell).encObjOrLoc, encNumSpikes) == 2;
        recRecallIdx    = repelem(phaseRes(iCell).recObjOrLoc, recNumSpikes) == 2;
    end

    % extract only spikes with specific power
    if strcmp(params.powLvl, 'all')
        encPowerIdx     = ones(size(encPower));
        recPowerIdx     = ones(size(recPower));
    elseif strcmp(params.powLvl, 'high')
        encPowerIdx     = encPower > median(encPower);
        recPowerIdx     = recPower > median(recPower);
    elseif strcmp(params.powLvl, 'low')
        encPowerIdx     = encPower <= median(encPower);
        recPowerIdx     = recPower <= median(recPower);
    end

    % extract only spikes with specific slope
    if strcmp(params.slopeLvl, 'all')
        encSlopeIdx     = ones(size(encSlope));
        recSlopeIdx     = ones(size(recSlope));
    elseif strcmp(params.slopeLvl, 'high')
        encSlopeIdx     = encSlope > median(encSlope, 'omitnan');
        recSlopeIdx     = recSlope > median(recSlope, 'omitnan');
    elseif strcmp(params.slopeLvl, 'low')
        encSlopeIdx     = encSlope <= median(encSlope, 'omitnan');
        recSlopeIdx     = recSlope <= median(recSlope, 'omitnan');
    end

    % extract only spikes with/without theta oscillations
    if strcmp(params.thetaLvl, 'all')
        encThetaIdx     = ones(size(allEncPhase));
        recThetaIdx     = ones(size(allRecPhase));
    elseif strcmp(params.thetaLvl, 'theta')
        encThetaIdx     = encTheta == 1;
        recThetaIdx     = recTheta == 1;
    elseif strcmp(params.thetaLvl, 'notheta')
        encThetaIdx     = encTheta == 0;
        recThetaIdx     = recTheta == 0;
    end

    % inclusion index
    encIncIdx           = encRecallIdx & encPowerIdx & encSlopeIdx & encThetaIdx;
    recIncIdx           = recRecallIdx & recPowerIdx & recSlopeIdx & recThetaIdx;

    %% separate successful and unsuccessful segments

    % successful and unsuccessful encoding
    encPhase            = allEncPhase(encIncIdx);
    encPhaseSucc        = allEncPhase(encIncIdx & bGoodMemEnc);
    encPhaseFail        = allEncPhase(encIncIdx & ~bGoodMemEnc);

    % successful and unsuccessful recall
    recPhase            = allRecPhase(recIncIdx);
    recPhaseSucc        = allRecPhase(recIncIdx & bGoodMemRec);
    recPhaseFail        = allRecPhase(recIncIdx & ~bGoodMemRec);

    % spike counts during successful and unsuccessful encoding
    allSizeEncSucc(iCell, 1)    = size(encPhaseSucc, 1);
    allSizeEncFail(iCell, 1)    = size(encPhaseFail, 1);

    % spike counts during successful and unsuccessful recall
    allSizeRecSucc(iCell, 1)    = size(recPhaseSucc, 1);
    allSizeRecFail(iCell, 1)    = size(recPhaseFail, 1);

    % concatenate successful and failed spikes
    encPhaseSuccFail    = cat(1, encPhaseSucc, encPhaseFail);
    recPhaseSuccFail    = cat(1, recPhaseSucc, recPhaseFail);

    % class variables
    encSuccFailClass    = cat(1 , ones(size(encPhaseSucc)), zeros(size(encPhaseFail)));
    recSuccFailClass    = cat(1 , ones(size(recPhaseSucc)), zeros(size(recPhaseFail)));

    try
        % get mean resultant vector length and mean angle
        [encMrv(iCell, 1), encAngle(iCell, 1)]         = circ_axialmean(encPhase);
        [recMrv(iCell, 1), recAngle(iCell, 1)]         = circ_axialmean(recPhase);
        [encSuccMrv(iCell, 1), encSuccAngle(iCell, 1)] = circ_axialmean(encPhaseSucc);
        [encFailMrv(iCell, 1), encFailAngle(iCell, 1)] = circ_axialmean(encPhaseFail);
        [recSuccMrv(iCell, 1), recSuccAngle(iCell, 1)] = circ_axialmean(recPhaseSucc);
        [recFailMrv(iCell, 1), recFailAngle(iCell, 1)] = circ_axialmean(recPhaseFail);

        % get PPC
        encPPC(iCell, 1)        = TG_PPC_20241128(encPhase);
        recPPC(iCell, 1)        = TG_PPC_20241128(recPhase);
        encSuccPPC(iCell, 1)    = TG_PPC_20241128(encPhaseSucc);
        encFailPPC(iCell, 1)    = TG_PPC_20241128(encPhaseFail);
        recSuccPPC(iCell, 1)    = TG_PPC_20241128(recPhaseSucc);
        recFailPPC(iCell, 1)    = TG_PPC_20241128(recPhaseFail);

        %% surrogate data

        % preallocate
        cellSurEncSuccAngle         = nan(1, params.nSur);
        cellSurEncFailAngle         = nan(1, params.nSur);
        cellSurRecSuccAngle         = nan(1, params.nSur);
        cellSurRecFailAngle         = nan(1, params.nSur);

        % loop through surrogates
        for iSur = 1:params.nSur

            % random shuffle of spike class assignment
            surEncSuccFailClass = datasample(encSuccFailClass, size(encSuccFailClass, 1), 'Replace', false);
            surRecSuccFailClass = datasample(recSuccFailClass, size(recSuccFailClass, 1), 'Replace', false);

            % create surrogates
            surEncPhaseSucc     = encPhaseSuccFail(surEncSuccFailClass == 1);
            surEncPhaseFail     = encPhaseSuccFail(surEncSuccFailClass == 0);
            surRecPhaseSucc     = recPhaseSuccFail(surRecSuccFailClass == 1);
            surRecPhaseFail     = recPhaseSuccFail(surRecSuccFailClass == 0);

            % get mean angle
            [~, cellSurEncSuccAngle(1, iSur)] = circ_axialmean(surEncPhaseSucc);
            [~, cellSurEncFailAngle(1, iSur)] = circ_axialmean(surEncPhaseFail);
            [~, cellSurRecSuccAngle(1, iSur)] = circ_axialmean(surRecPhaseSucc);
            [~, cellSurRecFailAngle(1, iSur)] = circ_axialmean(surRecPhaseFail);
        end

        %% collect surrogates across cells
        surEncSuccAngle(iCell, :)   = cellSurEncSuccAngle;
        surEncFailAngle(iCell, :)   = cellSurEncFailAngle;
        surRecSuccAngle(iCell, :)   = cellSurRecSuccAngle;
        surRecFailAngle(iCell, :)   = cellSurRecFailAngle;
    catch
        continue;
    end
end

%% exclude units with NaNs

% exclusion index
encExcludeIdx               = isnan(encSuccPPC) | isnan(encFailPPC);
recExcludeIdx               = isnan(recSuccPPC) | isnan(recFailPPC);

% exclude cells
encSuccPPC                  = encSuccPPC(~encExcludeIdx);
encFailPPC                  = encFailPPC(~encExcludeIdx);
recSuccPPC                  = recSuccPPC(~recExcludeIdx);
recFailPPC                  = recFailPPC(~recExcludeIdx);

%% t-tests for phase locking during encoding and recall across units

% permutation tests
[encRank, encT, surEncT]    = TG_PermTest1D_2PS_20230727(encSuccPPC, encFailPPC, params.nSur, randSeedNum);
[recRank, recT, surRecT]    = TG_PermTest1D_2PS_20230727(recSuccPPC, recFailPPC, params.nSur, randSeedNum);

% surrogate test results (Bonferroni correction)
encPval                     = (1 - encRank) * 2;
recPval                     = (1 - recRank) * 2;

%% Watson-William test

% mean difference
meanDiffEnc             = circ_mean(angdiff(encSuccAngle(~encExcludeIdx), encFailAngle(~encExcludeIdx)));
meanDiffRec             = circ_mean(angdiff(recSuccAngle(~recExcludeIdx), recFailAngle(~recExcludeIdx)));

% empirical Watson-William test
[~, wwEncTable]         = circ_wwtest(encSuccAngle(~encExcludeIdx), encFailAngle(~encExcludeIdx));
[~, wwRecTable]         = circ_wwtest(recSuccAngle(~recExcludeIdx), recFailAngle(~recExcludeIdx));
encWw                   = wwEncTable{2, 5};
recWw                   = wwRecTable{2, 5};

% surrogate Watson-William test
surEncWw = nan(params.nSur, 1);
surRecWw = nan(params.nSur, 1);
for iSur = 1:params.nSur

    % inclusion index
    surIncIdxEnc        = ~isnan(surEncSuccAngle(:, iSur)) & ~isnan(surEncFailAngle(:, iSur));
    surIncIdxRec        = ~isnan(surRecSuccAngle(:, iSur)) & ~isnan(surRecFailAngle(:, iSur));

    % surrogate Watson-William test
    [~, surEncWwTable]  = circ_wwtest(surEncSuccAngle(surIncIdxEnc, iSur), surEncFailAngle(surIncIdxEnc, iSur));
    [~, surRecWwTable]  = circ_wwtest(surRecSuccAngle(surIncIdxRec, iSur), surRecFailAngle(surIncIdxRec, iSur));
    surEncWw(iSur, 1)   = surEncWwTable{2, 5};
    surRecWw(iSur, 1)   = surRecWwTable{2, 5};
end

% surrogate test results - Watson-Williams tests
encWwRank       = sum(encWw > surEncWw) / params.nSur;
encWwPval       = (1 - encWwRank) * 2;
recWwRank       = sum(recWw > surRecWw) / params.nSur;
recWwPval       = (1 - recWwRank) * 2;

% additional Bonferroni correction, if two groups are compared
if ~strcmp(params.recallType, 'all') || ...
        ~strcmp(params.phaseLocking, 'all') || ...
        ~strcmp(params.unitType, 'all') || ...
        ~strcmp(params.cellType, 'all') || ...
        ~strcmp(params.powLvl, 'all') || ...
        ~strcmp(params.slopeLvl, 'all') || ...
        ~strcmp(params.thetaLvl, 'all')
    encPval     = encPval * 2;
    recPval     = recPval * 2;
    encWwPval   = encWwPval * 2;
    recWwPval   = recWwPval * 2;
end

%% plots for all spikes
if strcmp(params.recallType, 'all') && ...
        strcmp(params.phaseLocking, 'all') && ...
        strcmp(params.unitType, 'all') && ...
        strcmp(params.cellType, 'all') && ...
        strcmp(params.powLvl, 'all') && ...
        strcmp(params.slopeLvl, 'all') && ...
        strcmp(params.thetaLvl, 'all')

    %% plot mean angle and PPC distributions

    % plot all phases for successful encoding units
    encSuccPolar = figure;
    polarplot([encSuccAngle, encSuccAngle]', [zeros(size(encSuccAngle, 1), 1), encSuccMrv]', 'Color', [0.1, 0.6, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(encSuccAngle), circ_mean(encSuccAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preferred theta phase during successful encoding');
    print(encSuccPolar, fullfile(paths.save, 'allCells_encSucc_PolarHistogram'), '-dsvg', '-r300');

    % plot all phases for unsuccessful encoding units
    encFailPolar = figure;
    polarplot([encFailAngle, encFailAngle]', [zeros(size(encFailAngle, 1), 1), encFailMrv]', 'Color', [0.6, 0.1, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(encFailAngle), circ_mean(encFailAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preferred theta phase during unsuccessful encoding');
    print(encFailPolar, fullfile(paths.save, 'allCells_encFail_PolarHistogram'), '-dsvg', '-r300');

    % plot all phases for successful recall units
    recSuccPolar = figure;
    polarplot([recSuccAngle, recSuccAngle]', [zeros(size(recSuccAngle, 1), 1), recSuccMrv]', 'Color', [0.1, 0.6, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(recSuccAngle), circ_mean(recSuccAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preferred theta phase during successful recall');
    print(recSuccPolar, fullfile(paths.save, 'allCells_recSucc_PolarHistogram'), '-dsvg', '-r300');

    % plot all phases for unsuccessful recall units
    recFailPolar = figure;
    polarplot([recFailAngle, recFailAngle]', [zeros(size(recFailAngle, 1), 1), recFailMrv]', 'Color', [0.6, 0.1, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(recFailAngle), circ_mean(recFailAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', params.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preferred theta phase during unsuccessful recall');
    print(recFailPolar, fullfile(paths.save, 'allCells_recFail_PolarHistogram'), '-dsvg', '-r300');

    %% plot histogram of surrogate and empirical t-value - encoding
    encSurFigure       = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    encSurHist         = histogram(surEncT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(encT, 'color', 'k', 'LineWidth', 1);
    xlim([-10, 10]);
    ylim([0, 1000]);
    set(gca,'TickDir','out');
    box off;
    saveas(encSurFigure, fullfile(paths.save, strcat(params.resultName, '_EncodingSurDistr.svg')));

    %% plot histogram of surrogate and empirical t-value - recall
    recSurFigure       = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    recSurHist         = histogram(surRecT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(recT, 'color', 'k', 'LineWidth', 1);
    xlim([-10, 10]);
    ylim([0, 1000]);
    set(gca,'TickDir','out');
    box off;
    saveas(recSurFigure, fullfile(paths.save, strcat(params.resultName, '_RecallSurDistr.svg')));
end

%% plot mean PPC of successful and unsuccessful encoding

% figure
encPPCFig       = figure('units', 'centimeters', 'position', [2, 2, 10, 15]);

% mean
encSuccPPCMean  = mean(encSuccPPC);
encFailPPCMean  = mean(encFailPPC);

% standard error of the mean
semEncSucc      = std(encSuccPPC) / sqrt(size(encSuccPPC, 1));
semEncFail      = std(encFailPPC) / sqrt(size(encFailPPC, 1));

% bar plot
b1              = bar([encSuccPPCMean; encFailPPCMean], 'BarWidth', 0.6);
b1.FaceColor    = 'flat';
b1.CData        = [0.1, 0.6, 0.1; 0.6, 0.1, 0.1];
hold on;
e1              = errorbar([1; 2], [encSuccPPCMean; encFailPPCMean], [semEncSucc; semEncFail]);
e1.Color        = [0 0 0];
e1.LineStyle    = 'none';
if encPval < 0.05
    hold on;
    l1          = line([1, 2], [0.4, 0.4], 'LineWidth', 2, 'Color', 'k');
    a1          = plot(1.5, 0.41, '*', 'Color', 'k');
end
xlim([0.3, 2.7]);
if strcmp(params.recallType, 'all') && ...
        strcmp(params.phaseLocking, 'all') && ...
        strcmp(params.unitType, 'all') && ...
        strcmp(params.cellType, 'all') && ...
        strcmp(params.powLvl, 'all') && ...
        strcmp(params.slopeLvl, 'all') && ...
        strcmp(params.thetaLvl, 'all')
    ylim([0, 0.04]);
else
    ylim([0, 0.08]);
end
set(gca, 'xticklabel', {'Successful', 'Unsuccessful'}, 'TickDir', 'out', 'box', 'off');
ylabel('PPC');
title({'PPC (mean across units)',  'during encoding', ...
    strcat('(p = ', num2str(encPval), ')'), ...
    strcat('(n = ', num2str(sum(~encExcludeIdx)), ')')});

% print figure
print(encPPCFig, fullfile(paths.save, strcat(params.resultName, '_EncodingPPC')), '-dsvg', '-r300');

%% plot mean PPC of successful and unsuccessful recall

% figure
recPPCFig      = figure('units', 'centimeters', 'position', [2, 2, 10, 15]);

% mean
recSuccPPCMean  = mean(recSuccPPC);
recFailPPCMean  = mean(recFailPPC);

% standard error of the mean
semRecSucc      = std(recSuccPPC) / sqrt(size(recSuccPPC, 1));
semRecFail      = std(recFailPPC) / sqrt(size(recFailPPC, 1));

% bar plot
b1              = bar([recSuccPPCMean; recFailPPCMean], 'BarWidth', 0.6);
b1.FaceColor    = 'flat';
b1.CData        = [0.1, 0.6, 0.1; 0.6, 0.1, 0.1];
hold on;
e1              = errorbar([1; 2], [recSuccPPCMean; recFailPPCMean], [semRecSucc; semRecFail]);
e1.Color        = [0 0 0];
e1.LineStyle    = 'none';
if recPval < 0.05
    hold on;
    l1          = line([1, 2], [0.4, 0.4], 'LineWidth', 2, 'Color', 'k');
    a1          = plot(1.5, 0.41, '*', 'Color', 'k');
end
xlim([0.3, 2.7]);
if strcmp(params.recallType, 'all') && ...
        strcmp(params.phaseLocking, 'all') && ...
        strcmp(params.unitType, 'all') && ...
        strcmp(params.cellType, 'all') && ...
        strcmp(params.powLvl, 'all') && ...
        strcmp(params.slopeLvl, 'all') && ...
        strcmp(params.thetaLvl, 'all')
    ylim([0, 0.04]);
else
    ylim([0, 0.08]);
end
set(gca, 'xticklabel', {'Successful', 'Unsuccessful'}, 'TickDir', 'out', 'box', 'off');
ylabel('PPC');
title({'PPC (mean across units)',  'during recall', ...
    strcat('(p = ', num2str(recPval), ')'), ...
    strcat('(n = ', num2str(sum(~recExcludeIdx)), ')')});

% print figure
print(recPPCFig, fullfile(paths.save, strcat(params.resultName, '_RecallPPC')), '-dsvg', '-r300');