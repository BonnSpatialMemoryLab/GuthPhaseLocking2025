%=======================================================================
% This script tests for different phase-locking angles between
% encoding and retrieval.
%
% Tim Guth, 2025
%=======================================================================

% start
clear; close all; clc;

% add paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_PPC_Results'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisShifts_20231011\1_10_Hz_Results';
if ~isfolder(paths.save)
    mkdir(paths.save);
end

% add functions
addpath(genpath('D:\External\Functions\'));
addpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions');

% load result
additionalResults       = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

%% settings

% set random seed
randSeedNum             = 444;
rng(randSeedNum, 'twister');
randSeedNum             = randi(100000, 100000, 1); % for randomseed

% parameters
param                   = [];
param.nSur              = 10001;
param.polarHistEdges    = linspace(-pi, pi, 24);
param.polarHistLabels   = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.unitType          = 'all'; % 'all', 'single', or 'multi'
param.cellType          = 'all'; % 'all', 'object', or 'noobject'
param.powLvl            = 'all'; % 'all', 'high', or 'low'
param.slopeLvl          = 'all'; % 'all', 'high', or 'low'
param.freqLvl           = 'all'; % 'all', 'high', or 'low'
param.thetaLvl          = 'all'; % 'all', 'theta' or 'notheta'
param.plotVersion       = 'histogram'; % 'histogram', or 'line'
param.resultName        = strcat(param.unitType, 'Units', '_', ...
    param.cellType, 'Cells', '_', ...
    param.powLvl, 'Powers', '_', ...
    param.slopeLvl, 'Slopes', '_', ...
    param.thetaLvl, 'Osc', '_', ...
    param.freqLvl, 'Freq', '_', ...
    param.plotVersion, 'Plot');

% include only specified units (all, single units or multi units)
if strcmp(param.unitType, 'all') % all cells
    bUnit = true(1, size([additionalResults.phaseRes], 1));
elseif strcmp(param.unitType, 'single')
    bUnit = [additionalResults.phaseRes.bSingleUnit];
elseif strcmp(param.unitType, 'multi')
    bUnit = ~[additionalResults.phaseRes.bSingleUnit];
end

% include only specified cells
if strcmp(param.cellType, 'all') % all cells
    bSel = true(1, size([additionalResults.phaseRes], 1));
elseif strcmp(param.cellType, 'object') % cells with high firing rate during encoding
    bSel = [additionalResults.phaseRes.bObjectCell];
elseif strcmp(param.cellType, 'noobject')
    bSel = ~[additionalResults.phaseRes.bObjectCell];
end

% figure names
figureDir               = dir(fullfile(paths.phaseRes, 'UnitFigures'));
allFigureNames          = {figureDir.name}';

% extract phase results
bInc                    = bUnit & bSel;
phaseRes                = additionalResults.phaseRes(bInc);

% extract structure fields of interest
allIdx                  = cat(1, phaseRes.idx);

% preallocate mean phases for succ. and unsucc. encoding and recall
exclusionIndex          = false(size(phaseRes, 1), 1);
encAngle                = nan(size(phaseRes, 1), 1);
encSuccAngle            = nan(size(phaseRes, 1), 1);
encFailAngle            = nan(size(phaseRes, 1), 1);
recAngle                = nan(size(phaseRes, 1), 1);
recSuccAngle            = nan(size(phaseRes, 1), 1);
recFailAngle            = nan(size(phaseRes, 1), 1);
wwEncRec                = nan(size(phaseRes, 1), 1);
wwEncRecSucc            = nan(size(phaseRes, 1), 1);
wwEncRecFail            = nan(size(phaseRes, 1), 1);
surWwEncRec             = nan(size(phaseRes, 1), param.nSur);
surWwEncRecSucc         = nan(size(phaseRes, 1), param.nSur);
surWwEncRecFail         = nan(size(phaseRes, 1), param.nSur);
encRecRank              = nan(size(phaseRes, 1), 1);
encRecRankSucc          = nan(size(phaseRes, 1), 1);
encRecRankFail          = nan(size(phaseRes, 1), 1);

%% loop through units and unfold the trial-wise information
parfor iCell = 1:size(phaseRes, 1)

        % initialize variables
        encRecallIdx    = [];
        recRecallIdx    = [];
        encCatPowerIdx  = [];
        recCatPowerIdx  = [];
        encCatSlopeIdx  = [];
        recCatSlopeIdx  = [];
        encCatThetaIdx  = [];
        recCatThetaIdx  = [];
        encCatFreqIdx   = [];
        recCatFreqIdx   = [];
        smpEncPhaseSucc = [];
        smpEncPhaseFail = [];
        smpRecPhaseSucc = [];
        smpRecPhaseFail = [];

        % set random seed
        rng(randSeedNum(iCell, 1), 'twister');

        % display progress
        disp(iCell);

        %% number of spikes

        % number of spikes per trial
        encNumSpikes        = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).encSpikePower, 'UniformOutput', 0));
        recNumSpikes        = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).recSpikePower, 'UniformOutput', 0));

        %% extract phases, power and slope

        % phases
        allEncPhase         = phaseRes(iCell).encPhase;
        allRecPhase         = phaseRes(iCell).recPhase;

        % concatenate phases
        encCatPhase         = cat(1, allEncPhase{:});
        recCatPhase         = cat(1, allRecPhase{:});

        % concatenate power
        encCatPower         = cat(1, phaseRes(iCell).encSpikePower{:});
        recCatPower         = cat(1, phaseRes(iCell).recSpikePower{:});

        % theta
        encCatTheta         = cat(1, phaseRes(iCell).encBurstIdx{:});
        recCatTheta         = cat(1, phaseRes(iCell).recBurstIdx{:});

        % slope
        encCatSlope         = cat(1, phaseRes(iCell).encSpikeSlope{:});
        recCatSlope         = cat(1, phaseRes(iCell).recSpikeSlope{:});

        % instantaneous frequency
        encCatFreq         = cat(1, phaseRes(iCell).encSpikeFreq{:});
        recCatFreq         = cat(1, phaseRes(iCell).recSpikeFreq{:});

        %% inclusion index for spikes according to settings

        % extract only spikes with specific power
        if strcmp(param.powLvl, 'all')
            encCatPowerIdx     = ones(size(encCatPower));
            recCatPowerIdx     = ones(size(recCatPower));
        elseif strcmp(param.powLvl, 'high')
            encCatPowerIdx     = encCatPower > median(encCatPower);
            recCatPowerIdx     = recCatPower > median(recCatPower);
        elseif strcmp(param.powLvl, 'low')
            encCatPowerIdx     = encCatPower <= median(encCatPower);
            recCatPowerIdx     = recCatPower <= median(recCatPower);
        end

        % extract only spikes with specific slope
        if strcmp(param.slopeLvl, 'all')
            encCatSlopeIdx     = ones(size(encCatSlope));
            recCatSlopeIdx     = ones(size(recCatSlope));
        elseif strcmp(param.slopeLvl, 'high')
            encCatSlopeIdx     = encCatSlope > median(encCatSlope, 'omitnan');
            recCatSlopeIdx     = recCatSlope > median(recCatSlope, 'omitnan');
        elseif strcmp(param.slopeLvl, 'low')
            encCatSlopeIdx     = encCatSlope <= median(encCatSlope, 'omitnan');
            recCatSlopeIdx     = recCatSlope <= median(recCatSlope, 'omitnan');
        end

        % extract only spikes with specific instantaneous frequency
        if strcmp(param.freqLvl, 'all')
            encCatFreqIdx       = ones(size(encCatFreq));
            recCatFreqIdx       = ones(size(recCatFreq));
        elseif strcmp(param.freqLvl, 'high')
            encCatFreqIdx       = encCatFreq > median(encCatFreq, 'omitnan');
            recCatFreqIdx       = recCatFreq > median(recCatFreq, 'omitnan');
        elseif strcmp(param.freqLvl, 'low')
            encCatFreqIdx       = encCatFreq <= median(encCatFreq, 'omitnan');
            recCatFreqIdx       = recCatFreq <= median(recCatFreq, 'omitnan');
        end

        % extract only spikes with/without theta oscillations
        if strcmp(param.thetaLvl, 'all')
            encCatThetaIdx     = ones(size(encCatPhase));
            recCatThetaIdx     = ones(size(recCatPhase));
        elseif strcmp(param.thetaLvl, 'theta')
            encCatThetaIdx     = encCatTheta == 1;
            recCatThetaIdx     = recCatTheta == 1;
        elseif strcmp(param.thetaLvl, 'notheta')
            encCatThetaIdx     = encCatTheta == 0;
            recCatThetaIdx     = recCatTheta == 0;
        end

        %% include only selected results

        % concatenated inclusion index
        encCatIncIdx        = encCatPowerIdx & encCatSlopeIdx & encCatFreqIdx & encCatThetaIdx;
        recCatIncIdx        = recCatPowerIdx & recCatSlopeIdx & recCatFreqIdx & recCatThetaIdx;

        % recreate segment-wise results
        encIncIdx           = mat2cell(encCatIncIdx, encNumSpikes, 1);
        recIncIdx           = mat2cell(recCatIncIdx, recNumSpikes, 1);
        
        % include only specific phases
        encPhase            = cellfun(@(x, y) x(y == 1), allEncPhase, encIncIdx, 'Uni', 0);
        recPhase            = cellfun(@(x, y) x(y == 1), allRecPhase, recIncIdx, 'Uni', 0);

        %% separate successful and unsuccessful segments

        % memory performance index per spike
        bGoodMemEnc         = phaseRes(iCell).bGoodMemEnc;
        bGoodMemRec         = phaseRes(iCell).bGoodMemRec;

        % successful and unsuccessful encoding
        encPhaseSucc        = encPhase(bGoodMemEnc);
        encPhaseFail        = encPhase(~bGoodMemEnc);

        % successful and unsuccessful recall
        recPhaseSucc        = recPhase(bGoodMemRec);
        recPhaseFail        = recPhase(~bGoodMemRec);

        % get mean resultant vector length
        [encMrv, encAngle(iCell, 1)]            = circ_axialmean(cat(1, encPhase{:}));
        [encMrvSucc, encSuccAngle(iCell, 1)]    = circ_axialmean(cat(1, encPhaseSucc{:}));
        [encMrvFail, encFailAngle(iCell, 1)]    = circ_axialmean(cat(1, encPhaseFail{:}));
        [recMrv, recAngle(iCell, 1)]            = circ_axialmean(cat(1, recPhase{:}));
        [recMrvSucc, recSuccAngle(iCell, 1)]    = circ_axialmean(cat(1, recPhaseSucc{:}));
        [recMrvFail, recFailAngle(iCell, 1)]    = circ_axialmean(cat(1, recPhaseFail{:}));

        %% Watson-Williams test for phase shifts

        % concatenate phases of two groups
        encRecPhases                = cat(1, encPhase, recPhase);
        encRecPhasesSucc            = cat(1, encPhaseSucc, recPhaseSucc);
        encRecPhasesFail            = cat(1, encPhaseFail, recPhaseFail);

        % class variable
        classEncRec                 = [ones(size(encPhase, 1), 1); ones(size(recPhase, 1), 1) * 2];
        classEncRecSucc             = [ones(size(encPhaseSucc, 1), 1); ones(size(recPhaseSucc, 1), 1) * 2];
        classEncRecFail             = [ones(size(encPhaseFail, 1), 1); ones(size(recPhaseFail, 1), 1) * 2];

        % empirical Watson-Williams test
        try
            [~, wwEncRecTable]          = circ_wwtest(cat(1, encRecPhases{classEncRec == 1}), cat(1, encRecPhases{classEncRec == 2}));
            [~, wwEncRecTableSucc]      = circ_wwtest(cat(1, encRecPhasesSucc{classEncRecSucc == 1}), cat(1, encRecPhasesSucc{classEncRecSucc == 2}));
            [~, wwEncRecTableFail]      = circ_wwtest(cat(1, encRecPhasesFail{classEncRecFail == 1}), cat(1, encRecPhasesFail{classEncRecFail == 2}));
        catch
            exclusionIndex(iCell, 1) = true;
            continue;
        end

        % F-values
        wwEncRec(iCell, 1)          = wwEncRecTable{2, 5};
        wwEncRecSucc(iCell, 1)      = wwEncRecTableSucc{2, 5};
        wwEncRecFail(iCell, 1)      = wwEncRecTableFail{2, 5};

        % surrogate tests
        cellSurWwEncRec             = nan(param.nSur, 1);
        cellSurWwEncRecSucc         = nan(param.nSur, 1);
        cellSurWwEncRecFail         = nan(param.nSur, 1);
        for iSur = 1:param.nSur

            % shuffle class for surrogate dataset
            surClassEncRec                  = datasample(classEncRec, numel(classEncRec), 'replace', false);
            surClassEncRecSucc              = datasample(classEncRecSucc, numel(classEncRecSucc), 'replace', false);
            surClassEncRecFail              = datasample(classEncRecFail, numel(classEncRecFail), 'replace', false);

            % surrogate Watson-Williams test and surrogate F-values
            try
                [~, surWwEncRecTable]           = circ_wwtest(cat(1, encRecPhases{surClassEncRec == 1}), cat(1, encRecPhases{surClassEncRec == 2}));
                cellSurWwEncRec(iSur, 1)        = surWwEncRecTable{2, 5};
            catch
                cellSurWwEncRec(iSur, 1)        = NaN;
            end

            try
                [~, surWwEncRecTableSucc]       = circ_wwtest(cat(1, encRecPhasesSucc{surClassEncRecSucc == 1}), cat(1, encRecPhasesSucc{surClassEncRecSucc == 2}));
                cellSurWwEncRecSucc(iSur, 1)    = surWwEncRecTableSucc{2, 5};
            catch
                cellSurWwEncRecSucc(iSur, 1)    = NaN;
            end

            try
                [~, surWwEncRecTableFail]       = circ_wwtest(cat(1, encRecPhasesFail{surClassEncRecFail == 1}), cat(1, encRecPhasesFail{surClassEncRecFail == 2}));
                cellSurWwEncRecFail(iSur, 1)    = surWwEncRecTableFail{2, 5};
            catch
                cellSurWwEncRecFail(iSur, 1)    = NaN;
            end
        end

        % collect surrogates across cells
        surWwEncRec(iCell, :)       = cellSurWwEncRec;
        surWwEncRecSucc(iCell, :)   = cellSurWwEncRecSucc;
        surWwEncRecFail(iCell, :)   = cellSurWwEncRecFail;

        % rank of empirical WW-test F-value in surrogate dataset
        encRecRank(iCell, 1)        = sum(wwEncRec(iCell, 1) > cellSurWwEncRec) / param.nSur;
        encRecRankSucc(iCell, 1)    = sum(wwEncRecSucc(iCell, 1) > cellSurWwEncRecSucc) / param.nSur;
        encRecRankFail(iCell, 1)    = sum(wwEncRecFail(iCell, 1) > cellSurWwEncRecFail) / param.nSur;

        %% plot distribution of phase shifts

        % phase shift examples for figure
        examples = [5, 0, 32, 2; ...
            5, 0, 46, 2; ...
            10, 1, 34, 2; ...
            10, 1, 40, 4];

        % plot distribution
        if ismember(phaseRes(iCell).idx, examples, 'rows')  && ...
                strcmp(param.powLvl, 'all') && ...
                strcmp(param.slopeLvl, 'all') && ...
                strcmp(param.freqLvl, 'all') && ...
                strcmp(param.thetaLvl, 'all') && ...
                strcmp(param.unitType, 'all') && ...
                strcmp(param.cellType, 'all')

            % p-value
            wwP         = 1 - (sum(wwEncRec(iCell, 1) > surWwEncRec(iCell, :)) / param.nSur);

            % correct low p-values of permutation tests (see Phipson and Smyth, 2010)
            wwP(wwP < (1 / param.nSur)) = 1 / param.nSur;

            % plot
            permFig     = figure;
            permDistr   = histogram(surWwEncRec(iCell, :), 'FaceColor', [0.5, 0.5, 0.5]);
            hold on;
            empVal      = xline(wwEncRec(iCell, 1), 'k');
            title({strcat('{Permutation test (P = }', num2str(wwP, 3), ')')});
            xlabel('f-value Watson-Williams test');
            ylabel('Number of surrogates');
            box off;
            set(gca, 'tickDir', 'out');

            % print figure
            unitName    = strjoin(string(phaseRes(iCell).idx), '_');
            figureName  = strcat(unitName, '_SurrogateDistribution');
            set(permFig, 'renderer', 'painters');
            saveas(permFig, fullfile(paths.save, strcat(figureName, '.svg')), 'svg');

            % move figures of significant cells to separate folder
            if encRecRank(iCell, 1) > 0.95
                figIdx = strcmp(allFigureNames, strcat(regexprep(num2str(phaseRes(iCell).idx),'\s+','_'), '_both_PhaseAngleFigure.jpg'));
                copyfile(fullfile(figureDir(figIdx).folder, figureDir(figIdx).name), fullfile(paths.save, 'ShiftUnitFigures', figureDir(figIdx).name));
            elseif phaseRes(iCell).encBothRank > 0.95 && phaseRes(iCell).recBothRank > 0.95
                figIdx = strcmp(allFigureNames, strcat(regexprep(num2str(phaseRes(iCell).idx),'\s+','_'), '_both_PhaseAngleFigure.jpg'));
                copyfile(fullfile(figureDir(figIdx).folder, figureDir(figIdx).name), fullfile(paths.save, 'PhaseLockingUnitFigures', figureDir(figIdx).name));
            end
        end
end

%% find cells with significant phase shift
if strcmp(param.powLvl, 'all') && ...
        strcmp(param.slopeLvl, 'all') && ...
        strcmp(param.freqLvl, 'all') && ...
        strcmp(param.thetaLvl, 'all') && ...
        strcmp(param.unitType, 'all') && ...
        strcmp(param.cellType, 'all')

    %% identify cells with significant phase locking and phase shifts

    % significant phase locking - encoding
    bEncPL              = cat(1, phaseRes.encBothRank) > 0.95;
    bEncSuccPL          = cat(1, phaseRes.encSuccRank) > 0.95;
    bEncFailPL          = cat(1, phaseRes.encFailRank) > 0.95;

    % significant phase locking - recall
    bRecPL              = cat(1, phaseRes.recBothRank) > 0.95;
    bRecSuccPL          = cat(1, phaseRes.recSuccRank) > 0.95;
    bRecFailPL          = cat(1, phaseRes.recFailRank) > 0.95;

    % phase locking during encoding and recall
    bBothPL             = bEncPL & bRecPL;
    bSuccPL             = bEncSuccPL & bRecSuccPL;
    bFailPL             = bEncFailPL & bRecFailPL;

    % phase difference
    bShift              = encRecRank > 0.95;
    bShiftSucc          = encRecRankSucc > 0.95;
    bShiftFail          = encRecRankFail > 0.95;

    % phase locking and phase difference
    bPlShift            = bBothPL & bShift;
    bPlShiftSucc        = bSuccPL & bShiftSucc;
    bPlShiftFail        = bFailPL & bShiftFail;

    % index only for phase locking neurons
    bOnlyPlShift        = bShift(bBothPL);
    bOnlyPlShiftSucc    = bShiftSucc(bSuccPL);
    bOnlyPlShiftFail    = bShiftFail(bFailPL);

    % correlation between WW-test ranks of succ and unsucc memory performance
    [succFailShiftRho, succFailShiftPval]       = corr(encRecRankSucc, encRecRankFail);

    % binomial tests
    allPval             = myBinomTest(sum(bShift), size(bShift, 1), 0.05, 'one');
    shiftPval           = myBinomTest(sum(bPlShift), sum(bBothPL), 0.05, 'one');
    succPval            = myBinomTest(sum(bShiftSucc), size(encRecRankSucc, 1), 0.05, 'one') * 2; % Bonferroni corrected
    succShiftPval       = myBinomTest(sum(bPlShiftSucc), sum(bSuccPL), 0.05, 'one') * 2;
    failPval            = myBinomTest(sum(bShiftFail), size(encRecRankFail, 1), 0.05, 'one') * 2;
    failShiftPval       = myBinomTest(sum(bPlShiftFail), sum(bFailPL), 0.05, 'one') * 2;

    % all angles
    encAngleShift       = encAngle(bShift);
    encAnglePlShift     = encAngle(bPlShift);
    recAngleShift       = recAngle(bShift);
    recAnglePlShift     = recAngle(bPlShift);

    % successful memory angles
    encSuccAngleShift   = encSuccAngle(bShiftSucc);
    encSuccAnglePlShift = encSuccAngle(bPlShiftSucc);
    recSuccAngleShift   = recSuccAngle(bShiftSucc);
    recSuccAnglePlShift = recSuccAngle(bPlShiftSucc);

    % unsuccessful memory angles
    encFailAngleShift   = encFailAngle(bShiftFail);
    encFailAnglePlShift = encFailAngle(bPlShiftFail);
    recFailAngleShift   = recFailAngle(bShiftFail);
    recFailAnglePlShift = recFailAngle(bPlShiftFail);

    % calculate phase shifts
    circDiff            = angdiff(encAngleShift, recAngleShift);
    circDiffSucc        = angdiff(encSuccAngleShift, recSuccAngleShift);
    circDiffFail        = angdiff(encFailAngleShift, recFailAngleShift);
    circDiffSig         = angdiff(encAnglePlShift, recAnglePlShift);
    circDiffSigSucc     = angdiff(encSuccAnglePlShift, recSuccAnglePlShift);
    circDiffSigFail     = angdiff(encFailAnglePlShift, recFailAnglePlShift);

    %% permutation test for different number of significant units between successful and unsuccessful performance
    
    % difference in shifting units
    shiftDiff           = (sum(bShiftSucc) / size(bShiftSucc, 1)) - (sum(bShiftFail) / size(bShiftSucc, 1));
    shiftPlDiff         = (sum(bPlShiftSucc) / size(bOnlyPlShiftSucc, 1)) - (sum(bPlShiftFail) / size(bOnlyPlShiftFail, 1));
    
    % class variable
    succFailClass       = [ones(size(bShiftSucc)); ones(size(bShiftFail)) * 2];
    succFailPlClass     = [ones(size(bOnlyPlShiftSucc)); ones(size(bOnlyPlShiftFail)) * 2];
    
    % create surrogates
    surShiftDiff        = nan(param.nSur, 1);
    surShiftPlDiff      = nan(param.nSur, 1);
    for iSur = 1:param.nSur
        
        % concatenate successful and unsuccessful logicals
        bCatSuccFail            = cat(1, bShiftSucc, bShiftFail);
        bCatPlSuccFail          = cat(1, bOnlyPlShiftSucc, bOnlyPlShiftFail);

        % randomly permute class variable
        surSuccFailClass        = succFailClass(randperm(size(bCatSuccFail, 1)));
        surPlSuccFailClass      = succFailPlClass(randperm(size(bCatPlSuccFail, 1)));

        % assign shift indices to surrogate class variables
        surSucc                 = bCatSuccFail(surSuccFailClass == 1);
        surFail                 = bCatSuccFail(surSuccFailClass == 2);
        surPlSucc               = bCatPlSuccFail(surPlSuccFailClass == 1);
        surPlFail               = bCatPlSuccFail(surPlSuccFailClass == 2);

        % difference in shifting units
        surShiftDiff(iSur, 1)   = (sum(surSucc) / size(bShiftSucc, 1)) - (sum(surFail) / size(bShiftSucc, 1));
        surShiftPlDiff(iSur, 1) = (sum(surPlSucc) / size(bOnlyPlShiftSucc, 1)) - (sum(surPlFail) / size(bOnlyPlShiftFail, 1));
    end

    % rank
    permShiftRank       = sum(shiftDiff > surShiftDiff) / param.nSur;
    permShiftPlRank     = sum(shiftPlDiff > surShiftPlDiff) / param.nSur;

    % p-value
    permShiftPval       = 1 - permShiftRank;
    permShiftPlPval     = 1 - permShiftPlRank;

    %% permutation test for different percentage of positive shifts

    % percentage of units with positive shift
    posShifts           = sum((circDiffSucc > 0)) / size(circDiffSucc, 1);
    posShiftsPl         = sum((circDiffSigSucc > 0)) / size(circDiffSigSucc, 1);

    % class variable
    succFailClass       = [ones(size(circDiffSucc)); ones(size(circDiffFail)) * 2];
    succFailPlClass     = [ones(size(circDiffSigSucc)); ones(size(circDiffSigFail)) * 2];

    % create surrogates
    surPosShifts        = nan(param.nSur, 1);
    surPosShiftsPl      = nan(param.nSur, 1);
    for iSur = 1:param.nSur

        % concatenate successful and unsuccessful angle differences
        catCircDiffSuccFail     = cat(1, circDiffSucc, circDiffFail);
        catCircDiffSigSuccFail  = cat(1, circDiffSigSucc, circDiffSigFail);

        % randomly permute class variable
        surSuccFailClass        = succFailClass(randperm(size(catCircDiffSuccFail, 1)));
        surPlSuccFailClass      = succFailPlClass(randperm(size(catCircDiffSigSuccFail, 1)));

        % surrogate angle differences
        surCircDiffSucc          = catCircDiffSuccFail(surSuccFailClass == 1);
        surCircDiffFail          = catCircDiffSuccFail(surSuccFailClass == 2);
        surCircDiffSigSucc       = catCircDiffSigSuccFail(surPlSuccFailClass == 1);
        surCircDiffSigFail       = catCircDiffSigSuccFail(surPlSuccFailClass == 2);

        % percentage of units with positive shift
        surPosShifts(iSur, 1)   = sum((surCircDiffSucc > 0)) / size(surCircDiffSucc, 1);
        surPosShiftsPl(iSur, 1) = sum((surCircDiffSigSucc > 0)) / size(surCircDiffSigSucc, 1);
    end

    % rank
    permPosShiftRank    = sum(posShifts > surPosShifts) / param.nSur;
    permPosShiftPlRank  = sum(posShiftsPl > surPosShiftsPl) / param.nSur;

    % p-value
    permPosShiftPval    = 1 - permPosShiftRank;
    permPosShiftPlPval  = 1 - permPosShiftPlRank;

    %% figure of mean angle differences of cells with significant phase difference

    % all
    circDistSig         = figure;
    if strcmp(param.plotVersion, 'line')
        TG_PolarplotShifts_20241015(encAngleShift, recAngleShift, [0, 0.447, 0.741], [0.85, 0.325, 0.098], [0.16, 0.16, 0.16, 0.2], param);
        shiftLines                      = findall(gca, 'Type', 'line');
        for iLine = 1:length(shiftLines)
            bSigLines = bPlShift(bShift);
            if bSigLines(iLine)
                shiftLines(iLine).Color = [0.16, 0.16, 0.16, 1];
            end
        end
    elseif strcmp(param.plotVersion, 'histogram')
        circDiffHist        = polarhistogram(circDiff, param.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
        hold on;
        circDiffSigHist     = polarhistogram(circDiffSig, param.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16]);
        rlim([0, 10]);
        set(gca, 'ThetaTickLabel', param.polarHistLabels);
    end
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bPlShift)), 32, 'out of', 32, num2str(sum(bBothPL)), ')'));

    % print figure
    print(circDistSig, fullfile(paths.save, strcat(param.resultName, '_BothPhaseDifference')), '-dsvg', '-r300');

    % successful memory performance
    circDistSigSucc     = figure;
    if strcmp(param.plotVersion, 'line')
        TG_PolarplotShifts_20241015(encSuccAngleShift, recSuccAngleShift, [0, 0.447, 0.741], [0.85, 0.325, 0.098], [0.1, 0.6, 0.1, 0.2], param);
        shiftLinesSucc      = findall(gca, 'Type', 'line');
        for iLine = 1:length(shiftLinesSucc)
            bSigLinesSucc = bPlShiftSucc(bShiftSucc);
            if bSigLinesSucc(iLine)
                shiftLinesSucc(iLine).Color = [0.1, 0.6, 0.1, 1];
            end
        end
    elseif strcmp(param.plotVersion, 'histogram')
        circDiffHistSucc    = polarhistogram(circDiffSucc, param.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1], 'FaceAlpha', 0.2);
        hold on;
        circDiffSigHistSucc = polarhistogram(circDiffSigSucc, param.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1]);
        rlim([0, 10]);
        set(gca, 'ThetaTickLabel', param.polarHistLabels);
    end
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bPlShiftSucc)), 32, 'out of', 32, num2str(sum(bSuccPL)), ')'));

    % print figure
    print(circDistSigSucc, fullfile(paths.save, strcat(param.resultName, '_SuccPhaseDifference')), '-dsvg', '-r300');

    % unsuccessful memory performance
    circDistSigFail     = figure;
    if strcmp(param.plotVersion, 'line')
        TG_PolarplotShifts_20241015(encFailAngleShift, recFailAngleShift, [0, 0.447, 0.741], [0.85, 0.325, 0.098], [0.6, 0.1, 0.1, 0.2], param);
        shiftLinesFail      = findall(gca, 'Type', 'line');
        for iLine = 1:length(shiftLinesFail)
            bSigLinesFail = bPlShiftFail(bShiftFail);
            if bSigLinesFail(iLine)
                shiftLinesFail(iLine).Color = [0.6, 0.1, 0.1, 1];
            end
        end
    elseif strcmp(param.plotVersion, 'histogram')
        circDiffHistFail    = polarhistogram(circDiffFail, param.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1], 'FaceAlpha', 0.2);
        hold on;
        circDiffSigHistFail = polarhistogram(circDiffSigFail, param.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1]);
        rlim([0, 10]);
        set(gca, 'ThetaTickLabel', param.polarHistLabels);
    end
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bPlShiftFail)), 32, 'out of', 32, num2str(sum(bFailPL)), ')'));

    % print figure
    print(circDistSigFail, fullfile(paths.save, strcat(param.resultName, '_FailPhaseDifference')), '-dsvg', '-r300');

    %% circular mean differences

    % circular means
    meanCircDiff        = rad2deg(circ_mean(abs(circDiff)));
    meanCircDiffSig     = rad2deg(circ_mean(abs(circDiffSig)));
    meanCircDiffSucc    = rad2deg(circ_mean(abs(circDiffSucc)));
    meanCircDiffSigSucc = rad2deg(circ_mean(abs(circDiffSigSucc)));
    meanCircDiffFail    = rad2deg(circ_mean(abs(circDiffFail)));
    meanCircDiffSigFail = rad2deg(circ_mean(abs(circDiffSigFail)));
    
    % standard deviations
    stdCircDiff         = rad2deg(circ_std(abs(circDiff)));
    stdCircDiffSig      = rad2deg(circ_std(abs(circDiffSig)));
    stdCircDiffSucc     = rad2deg(circ_std(abs(circDiffSucc)));
    stdCircDiffSigSucc  = rad2deg(circ_std(abs(circDiffSigSucc)));
    stdCircDiffFail     = rad2deg(circ_std(abs(circDiffFail)));
    stdCircDiffSigFail  = rad2deg(circ_std(abs(circDiffSigFail)));

    %% bar plot

    % define figure for the bar plot
    shiftFig        = figure;

    % define categories for the x-axis
    shiftBounds     = categorical({'all shifts', 'PL shifts', 'succ shifts', 'succ PL shifts', 'fail shifts', 'fail PL shifts'});

    % define shiftBounds with explicit order
    shiftBounds     = categorical(shiftBounds, shiftBounds);

    % reorder the categories according to the desired sequence
    desiredOrder    = {'all shifts', 'PL shifts', 'succ shifts', 'succ PL shifts', 'fail shifts', 'fail PL shifts'};
    shiftBounds     = reordercats(shiftBounds, desiredOrder);

    % calculate ratio for each category
    shiftRatio      = [
        (sum(bShift) / size(bShift, 1)), ...
        (sum(bPlShift) / sum(bBothPL)), ...
        (sum(bShiftSucc) / size(bShiftSucc, 1)), ...
        (sum(bPlShiftSucc) / sum(bSuccPL)), ...
        (sum(bShiftFail) / size(bShiftFail, 1)), ...
        (sum(bPlShiftFail) / sum(bFailPL))
        ];

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
    b = gobjects(size(shiftBounds));

    % clear current axes and begin plotting
    cla();
    hold on;

    % loop through each category to create bars with corresponding properties
    for iLine = 1:numel(shiftBounds)
        b(iLine) = bar(shiftBounds(iLine), shiftRatio(iLine));
        set(b(iLine), 'FaceColor', shiftColors(iLine, :), 'FaceAlpha', shiftAlphas(iLine));
    end

    % adjustments
    ylabel('Ratio');
    ylim([0, 0.1801]);
    yticks(0:0.05:1);
    set(gca, 'TickDir', 'out');

    % print figure
    print(shiftFig, fullfile(paths.save, 'allUnits_allCells_ShiftCellPercentage'), '-dsvg', '-r300');
else
    
    % phase locking and phase difference
    bShift              = encRecRank > 0.95;

    % number of included cells
    numInc              = sum(exclusionIndex == 0);

    % binomial tests
    allPval             = myBinomTest(sum(bShift), numInc, 0.05, 'one') * 2; % Bonferroni correction

    %% figure of mean angle differences of cells with significant phase difference

    % all
    circDistSig         = figure;
    encAngleShift       = encAngle(bShift);
    recAngleShift       = recAngle(bShift);
    circDiff            = angdiff(encAngle(bShift), recAngle(bShift));
    if strcmp(param.plotVersion, 'line')
        TG_PolarplotShifts_20241015(encAngleShift, recAngleShift, [0, 0.447, 0.741], [0.85, 0.325, 0.098], 'k', param);
    elseif strcmp(param.plotVersion, 'histogram')
        circDiffHist    = polarhistogram(circDiff, param.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
        rlim([0, 10]);
        set(gca, 'ThetaTickLabel', param.polarHistLabels);
    end
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bShift)), 32, 'out of', 32, num2str(numInc), ')'));

    % print figure
    print(circDistSig, fullfile(paths.save, strcat(param.resultName, '_BothPhaseDifference')), '-dsvg', '-r300');

    %% test, if there is a specific cluster of phase shifts for upper-frequency spikes
    if strcmp(param.powLvl, 'all') && ...
        strcmp(param.slopeLvl, 'all') && ...
        strcmp(param.freqLvl, 'all') && ...
        strcmp(param.thetaLvl, 'all') && ...
        strcmp(param.unitType, 'all') && ...
        strcmp(param.cellType, 'all')

        % create surrogate phase shifts
        surCircDiff     = 2 * pi * rand(size(circDiff, 1), param.nSur) - pi;

        % empirical PPC
        shiftPPC        = TG_PPC_20241128(circDiff);

        % surrogate PPCs
        surShiftPPC     = nan(param.nSur, 1);
        for iSur = 1:param.nSur
            surShiftPPC(iSur, 1) = TG_PPC_20241128(surCircDiff(:, iSur));
        end

        % p-value
        shiftClusPval   = 1 - (sum(shiftPPC > surShiftPPC) / param.nSur);
    end
end