%==========================================================================
% This script quantifies the occurrence of SPRiNT oscillation peaks and
% the aperiodic slope steepness during baseline, encoding and retrieval.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths               = struct();
paths.behlog        = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.sprint        = 'D:\TreasureHunt\MicroSPRiNT_20230713'; % folder with SPRiNT results
paths.bycycle       = 'D:\TreasureHunt\MicroBycycle_20241216'; % folder with Bycycle results
paths.artifact      = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.phaseData     = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.SUData        = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % folder with single unit data
paths.phaseRes	    = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_Results'; % phase analysis folder
paths.save    	    = 'D:\TreasureHunt\SPRiNTAnalysis_20231212'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

%% set configurations

% parameters
params              = struct();
params.filterBand   = [1, 10];
params.sprintName   = 'datacutSprint.mat';
params.nSur         = 10001;
params.nSub         = 100;
params.surHistEdges = linspace(-4, 4, 201);

% set random seed
randSeed    = 444;
rng(randSeed, 'twister');
randSeedNum = randi(100000, params.nSur, 1);

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


%% load phase results to get list of included wires
phaseRes            = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

%% wire index
unitIdx             = cat(1, phaseRes.phaseRes.idx);
[wireIdx, selIdx]   = unique(unitIdx(:, 1:3), 'rows'); % all included wires

%% loop through wires
encTheta            = nan(size(wireIdx, 1), 1);
recTheta            = nan(size(wireIdx, 1), 1);
blnTheta            = nan(size(wireIdx, 1), 1);
encSlope            = nan(size(wireIdx, 1), 1);
recSlope            = nan(size(wireIdx, 1), 1);
blnSlope            = nan(size(wireIdx, 1), 1);
parfor iWire = 1:size(wireIdx, 1)
    
    % report progress
    disp(iWire);
    
    %% get data

    % get path of sprint data
    subName             = subjects{wireIdx(iWire, 1)};
    sessName            = strcat('session_', num2str(wireIdx(iWire, 2)));
    chanDir             = TG_GetChanDir_20210812(paths.SUData, subName, sessName);
    
    % load sprint data
    sprintData          = load(fullfile(paths.sprint, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.sprintName));

    % load bycycle data
    bycycleData         = load(fullfile(paths.bycycle, subName, sessName, chanDir(wireIdx(iWire, 3)).name, 'datacutFt2000HzBP_bycycle_1_10.mat'));
    
    % load data for artifact removal
    artifactData        = load(fullfile(paths.artifact, subName, sessName, chanDir(wireIdx(iWire, 3)).name, 'detectedArtifacts.mat'), 'bArtifact');
    
    % load sampling rate of artifact index
    phaseData           = load(fullfile(paths.phaseData, subName, sessName, chanDir(wireIdx(iWire, 3)).name, 'datacutFt2000HzBP_generalized_1_10.mat'), 'fsample', 'sampleinfo');

    % load segment information
    segmentInfo         = load(fullfile(paths.behlog, subName, sessName, 'segmentInfo_20230208.mat'));

    %% extract slopes

    % SPRiNT times
    sprintTime          = cat(2, sprintData.SPRiNT.channel.aperiodics.time);

    % overall time coverage of SPRiNT windows
    sprintWindow        = sprintData.SPRiNT.options.WinLength + (sprintData.SPRiNT.options.WinAverage - 1) * (sprintData.SPRiNT.options.WinLength * (1 - (sprintData.SPRiNT.options.WinOverlap / 100)));

    % SPRiNT times in samples
    sprintSamples       = sprintTime * phaseData.fsample;

    % samples of overall SPRiNT windows
    sprintWindowsSmp    = sprintSamples' + ((-(sprintWindow / 2 * phaseData.fsample) + 1):1:((sprintWindow / 2 * phaseData.fsample)));

    % find artifacts
    sprintBArtifact     = any(artifactData.bArtifact(sprintWindowsSmp), 2)';

    % extract slopes
    slopeOrig           = cat(2, sprintData.SPRiNT.channel.aperiodics.exponent);

    % upsample slope and artifact index
    sprintTimeWindow                = mode(diff(sprintTime));
    sprintFsample                   = 1 / sprintTimeWindow;
    upsamplingRatio                 = phaseData.fsample / sprintFsample; % upsampling factor
    slope                           = repelem(slopeOrig, upsamplingRatio);
    sprintArtifactIdx               = repelem(sprintBArtifact, upsamplingRatio);

    % add NaNs to SPRiNT data to start and end borders for padding
    missingSprint                   = size(artifactData.bArtifact, 2) - size(slope, 2);
    slope                           = cat(2, repelem(NaN, missingSprint / 2), slope, repelem(NaN, missingSprint / 2));
    sprintArtifactIdx               = cat(2, repelem(1, missingSprint / 2), sprintArtifactIdx, repelem(1, missingSprint / 2));

    % extract theta index
    thetaIdx                        = double(bycycleData.bBurst);

    % remove artifacts
    thetaIdx(artifactData.bArtifact == 1)   = NaN;
    slope(sprintArtifactIdx == 1)           = NaN;

    %% calculate oscillation percentage and slope during encoding

    % encoding period
    encodingSamples                 = round(segmentInfo.encoding ./ (segmentInfo.fsample / phaseData.fsample));
    encIdx                          = [];
    for iEnc = 1:size(encodingSamples, 1)
        encIdx      = cat(2, encIdx, encodingSamples(iEnc, 1):encodingSamples(iEnc, 2));
    end

    % oscillation percentage
    thetaIdxEnc                     = thetaIdx(encIdx);
    thetaIdxEnc(isnan(thetaIdxEnc)) = [];
    encTheta(iWire, 1)              = sum(thetaIdxEnc) / size(thetaIdxEnc, 2);

    % slope
    slopeEnc                        = slope(encIdx);
    encSlope(iWire, 1)              = mean(slopeEnc, 'omitnan');

    %% calculate oscillation percentage and slope during recall

    % recall index
    recallSamples                   = round(cat(1, segmentInfo.objRecall, segmentInfo.locRecall) ./ (segmentInfo.fsample / phaseData.fsample));
    allRecIdx                       = [];
    for iRec = 1:size(recallSamples, 1)
        allRecIdx   = cat(2, allRecIdx, recallSamples(iRec, 1):recallSamples(iRec, 2));
    end

    % subsample to encoding period and compute oscillation percentage and slope
    rng(randSeedNum(iWire));
    recThetaSub     = nan(params.nSub, 1);
    recSlopeSub     = nan(params.nSub, 1);
    for iSub = 1:params.nSub

        % recall index
        recIdx                          = allRecIdx(randperm(length(allRecIdx), length(encIdx)));

        % oscillation
        thetaIdxRec                     = thetaIdx(recIdx);
        thetaIdxRec(isnan(thetaIdxRec)) = [];
        recThetaSub(iSub, 1)            = sum(thetaIdxRec) / size(thetaIdxRec, 2);

        % slope
        slopeRec                        = slope(recIdx);
        recSlopeSub(iSub, 1)            = mean(slopeRec, 'omitnan');
    end

    % mean
    recTheta(iWire, 1)      = mean(recThetaSub);
    recSlope(iWire, 1)      = mean(recSlopeSub);

    %% calculate oscillation percentage and slope during baseline

    % baseline index
    rng(randSeedNum(iWire));
    bBln                            = ones(phaseData.sampleinfo);
    bBln(allRecIdx)                 = 0;
    bBln(encIdx)                    = 0;
    allBlnIdx                       = find(bBln);
    
    % subsample to encoding period
    rng(randSeedNum(iWire));
    blnThetaSub     = nan(params.nSub, 1);
    blnSlopeSub     = nan(params.nSub, 1);
    for iSub = 1:params.nSub

        % baseline index
        blnIdx                          = allBlnIdx(randperm(length(allBlnIdx), length(encIdx)));

        % oscillation
        thetaIdxBln                     = thetaIdx(blnIdx);
        thetaIdxBln(isnan(thetaIdxBln)) = [];
        blnThetaSub(iSub, 1)            = sum(thetaIdxBln) / size(thetaIdxBln, 2);

        % slope
        slopeBln                        = slope(blnIdx);
        blnSlopeSub(iSub, 1)            = mean(slopeBln, 'omitnan');
    end

    % mean
    blnTheta(iWire, 1)      = mean(blnThetaSub);
    blnSlope(iWire, 1)      = mean(blnSlopeSub);
end

%% average over sessions

% session index
[sessNum, ~, sessIdx] = unique(wireIdx(:, 1:2), 'row');

% loop through sessions
thisSessEncTheta    = nan(size(sessNum, 1), 1);
thisSessRecTheta    = nan(size(sessNum, 1), 1);
thisSessBlnTheta    = nan(size(sessNum, 1), 1);
thisSessEncSlope    = nan(size(sessNum, 1), 1);
thisSessRecSlope    = nan(size(sessNum, 1), 1);
thisSessBlnSlope    = nan(size(sessNum, 1), 1);
for iSess = 1:size(sessNum, 1)

    % get mean theta values for this session
    thisSessEncTheta(iSess, 1)      = mean(encTheta(sessIdx == iSess));
    thisSessRecTheta(iSess, 1)      = mean(recTheta(sessIdx == iSess));
    thisSessBlnTheta(iSess, 1)      = mean(blnTheta(sessIdx == iSess));

    % get mean slope values for this session
    thisSessEncSlope(iSess, 1)      = mean(encSlope(sessIdx == iSess));
    thisSessRecSlope(iSess, 1)      = mean(recSlope(sessIdx == iSess));
    thisSessBlnSlope(iSess, 1)      = mean(blnSlope(sessIdx == iSess));
end

%% means and standard erros of the mean

% theta means
meanEncTheta    = mean(thisSessEncTheta);
meanRecTheta    = mean(thisSessRecTheta);
meanBlnTheta    = mean(thisSessBlnTheta);

% theta sem
semEncTheta     = std(thisSessEncTheta) / sqrt(length(thisSessEncTheta));
semRecTheta     = std(thisSessRecTheta) / sqrt(length(thisSessRecTheta));
semBlnTheta     = std(thisSessBlnTheta) / sqrt(length(thisSessBlnTheta));

% slope means
meanEncSlope    = mean(thisSessEncSlope);
meanRecSlope    = mean(thisSessRecSlope);
meanBlnSlope    = mean(thisSessBlnSlope);

% slope sem
semEncSlope     = std(thisSessEncSlope) / sqrt(length(thisSessEncSlope));
semRecSlope     = std(thisSessRecSlope) / sqrt(length(thisSessRecSlope));
semBlnSlope     = std(thisSessBlnSlope) / sqrt(length(thisSessBlnSlope));

%% statistical tests

% theta
[encRecThetaRank, encRecThetaT, encRecThetaSurT]     = TG_PermTest1D_2PS_20230727(thisSessEncTheta, thisSessRecTheta, params.nSur, randSeedNum);
[encBlnThetaRank, encBlnThetaT, encBlnThetaSurT]     = TG_PermTest1D_2PS_20230727(thisSessEncTheta, thisSessBlnTheta, params.nSur, randSeedNum);
[recBlnThetaRank, recBlnThetaT, recBlnThetaSurT]     = TG_PermTest1D_2PS_20230727(thisSessRecTheta, thisSessBlnTheta, params.nSur, randSeedNum);

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

% slope
[encRecSlopeRank, encRecSlopeT, encRecSlopeSurT]     = TG_PermTest1D_2PS_20230727(thisSessEncSlope, thisSessRecSlope, params.nSur, randSeedNum);
[encBlnSlopeRank, encBlnSlopeT, encBlnSlopeSurT]     = TG_PermTest1D_2PS_20230727(thisSessEncSlope, thisSessBlnSlope, params.nSur, randSeedNum);
[recBlnSlopeRank, recBlnSlopeT, recBlnSlopeSurT]     = TG_PermTest1D_2PS_20230727(thisSessRecSlope, thisSessBlnSlope, params.nSur, randSeedNum);

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

%% plots

% theta figure
thetaFig        = figure;
bTheta          = bar([meanEncTheta, meanRecTheta, meanBlnTheta] * 100, 'FaceColor', 'flat');
bTheta.CData    = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.6, 0.4, 0.2];
hold on;
errorbar([meanEncTheta, meanRecTheta, meanBlnTheta] * 100, [semEncTheta, semRecTheta, semBlnTheta] * 100, 'LineStyle', 'none', 'Color', 'k');
set(gca, 'XTickLabel', {'Encoding', 'Recall', 'Baseline'});
set(gca,'TickDir', 'out');
box off;
ylim([0, 50]);
ylabel('Oscillation detected (%)');
set(gcf, 'Renderer', 'painters');
saveas(thetaFig, fullfile(paths.save, 'blnEncRecOscillationPercentage.svg'));

% slope figure
slopeFig        = figure;
bSlope          = bar([meanEncSlope, meanRecSlope, meanBlnSlope], 'FaceColor', 'flat');
bSlope.CData    = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.6, 0.4, 0.2];
hold on;
errorbar([meanEncSlope, meanRecSlope, meanBlnSlope], [semEncSlope, semRecSlope, semBlnSlope], 'LineStyle', 'none', 'Color', 'k');
set(gca, 'XTickLabel', {'Encoding', 'Recall', 'Baseline'});
set(gca,'TickDir', 'out');
box off;
ylim([0, 2]);
yticks([0, 1, 2]);
ylabel('Aperiodic slope');
set(gcf, 'Renderer', 'painters');
saveas(slopeFig, fullfile(paths.save, 'blnEncRecSlope.svg'));

%% plot theta surrogate distributions

% encoding versus recall theta
encRecThetaFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encRecThetaTCounts  = histcounts(encRecThetaSurT, params.surHistEdges);
encRecThetaHist     = histogram('BinCounts', encRecThetaTCounts, 'BinEdges', params.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encRecThetaT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encRecThetaFigure, fullfile(paths.save, 'encRecThetaSurDistr.svg'));

% encoding versus baseline theta
encBlnThetaFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encBlnThetaTCounts  = histcounts(encBlnThetaSurT, params.surHistEdges);
encBlnThetaHist     = histogram('BinCounts', encBlnThetaTCounts, 'BinEdges', params.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encBlnThetaT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encBlnThetaFigure, fullfile(paths.save, 'encBlnThetaSurDistr.svg'));

% recall versus baseline theta
recBlnThetaFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
recBlnThetaTCounts  = histcounts(recBlnThetaSurT, params.surHistEdges);
recBlnThetaHist     = histogram('BinCounts', recBlnThetaTCounts, 'BinEdges', params.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(recBlnThetaT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(recBlnThetaFigure, fullfile(paths.save, 'recBlnThetaSurDistr.svg'));

%% plot slope surrogate distributions

% encoding versus recall slope
encRecSlopeFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encRecSlopeTCounts  = histcounts(encRecSlopeSurT, params.surHistEdges);
encRecSlopeHist     = histogram('BinCounts', encRecSlopeTCounts, 'BinEdges', params.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encRecSlopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encRecSlopeFigure, fullfile(paths.save, 'encRecSlopeSurDistr.svg'));

% encoding versus baseline slope
encBlnSlopeFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
encBlnSlopeTCounts  = histcounts(encBlnSlopeSurT, params.surHistEdges);
encBlnSlopeHist     = histogram('BinCounts', encBlnSlopeTCounts, 'BinEdges', params.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(encBlnSlopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(encBlnSlopeFigure, fullfile(paths.save, 'encBlnSlopeSurDistr.svg'));

% recall versus baseline slope
recBlnSlopeFigure   = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
recBlnSlopeTCounts  = histcounts(recBlnSlopeSurT, params.surHistEdges);
recBlnSlopeHist     = histogram('BinCounts', recBlnSlopeTCounts, 'BinEdges', params.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(recBlnSlopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-4, 4]);
xticks([-4, 0, 4]);
ylim([0, 300]);
yticks([0, 300]);
set(gca,'TickDir','out');
box off;
saveas(recBlnSlopeFigure, fullfile(paths.save, 'recBlnSlopeSurDistr.svg'));

