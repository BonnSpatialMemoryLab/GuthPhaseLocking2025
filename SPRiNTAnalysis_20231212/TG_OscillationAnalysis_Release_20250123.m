%==========================================================================
% This script analyzes the occurrence and frequency of oscillations
% between 1 and 10 Hz during the whole experiment.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.bycycle   = 'D:\TreasureHunt\MicroBycycle_20241216'; % folder with Bycycle results
paths.artifact  = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.phaseData = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.SUData    = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % folder with single unit data
paths.phaseRes	= 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_Results'; % phase analysis folder
paths.save    	= 'D:\TreasureHunt\SPRiNTAnalysis_20231212'; % save folder
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
params                 = [];
params.percEdges       = linspace(0, 1, 11);
params.binEdges        = linspace(1, 10, 19);
params.binCenters      = (params.binEdges(1 : end - 1) + params.binEdges(2 : end)) / 2;
params.wireEdges       = linspace(1, 10, 10);
params.filterBand      = [1, 10];
params.bycycleName     = 'datacutFt2000HzBP_bycycle_1_10.mat';

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

%% brain region

% brain region for all included units
brainRegIdx          = {phaseRes.phaseRes.brainRegion}';
brainRegIdxSplit     = split(brainRegIdx, '_');

% brain region for all wires
brainRegWireIdx      = brainRegIdxSplit(selIdx, 1);

% names of unique brain regions
uniqueBrainReg       = unique(brainRegIdxSplit(:, 1));

%% loop through wires

% define variables
burstPerc               = nan(size(wireIdx(:, 1)));
burstFreq               = cell(size(wireIdx(:, 1)));
cycleFreq               = cell(size(wireIdx(:, 1)));
thetaDuringBurstsPerc   = nan(size(wireIdx(:, 1)));
thetaAbsenceBurstsPerc  = nan(size(wireIdx(:, 1)));

% loop through wires
parfor iWire = 1:size(wireIdx, 1)
    
    % report progress
    disp(iWire);
    
    % get path of data
    subName             = subjects{wireIdx(iWire, 1)};
    sessName            = strcat('session_', num2str(wireIdx(iWire, 2)));
    chanDir             = TG_GetChanDir_20210812(paths.SUData, subName, sessName);

    % load bycycle data
    bycycleData         = load(fullfile(paths.bycycle, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.bycycleName));
    
    % load data for artifact removal
    artifactData        = load(fullfile(paths.artifact, subName, sessName, chanDir(wireIdx(iWire, 3)).name, 'detectedArtifacts.mat'), 'bArtifact');
    
    % load sampling rate of artifact index
    phaseData           = load(fullfile(paths.phaseData, subName, sessName, chanDir(wireIdx(iWire, 3)).name, 'datacutFt2000HzBP_generalized_1_10.mat'), 'fsample', 'sampleinfo');

    %% cycle frequencies

    % cycle start and end times
    cycleStart          = bycycleData.bycycleTable.sample_last_peak;
    cycleEnd            = bycycleData.bycycleTable.sample_next_peak;

    % cycle index
    cycleIdx                                = zeros(phaseData.sampleinfo);
    cycleIdx([cycleStart; cycleEnd(end)])   = 1;
    cycleIdx                                = cumsum(cycleIdx);
    cycleIdx(cycleEnd(end):end)             = NaN;

    % find cycles during artifacts (have to be removed entirely)
    cycleArtifactIdx    = unique(cycleIdx(artifactData.bArtifact))';
    cycleArtifactIdx    = cycleArtifactIdx(~isnan(cycleArtifactIdx) & cycleArtifactIdx > 0);

    % remove from cycle start and end arrays
    cycleStart(cycleArtifactIdx)    = [];
    cycleEnd(cycleArtifactIdx)      = [];

    % cycle frequency
    cycleFreq{iWire, 1} = phaseData.fsample ./ (cycleEnd - cycleStart);

    %% burst frequencies

    % burst index
    burstIdx            = double(bycycleData.bBurst);
    burstIdxLabeled     = bwlabel(burstIdx);

    % find bursts during artifacts (have to be removed entirely)
    artifactIdx         = unique(burstIdxLabeled(artifactData.bArtifact));
    artifactIdx         = artifactIdx(artifactIdx > 0);

    % adjust burst index
    burstIdx(ismember(burstIdxLabeled, artifactIdx))        = NaN;
    burstIdxLabeled(ismember(burstIdxLabeled, artifactIdx)) = NaN;

    % included bursts
    includedBurstIdx    = unique(burstIdxLabeled(burstIdxLabeled > 0));

    % burst information
    burstFreq{iWire, 1} = bycycleData.burstFrequency(includedBurstIdx);

    %% burst percentages

    % percentage
    burstPerc(iWire, 1) = sum(burstIdx(~isnan(burstIdx) & ~artifactData.bArtifact)) / size(burstIdx(~isnan(burstIdx) & ~artifactData.bArtifact), 2);
end

%% results for different brain regions

% add 'all' to loop
brainRegForLoop = cat(1, {'ALL'}, uniqueBrainReg);

% loop through brain regions
for iReg = 1:size(brainRegForLoop, 1)

    %  selection index for specific regions
    if strcmp(brainRegForLoop{iReg}, 'ALL')
        bSel = true(size(brainRegWireIdx));
    else
        bSel = strcmp(brainRegForLoop{iReg}, brainRegWireIdx);
    end
    
    % only include brain areas with data from at least 5 sessions
    subAndSess          = wireIdx(bSel, 1:2);
    uniqueSess          = unique(subAndSess, 'rows');
    if size(uniqueSess, 1) < 5
        continue;
    end

    % select results
    thisRegCycleFreq    = cycleFreq(bSel);
    thisRegBurstPerc    = burstPerc(bSel);
    thisRegBurstFreq    = burstFreq(bSel);

    %% print results for cycles

    % print histogram of theta frequencies for all cycles (normalized)
    binBurstFig = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    burstCounts = cell2mat(cellfun(@(x) histcounts(x, params.binEdges, 'Normalization', 'probability'), thisRegCycleFreq, 'Uni', 0));
    TG_ShadeSEM_20210714(params.binCenters, burstCounts * 100, 'k', 0.5);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Frequency (Hz)');
    ylabel('Probability (%)');
    xticks(0:2:10);
    ylim([0, 25]);
    grid on;
    title(strcat(brainRegForLoop{iReg}, 32, num2str(size(cat(1, thisRegCycleFreq{:}), 1)), 32, 'cycles'));
    saveas(binBurstFig, fullfile(paths.save, strcat(brainRegForLoop{iReg}, '_FrequenciesOfCycles.svg')));

    %% print results for bursts

    % print histogram of burst percentages
    prctBurstFig    = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    percCounts      = histcounts(thisRegBurstPerc, params.percEdges);
    histogram('BinCounts', percCounts, 'BinEdges', params.percEdges * 100, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Burst presence (%)');
    ylabel('Number of wires');
    ylim([0, 400]);
    title(strcat(brainRegForLoop{iReg}, 32, num2str(sum(bSel)), 32, 'wires'));
    saveas(prctBurstFig, fullfile(paths.save, strcat(brainRegForLoop{iReg}, '_RatioOscillatoryBurst.svg')));

    % print histogram of theta frequencies for all wires
    wireBurstFig    = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    wireBurstFreq   = cell2mat(cellfun(@(x) mode(round(x)), thisRegBurstFreq, 'Uni', 0));
    wireBurstCounts = histcounts(wireBurstFreq, params.wireEdges);
    histogram('BinCounts', wireBurstCounts, 'BinEdges', params.wireEdges, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Frequency (Hz)');
    ylabel('Number of wires');
    title(strcat(brainRegForLoop{iReg}, 32, num2str(sum(bSel)), 32, 'wires'));
    saveas(wireBurstFig, fullfile(paths.save, strcat(brainRegForLoop{iReg}, '_BurstFrequenciesWires.svg')));

    % print histogram of theta frequencies for all bursts
    binBurstFig = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    burstCounts = cell2mat(cellfun(@(x) histcounts(x, params.binEdges, 'Normalization', 'probability'), thisRegBurstFreq, 'Uni', 0));
    TG_ShadeSEM_20210714(params.binCenters, burstCounts * 100, 'k', 0.5);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Frequency (Hz)');
    ylabel('Probability (%)');
    xticks(0:2:10);
    ylim([0, 25]);
    grid on;
    title(strcat(brainRegForLoop{iReg}, 32, num2str(size(cat(1, thisRegBurstFreq{:}), 1)), 32, 'oscillations'));
    saveas(binBurstFig, fullfile(paths.save, strcat(brainRegForLoop{iReg}, '_FrequenciesOfBursts.svg')));
end