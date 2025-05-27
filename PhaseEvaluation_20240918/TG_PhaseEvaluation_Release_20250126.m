%==========================================================================
% This script checks, if phase computations are valid. It also computes
% the cycle-by-cycle frequency of the generalized phase estimates and
% analyzes these frequencies during the distractor task and during times of
% different aperiodic slopes.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% set random seed
rng(444);

% paths
paths.filteredData          = 'D:\TreasureHunt\MicroFiltered_20210930'; % filtered data
paths.phaseData             = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.SUData                = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % folder with single unit data
paths.phaseRes	            = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % phase analysis folder
paths.artifact              = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.slope                 = 'D:\TreasureHunt\MicroSPRiNT_20230713'; % slope information
paths.bycycle               = 'D:\TreasureHunt\MicroBycycle_20241216'; % oscillation information
paths.behData               = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.save                  = 'D:\TreasureHunt\PhaseEvaluation_20240918';

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% parameters
params.maxSamples           = 1.2e7;
params.polarHistLabels      = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
params.windowSize           = 10; % in seconds
params.filterBand           = [1, 10];
params.filterBandStr        = strjoin(strsplit(num2str(params.filterBand)), '_');
params.lowThetaBand         = [2, 5];
params.highThetaBand        = [6, 9];
params.filteredDataName     = strcat('datacutFt2000HzBP_', params.filterBandStr, '.mat');
params.phaseDataName        = strcat('datacutFt2000HzBP_generalized_', params.filterBandStr, '.mat');
params.artifactName         = 'detectedArtifacts.mat';
params.slopeName            = 'datacutSprint.mat';
params.bycycleName          = 'datacutFt2000HzBP_bycycle_1_10.mat';
params.segmentInfoName      = 'segmentInfo_20230208.mat';
params.lowThetaDataName     = 'datacutFt2000HzBP_hilbert_2_5.mat';
params.highThetaDataName    = 'datacutFt2000HzBP_hilbert_6_9.mat';
params.phaseColor           = [0.3, 0.3, 0.3, 0.5];
params.heatmapEdges         = -1000:1000;
params.heatmapCenters       = (params.heatmapEdges(1 : end - 1) + params.heatmapEdges(2 : end)) / 2;

% parameters for Fourier transform
params.FouCfg               = [];
params.FouCfg.method        = 'mtmfft';
params.FouCfg.output        = 'pow';
params.FouCfg.pad           = 'nextpow2';
params.FouCfg.taper         = 'hanning';
params.FouCfg.foi           = round(logspace(log10(0.01), log10(20), 100), 2);

% phase analysis folder
paths.phaseRes	            = fullfile(paths.phaseRes, strcat(params.filterBandStr, '_Hz_Results'));

% subjects
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

% load phase results to get list of included wires
phaseRes            = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

% wire index
unitIdx             = cat(1, phaseRes.phaseRes.idx);
[wireIdx, selIdx]   = unique(unitIdx(:, 1:3), 'rows'); % all included wires

% preallocate
FreqanalysisData     = cell(size(wireIdx, 1), 1);
allReal             = cell(size(wireIdx, 1), 1);
allImag             = cell(size(wireIdx, 1), 1);
allAngle            = cell(size(wireIdx, 1), 1);
allFreqs            = cell(size(wireIdx, 1), 1);
allAngleCorr        = cell(size(wireIdx, 1), 1);
lowThetaDiff        = cell(size(wireIdx, 1), 1);
highThetaDiff       = cell(size(wireIdx, 1), 1);
oscHighSlopePerc    = nan(size(wireIdx, 1), 1);
oscLowSlopePerc     = nan(size(wireIdx, 1), 1);
noOscHighSlopePerc  = nan(size(wireIdx, 1), 1);
noOscLowSlopePerc   = nan(size(wireIdx, 1), 1);
succDistrFreqs      = nan(size(wireIdx, 1), 1);
failDistrFreqs      = nan(size(wireIdx, 1), 1);
corrSlope           = nan(size(wireIdx, 1), 1);
polyFitSlopeCycle   = nan(size(wireIdx, 1), 2);
meanAnalytic        = nan(size(wireIdx, 1), 1);
stdAnalytic         = nan(size(wireIdx, 1), 1);
mrvl                = nan(size(wireIdx, 1), 1);
meanAngle           = nan(size(wireIdx, 1), 1);
meanInstFreq        = nan(size(wireIdx, 1), 1);

%% loop through wires
parfor iWire = 1:size(wireIdx, 1)

    % report progress
    disp(iWire);

    % get path of data
    subName                                 = subjects{wireIdx(iWire, 1)};
    sessName                                = strcat('session_', num2str(wireIdx(iWire, 2)));
    chanDir                                 = TG_GetChanDir_20210812(paths.SUData, subName, sessName);

    %% load data

    % load filtered data
    filteredData                            = load(fullfile(paths.filteredData, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.filteredDataName));

    % load phase data
    phaseData                               = load(fullfile(paths.phaseData, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.phaseDataName));
    thisPhaseData                           = phaseData.trial{1, 1};

    % load artifact data
    artifactData                            = load(fullfile(paths.artifact, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.artifactName));
    bArtifact                               = artifactData.bArtifact;

    % load slope data
    sprintData                              = load(fullfile(paths.slope, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.slopeName));

    % load oscillation data
    bycycleData                             = load(fullfile(paths.bycycle, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.bycycleName));

    % load information about distractor task
    behData                                 = load(fullfile(paths.behData, subName, sessName, params.segmentInfoName));
    distrTask                               = round(behData.distrTask ./ behData.fsample .* phaseData.fsample);
    distrTaskPerf                           = behData.distrTaskPerf;

    %% power spectrum of filtered data

    % compute power spectrum
    FreqanalysisData{iWire, 1}              = ft_freqanalysis(params.FouCfg, filteredData);

    %% extract slopes and theta frequencies

    % SPRiNT times
    sprintTime                              = cat(2, sprintData.SPRiNT.channel.aperiodics.time);

    % extract slopes
    slopeOrig                               = cat(2, sprintData.SPRiNT.channel.aperiodics.exponent);

    % overall time coverage of SPRiNT windows
    sprintWindow                            = sprintData.SPRiNT.options.WinLength + (sprintData.SPRiNT.options.WinAverage - 1) * (sprintData.SPRiNT.options.WinLength * (1 - (sprintData.SPRiNT.options.WinOverlap / 100)));

    % SPRiNT times in samples
    sprintSamples                           = sprintTime * phaseData.fsample;

    % samples of overall SPRiNT windows
    sprintWindowsSmp                        = sprintSamples' + ((-(sprintWindow / 2 * phaseData.fsample) + 1):1:((sprintWindow / 2 * phaseData.fsample)));

    % find artifacts
    sprintBArtifact                         = any(artifactData.bArtifact(sprintWindowsSmp), 2)';

    % upsample slope, theta frequency, and artifact index
    sprintTimeWindow                        = mode(diff(sprintTime));
    sprintFsample                           = 1 / sprintTimeWindow;
    upsamplingRatio                         = phaseData.fsample / sprintFsample;
    slope                                   = repelem(slopeOrig, upsamplingRatio);
    sprintArtifactIdx                       = repelem(sprintBArtifact, upsamplingRatio);

    % add NaNs to SPRiNT data to start and end borders for padding
    missingSprint                           = size(artifactData.bArtifact, 2) - size(slope, 2);
    slope                                   = cat(2, repelem(NaN, missingSprint / 2), slope, repelem(NaN, missingSprint / 2));
    sprintArtifactIdx                       = cat(2, repelem(1, missingSprint / 2), sprintArtifactIdx, repelem(1, missingSprint / 2));

    % median split slope index
    slopeIdx                                = double(slope > median(slopeOrig(~sprintBArtifact))); 

    % remove artifacts
    slope(sprintArtifactIdx == 1)           = NaN;
    slopeIdx(sprintArtifactIdx == 1)        = NaN;

    %% cycle-by-cycle analysis

    % cycle start and end times
    cycleStart                              = bycycleData.bycycleTable.sample_last_peak;
    cycleEnd                                = bycycleData.bycycleTable.sample_next_peak;

    % cycle index
    cycleIdx                                = zeros(phaseData.sampleinfo);
    cycleIdx([cycleStart; cycleEnd(end)])   = 1;
    cycleIdx                                = cumsum(cycleIdx);
    cycleIdx(cycleEnd(end):end)             = NaN;

    % find cycles during artifacts
    cycleArtifactIdx                        = unique(cycleIdx(artifactData.bArtifact))';
    cycleArtifactIdx                        = cycleArtifactIdx(~isnan(cycleArtifactIdx) & cycleArtifactIdx > 0);

    % remove from cycle start and end arrays
    cycleStart(cycleArtifactIdx)            = NaN;
    cycleEnd(cycleArtifactIdx)              = NaN;

    % cycle frequency
    cycleFreq           = phaseData.fsample ./ (cycleEnd - cycleStart);

    % cycle frequency for each sample
    dataCycleFreq       = nan(phaseData.sampleinfo);
    dataCycleFreq(~isnan(cycleIdx) & cycleIdx > 0) = cycleFreq(cycleIdx((~isnan(cycleIdx) & cycleIdx > 0)));

    %% percentages oscillations and slope

    % burst index
    oscIdx                                  = double(bycycleData.bBurst);

    % oscillation and high slope
    oscHighSlopePerc(iWire, 1)              = sum(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx & slopeIdx == 1) == 1) / ...
        size(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx), 2);
    oscLowSlopePerc(iWire, 1)               = sum(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx & slopeIdx == 0) == 1) / ...
        size(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx), 2);
    noOscHighSlopePerc(iWire, 1)            = sum(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx & slopeIdx == 1) == 0) / ...
        size(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx), 2);
    noOscLowSlopePerc(iWire, 1)             = sum(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx & slopeIdx == 0) == 0) / ...
        size(oscIdx(~artifactData.bArtifact & ~sprintArtifactIdx), 2);

    %% extract further data

    % remove artifacts from phase data
    thisPhaseData(bArtifact)                = NaN;

    % assign data
    thisComplex                             = single(thisPhaseData);
    allReal{iWire, 1}                       = real(thisComplex);
    allImag{iWire, 1}                       = imag(thisComplex);
    allAngle{iWire, 1}                      = angle(thisComplex);

    % loop through distractor tasks
    allDistrFreqs = nan(size(distrTask, 1), 1);
    for iDistr = 1:size(distrTask, 1)
        allDistrFreqs(iDistr, 1)            = median(dataCycleFreq(distrTask(iDistr, 1):distrTask(iDistr, 2)), 'omitnan');
    end

    % mean frequency during successful and unsuccessful distractor tasks
    succDistrFreqs(iWire, 1)                = median(allDistrFreqs(distrTaskPerf), 'omitnan');
    failDistrFreqs(iWire, 1)                = median(allDistrFreqs(~distrTaskPerf), 'omitnan');

    % angular correction through centering
    realMean                                = movmean(real(thisComplex), params.windowSize * phaseData.fsample);
    imagMean                                = movmean(imag(thisComplex), params.windowSize * phaseData.fsample);
    recenteredData                          = thisComplex - (realMean + 1i * imagMean);
    allAngleCorr{iWire, 1}                  = angdiff(angle(recenteredData), angle(thisComplex));

    % correlations between generalized phase and narrow low and high theta
    lowThetaData                            = load(fullfile(paths.phaseData, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.lowThetaDataName));
    highThetaData                           = load(fullfile(paths.phaseData, subName, sessName, chanDir(wireIdx(iWire, 3)).name, params.highThetaDataName));
    bLowThetaCycles                         = dataCycleFreq > params.lowThetaBand(1, 1) & dataCycleFreq < params.lowThetaBand(1, 2);
    bHighThetaCycles                        = dataCycleFreq > params.highThetaBand(1, 1) & dataCycleFreq < params.highThetaBand(1, 2);
    lowThetaDiff{iWire, 1}                  = single(angdiff(angle(thisPhaseData(bLowThetaCycles)), angle(lowThetaData.trial{1, 1}(bLowThetaCycles))));
    highThetaDiff{iWire, 1}                 = single(angdiff(angle(thisPhaseData(bHighThetaCycles)), angle(highThetaData.trial{1, 1}(bHighThetaCycles))));

    %% correlation between slope and frequency

    % correlation
    [corrSlope(iWire, 1), ~]                = corr(slope(~isnan(slope))', dataCycleFreq(~isnan(slope))', 'Type', 'Spearman');

    % linear fit steepness
    polyFitSlopeCycle(iWire, :)             = polyfit(slope(~isnan(slope))', dataCycleFreq(~isnan(slope))', 1);

    %% collect further results

    % mean of analytic signal
    meanAnalytic(iWire, 1)                  = mean(thisPhaseData, 'omitnan');

    % standard deviation 
    stdAnalytic(iWire, 1)                   = std(thisPhaseData, 'omitnan');

    % mean instantaneous frequency
    meanInstFreq(iWire, 1)                  = mean(phaseData.frequency{1, 1});

    % % mean resultant vector length
    [mrvl(iWire, 1), meanAngle(iWire, 1)]   = circ_axialmean(angle(thisPhaseData), [], 2);
end

%% transform data
allAngle        = cat(2, allAngle{:});
allReal         = cat(2, allReal{:});
allImag         = cat(2, allImag{:});
allAngleCorr    = cat(2, allAngleCorr{:});

%% remove NaNs
allAngle        = allAngle(~isnan(allAngle));
allImag         = allImag(~isnan(allReal));
allReal         = allReal(~isnan(allReal));
allAngleCorr    = allAngleCorr(~isnan(allAngleCorr));

%% mean + sem

% mean
meanReal        = mean(allReal);
meanImag        = mean(allImag);

% sem
semReal         = std(allReal) / sqrt(size(allReal, 2));
semImag         = std(allImag) / sqrt(size(allImag, 2));

%% collect power spectra
powSpctrmFreqs  = FreqanalysisData{1, 1}.freq;
allPowspctra    = cell2mat(cellfun(@(x) x.powspctrm, FreqanalysisData, 'Uni', 0));

% plot mean powerspectrum
powerSpctrmFig  = figure;
TG_ShadeSEM_20210714(powSpctrmFreqs, allPowspctra, 'k', 0.5);
xlabel('Frequency (Hz)');
ylabel('Power');
xticks(0:10);
xlim([0, 11]);
ylim([0, 0.5]);
set(gca, 'tickDir', 'out', 'box', 'off');
print(powerSpctrmFig, fullfile(paths.save, '1_10_Hz_PowerSpectrum'), '-dsvg', '-r300');

%% get session-wise results

% session index
[sessNum, ~, sessIdx]   = unique(wireIdx(:, 1:2), 'row');

% loop through sessions
thissSessMeanReal       = nan(size(sessNum, 1), 1);
thissSessMeanImag       = nan(size(sessNum, 1), 1);
thisSessSuccDistrFreqs  = nan(size(sessNum, 1), 1);
thisSessFailDistrFreqs  = nan(size(sessNum, 1), 1);
thisSessCorrSlope       = nan(size(sessNum, 1), 1);
thisSessFitSlopeCycle   = nan(size(sessNum, 1), 2);
thisSessOscHighSlope    = nan(size(sessNum, 1), 1);
thisSessOscLowSlope     = nan(size(sessNum, 1), 1);
thisSessNoOscHighSlope  = nan(size(sessNum, 1), 1);
thisSessNoOscLowSlope   = nan(size(sessNum, 1), 1);
for iSess = 1:size(sessNum, 1)

    % get mean values for this session
    thissSessMeanReal(iSess, 1)         = mean(real(meanAnalytic(sessIdx == iSess)));
    thissSessMeanImag(iSess, 1)         = mean(imag(meanAnalytic(sessIdx == iSess)));
    thisSessSuccDistrFreqs(iSess, 1)    = median(succDistrFreqs(sessIdx == iSess), 'omitnan');
    thisSessFailDistrFreqs(iSess, 1)    = median(failDistrFreqs(sessIdx == iSess), 'omitnan');
    thisSessCorrSlope(iSess, 1)         = mean(corrSlope(sessIdx == iSess), 'omitnan');
    thisSessFitSlopeCycle(iSess, :)     = mean(polyFitSlopeCycle(sessIdx == iSess, :), 'omitnan');
    thisSessOscHighSlope(iSess, 1)      = mean(oscHighSlopePerc(sessIdx == iSess));
    thisSessOscLowSlope(iSess, 1)       = mean(oscLowSlopePerc(sessIdx == iSess));
    thisSessNoOscHighSlope(iSess, 1)    = mean(noOscHighSlopePerc(sessIdx == iSess));
    thisSessNoOscLowSlope(iSess, 1)     = mean(noOscLowSlopePerc(sessIdx == iSess));    
end

%% oscillation and slope 2x2 table

% means
oscillationHighSlope    = mean(thisSessOscHighSlope);
oscillationLowSlope     = mean(thisSessOscLowSlope);
noOscillationHighSlope  = mean(thisSessNoOscHighSlope);
noOscillationLowSlope   = mean(thisSessNoOscLowSlope);

% statistical test
[~, pSlope, ~, statsSlope]  = ttest(thisSessOscHighSlope, thisSessOscLowSlope);

% standard errors
semOscHighSlope         = std(thisSessOscHighSlope) ./ sqrt(size(thisSessOscHighSlope, 1));
semOscLowSlope          = std(thisSessOscLowSlope) ./ sqrt(size(thisSessOscLowSlope, 1));
semNoOscHighSlope       = std(thisSessNoOscHighSlope) ./ sqrt(size(thisSessNoOscHighSlope, 1));
semNoOscLowSlope        = std(thisSessNoOscLowSlope) ./ sqrt(size(thisSessNoOscLowSlope, 1));

% bar plot
oscSlopeFig     = figure;
b1              = bar([1, 2, 3, 4], [oscillationHighSlope, oscillationLowSlope, noOscillationHighSlope, noOscillationLowSlope], 'FaceColor', [0.5, 0.5, 0.5]);
hold on;
e1              = errorbar([1, 2, 3, 4], [oscillationHighSlope, oscillationLowSlope, noOscillationHighSlope, noOscillationLowSlope], ...
    [semOscHighSlope,  semOscLowSlope, semNoOscHighSlope, semNoOscLowSlope]);
e1.Color        = [0 0 0];
e1.LineStyle    = 'none';
ylim([0, 0.4]);
set(gca, 'tickDir', 'out', 'box', 'off');
print(oscSlopeFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_PercentagesOscillationsAndSlopes')), '-dsvg', '-r300');

%% statistical tests

% frequencies of successful and unsuccessful distractor task
[~, pvalDistr, ~, statsDistr]           = ttest(thisSessSuccDistrFreqs, thisSessFailDistrFreqs);

% correlation between slope and frequency
[~, pvalSlopeFreq, ~, statsSlopeFreq]   = ttest(thisSessCorrSlope);

%% plot cycle frequencies of successful vs. unsuccessful distractor task

% figure
distrFig        = figure;
bD              = bar([1, 2], [mean(thisSessSuccDistrFreqs, 'omitnan'), mean(thisSessFailDistrFreqs, 'omitnan')], 'FaceColor', [0.5, 0.5, 0.5]);
hold on;
semSuccDistr    = std(thisSessSuccDistrFreqs, 'omitnan') / sqrt(size(thisSessSuccDistrFreqs(~isnan(thisSessSuccDistrFreqs)), 1));
semFailDistr    = std(thisSessFailDistrFreqs, 'omitnan') / sqrt(size(thisSessFailDistrFreqs(~isnan(thisSessFailDistrFreqs)), 1));
e1              = errorbar([1, 2], [mean(thisSessSuccDistrFreqs, 'omitnan'), mean(thisSessFailDistrFreqs, 'omitnan')], ...
    [semSuccDistr, semFailDistr]);
e1.Color        = [0 0 0];
e1.LineStyle    = 'none';
set(gca, 'tickDir', 'out', 'box', 'off');
print(distrFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_OscillationsDistractorTask')), '-dsvg', '-r300');

%% plot cycle frequency correlation with slope
rhoSlopeCycleFig    = figure;
h1                  = histogram(thisSessCorrSlope, -0.3:0.01:-0.1, 'FaceColor', [0.5, 0.5, 0.5]);
set(gca, 'tickDir', 'out', 'box', 'off');
xticks([-0.3, -0.2, -0.1]);
print(rhoSlopeCycleFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_RhoSlopeCycleFrequency')), '-dsvg', '-r300');

%% plot steepness of linear fits
fitSlopeCycleFig    = figure;
h2                  = histogram(thisSessFitSlopeCycle(:, 1), -2:0.1:0, 'FaceColor', [0.5, 0.5, 0.5]);
set(gca, 'tickDir', 'out', 'box', 'off');
xticks([-2, -1, 0]);
print(fitSlopeCycleFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_FitSteepnessSlopeCycle')), '-dsvg', '-r300');

%% plot angles
[grandMrvl, grandMean] = circ_axialmean(allAngle, [], 2);
anglesFig       = figure;
polarhistogram(allAngle, 60, 'normalization', 'probability', 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1);
set(gca, 'ThetaTickLabel', params.polarHistLabels);
title(strcat('MRVL = ', num2str(grandMrvl)));
rlim([0, 0.02]);
print(anglesFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_allPhaseAngles')), '-dsvg', '-r300');

%% plot angle correction through centering
[corrMrvl, corrMean] = circ_axialmean(allAngleCorr, [], 2);
anglesCorrFig       = figure;
polarhistogram(allAngleCorr, 360, 'normalization', 'probability', 'FaceColor', params.phaseColor(1:3), 'FaceAlpha', params.phaseColor(4));
set(gca, 'ThetaTickLabel', params.polarHistLabels);
title(strcat('MRVL = ', num2str(corrMrvl)));
rlim([0, 0.4]);
print(anglesCorrFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_allPhaseAngleCorrectionsThroughCentering')), '-dsvg', '-r300');

%% complex plane samples

% heatmap
heatMapCounts       = histcounts2(allReal(:, 1:4:end), allImag(:, 1:4:end), params.heatmapEdges, params.heatmapEdges); % downsampled to 500 Hz
heatMapCounts       = heatMapCounts * 4; % upsample to 2000 Hz again (this procedure is an approximation to not run out of memory)
logHeatMapCounts    = log10(heatMapCounts);
logHeatMapCounts(isinf(logHeatMapCounts)) = NaN;

% plot the 2D histogram
complexFig  = figure();
pcolor(params.heatmapCenters, params.heatmapCenters, logHeatMapCounts);
shading interp;
axis equal;
box off;
set(gca, 'TickDir', 'out');
xlim([-500, 500]);
ylim([-500, 500]);
xline(0);
yline(0);

% colorbar
c1              = colorbar;

% print figure
print(complexFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_complexPlane_500_downsampled500Hz')), '-djpeg', '-r600');

%% plot differences between generalized phase and low theta
lowThetaDiffFig = figure();
lowThetaCircStd = rad2deg(circ_std(cat(2, lowThetaDiff{:}), [], [], 2));
lowThetaMRVL    = circ_axialmean(cat(2, lowThetaDiff{:}), [], 2);
polarhistogram(cat(2, lowThetaDiff{:}), 60, 'normalization', 'probability', 'FaceColor', params.phaseColor(1:3), 'FaceAlpha', params.phaseColor(4));
set(gca, 'ThetaTickLabel', params.polarHistLabels);
rlim([0, 0.1]);
title(['low theta; MRVL = ', num2str(lowThetaMRVL)]);
print(lowThetaDiffFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_lowThetaVsGeneralizedPhase')), '-dsvg', '-r300');

%% plot differences between generalized phase and high theta
highThetaDiffFig = figure();
highThetaCircStd = rad2deg(circ_std(cat(2, highThetaDiff{:}), [], [], 2));
highThetaMRVL    = circ_axialmean(cat(2, highThetaDiff{:}), [], 2);
polarhistogram(cat(2, highThetaDiff{:}), 60, 'normalization', 'probability', 'FaceColor', params.phaseColor(1:3), 'FaceAlpha', params.phaseColor(4));
set(gca, 'ThetaTickLabel', params.polarHistLabels);
rlim([0, 0.1]);
title(['high theta; MRVL = ', num2str(highThetaMRVL)]);
print(highThetaDiffFig, fullfile(paths.save, strcat(params.filterBandStr, '_Hz_highThetaVsGeneralizedPhase')), '-dsvg', '-r300');