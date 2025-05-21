%==========================================================================
% This script tests the generalized phase approach on different simulated
% signals.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% random seed
rng(1);

% paths
paths                   = struct();
paths.save              = 'D:\TreasureHunt\PhaseEvaluation_20240918';

% parameters
params                  = struct();
params.signalName       = 'FrequencyCouplingLowDominantBelow1Hz';
params.sr               = 2000; % sampling rate in Hz
params.frequency        = 2; % Hz
params.freqRange        = [1, 10]; % frequency band
params.timeStart        = 0;
params.timeEnd          = 10;
params.duration         = 3600;
params.toiStart         = 2.5;
params.toiEnd           = 5;
params.binEdges         = 0:0.5:10;
params.binCenters       = (params.binEdges(1:end - 1) + params.binEdges(2:end)) / 2;

% set configurations for bandpass filter
cfg                     = [];
cfg.bpfilter            = 'yes';
cfg.bpfilttype          = 'fir';
cfg.bpfreq              = params.freqRange;
cfg.demean              = 'yes';

% time axis
thisTimeOrig            = params.timeStart:(1/params.sr):params.duration;

% signal
thisSignalOrig          = sin(2 * pi * thisTimeOrig * params.frequency);

% oscillation with amplitude modulation
if strcmp(params.signalName, 'AmplitudeModulation')
    thisSignalOrig              = thisSignalOrig .* (0.9 .* sin(2 * pi * thisTimeOrig * params.frequency / 20) + 0.1);
elseif strcmp(params.signalName, 'Tilted')
    thisSignalOrig              = sin(2 * pi * thisTimeOrig * params.frequency + (thisSignalOrig / 2));
elseif strcmp(params.signalName, 'Sawtooth')
    thisSignalOrig              = sawtooth(2 * pi * (thisTimeOrig + (1  / params.frequency / 4)) * params.frequency, 0.5);
elseif strcmp(params.signalName, 'Skewed')
    newTime                     = 1.5 .^ thisTimeOrig;
    newToiStart                 = 2.1335;
    newToiEnd                   = 4.7545;
    newSignal                   = sin(2 * pi * newTime * params.frequency);
    newSignalCut                = newSignal(thisTimeOrig >= newToiStart & thisTimeOrig <= newToiEnd);
    newSignalFlip               = cat(2, newSignalCut, flip(newSignalCut));
    newSignalRep                = repmat(newSignalFlip, [1, ceil((params.duration * params.sr) / size(newSignalFlip, 2)) + 1]);
    thisSignalOrig              = newSignalRep(5740:(params.duration * params.sr) + 5740);
elseif strcmp(params.signalName, 'Artifact')
    numSpikes                   = 7200;
    spikeAmplitude              = (2 * randi([0, 1], numSpikes, 1) - 1) * 5;
    spikeSamples                = 0.05 * params.sr;
    spikeSignal                 = gausswin(spikeSamples)' .* spikeAmplitude;
    spikeTimes                  = randi([1, (params.duration * params.sr)], numSpikes, 1);
    spikeIdx                    = spikeTimes + (1:spikeSamples);
    thisSignalOrig(spikeIdx)    = thisSignalOrig(spikeIdx) + spikeSignal;
elseif strcmp(params.signalName, 'FrequencyCouplingLowDominant')
    thisSignalOrig              = thisSignalOrig + (sin(2 * pi * thisTimeOrig * 9) .* 0.5);
elseif strcmp(params.signalName, 'FrequencyCouplingLowDominantBelow1Hz')
    thisSignalOrig              = thisSignalOrig + (sin(2 * pi * thisTimeOrig * 0.5) .* 2);
end

% transform into Fieldtrip structure
ftSignalOrig                    = struct();
ftSignalOrig.fsample            = params.sr;
ftSignalOrig.sampleinfo         = [1, params.duration * params.sr + 1];
ftSignalOrig.label              = {params.signalName};
ftSignalOrig.trial              = {thisSignalOrig};
ftSignalOrig.time               = {thisTimeOrig};

% bandpass filter data
ftFilteredOrig                  = ft_preprocessing(cfg, ftSignalOrig);
thisFilteredOrig                = ftFilteredOrig.trial{1, 1};

% generalized phase
thisComplexOrig     = TG_generalized_phase_vector(thisFilteredOrig', params.sr, 1);

% cut signal to time window of interest
bCut                = thisTimeOrig >= params.toiStart & thisTimeOrig <= params.toiEnd;
thisTime            = (0:1:sum(bCut) - 1) / params.sr;
thisSignal          = thisSignalOrig(bCut);
thisFiltered        = thisFilteredOrig(bCut);
thisComplex         = thisComplexOrig(bCut);

% compute means
thisMean            = mean(thisSignal);
thisMeanComplex     = mean(thisComplex);

% Fourier transform
thisSignalLength    = length(thisSignalOrig);
thisFourier         = fft(thisSignalOrig);
powerSpectrum       = (abs(thisFourier) .^ 2) / thisSignalLength;
thisFreqs           = (0:(thisSignalLength - 1)) * (params.sr / thisSignalLength);

% loop through bins
binnedPower         = zeros(1, length(params.binEdges) - 1);
for iBin = 1:length(params.binCenters)
    binIdx              = thisFreqs >= params.binEdges(iBin) & thisFreqs < params.binEdges(iBin + 1);
    binnedPower(iBin)   = mean(powerSpectrum(binIdx));
end

% plot power spectrum
powerFig            = figure('Units', 'centimeters', 'Position', [14, 14, 3, 3]);
plot(params.binCenters, binnedPower, 'Color', 'k', 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Power');
yticks([]);
xlim([0, 10]);
xticks([0, 5, 10]);
set(gca, 'tickDir', 'out');
box off;
set(gcf, 'Renderer', 'painters');
saveas(powerFig, fullfile(paths.save, strcat('5_SimulatedPowerSpectrum', params.signalName, '.svg')));

% plot real signal
realFig             = figure('Units', 'centimeters', 'Position', [14, 14, 5, 3]);
plot(thisTime, thisSignal, 'Color', 'k', 'LineWidth', 1);
hold on;
plot(thisTime, repmat(thisMean, size(thisSignal)), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'tickDir', 'out');
box off;
xticks([0, 1, 2]);
yticks([]);
xlabel('Time (s)');
ylabel('Real');
set(gcf, 'Renderer', 'painters');
saveas(realFig, fullfile(paths.save, strcat('1_SimulatedRealSignal', params.signalName, '.svg')));

% plot analytic signal
analyticFig                 = figure('Units', 'centimeters', 'Position', [14, 14, 5, 3]);
analyticTime                = thisTime;
analyticReal                = real(thisComplex);
analyticImag                = imag(thisComplex);
analyticMeanReal            = repmat(real(thisMeanComplex), size(thisSignal));
analyticMeanImag            = repmat(imag(thisMeanComplex), size(thisSignal));
plot3(analyticTime, analyticReal, analyticImag, 'Color', 'k', 'LineWidth', 1);
hold on;
plot3(thisTime, analyticMeanReal, analyticMeanImag, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'tickDir', 'out');
xticks([0, 1, 2]);
yticks([]);
zticks([]);
xlabel('Time (s)');
ylabel('Real');
zlabel('Imaginary');
view(45, 45);
set(gcf, 'Renderer', 'painters');
saveas(analyticFig, fullfile(paths.save, strcat('2_SimulatedAnalyticSignal', params.signalName, '.svg')));

% plot in polar histogram
polarFig    = figure('Units', 'centimeters', 'Position', [14, 14, 3, 3]);
polarplot(angle(thisComplex), abs(thisComplex), 'Color', 'k', 'LineWidth', 1);
hold on;
polarplot(angle(thisMeanComplex), abs(thisMeanComplex), 'x', 'MarkerSize', 10, 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5]);
set(gca, 'RTick', [], 'ThetaTickLabel', {'0'; ''; ''; '90'; ''; ''; '\pm180'; ''; ''; '-90'; ''; ''});
set(gcf, 'Renderer', 'painters');
saveas(polarFig, fullfile(paths.save, strcat('3_SimulatedPhases', params.signalName, '.svg')));

% phase estimate
phaseFig             = figure('Units', 'centimeters', 'Position', [14, 14, 6.5, 3]);
plot(thisTime, thisFiltered, 'Color', 'k', 'LineWidth', 1);
xticks([0, 1, 2]);
yticks([]);
xlabel('Time (s)');
ylabel('Real');
yyaxis right;
plot(thisTime, angle(thisComplex), 'LineWidth', 1);
ylabel('Phase');
yticks([-pi, pi]);
ylim([-pi, pi]);
yticklabels({'-180°', '180°'});
set(gca, 'tickDir', 'out');
box off;
set(gcf, 'Renderer', 'painters');
saveas(phaseFig, fullfile(paths.save, strcat('4_SimulatedSignalPhase', params.signalName, '.svg')));