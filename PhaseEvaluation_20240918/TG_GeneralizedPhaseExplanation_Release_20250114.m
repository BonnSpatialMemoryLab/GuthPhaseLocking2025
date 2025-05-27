%==========================================================================
% This script demonstrates an example of the generalized phase procedure, 
% illustrating its key steps.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% paths
paths               = struct();
paths.LFPData       = 'D:\TreasureHunt\MicroDownsampled_20210910\TH11\session_0\chan3\datacutFt2000Hz.mat';
paths.filteredData  = 'D:\TreasureHunt\MicroFiltered_20210930\TH11\session_0\chan3\datacutFt2000HzBP_1_10.mat';
paths.phaseData     = 'D:\TreasureHunt\MicroPhases_20210930\TH11\session_0\chan3\datacutFt2000HzBP_generalized_1_10.mat';
paths.hilbertData   = 'D:\TreasureHunt\MicroPhases_20210930\TH11\session_0\chan3\datacutFt2000HzBP_hilbert_1_10.mat';
paths.save          = 'D:\TreasureHunt\PhaseEvaluation_20240918';

% parameters
params              = struct();
params.samplesToCut = 6585000:6588000;
params.interpolated = 1147:1357;

% load data
lfpData             = load(paths.LFPData);
filteredData        = load(paths.filteredData);
phaseData           = load(paths.phaseData);
hilbertData         = load(paths.hilbertData);

% cut segment
thisLFP             = lfpData.trial{1, 1}(params.samplesToCut);
thisFiltered        = filteredData.trial{1, 1}(params.samplesToCut);
thisPhase           = phaseData.trial{1, 1}(params.samplesToCut);
thisHilbert         = hilbertData.trial{1, 1}(params.samplesToCut);
thisTime            = (0:1:size(thisLFP, 2) - 1) / lfpData.fsample;

% plot raw signal
rawSignalFig        = figure;
plot(thisTime, thisLFP, 'Color', 'k', 'LineWidth', 2);
hold on;
plot(thisTime(params.interpolated), thisLFP(params.interpolated), 'Color', 'b', 'LineWidth', 2);
set(gca, 'tickDir', 'out');
ylim([-400, 400]);
yticks([-400, 0, 400]);
box off;
xlabel('Time (s)');
ylabel('Amplitude (µV)');
set(gcf, 'Renderer', 'painter');
saveas(rawSignalFig, fullfile(paths.save, 'ExampleRawSignal.svg'));

% plot filtered signal
filteredSignalFig   = figure;
plot(thisTime, thisFiltered, 'Color', 'k', 'LineWidth', 2);
hold on;
plot(thisTime(params.interpolated), thisFiltered(params.interpolated), 'Color', 'b', 'LineWidth', 2);
set(gca, 'tickDir', 'out');
ylim([-400, 400]);
yticks([-400, 0, 400]);
box off;
xlabel('Time (s)');
ylabel('Amplitude (µV)');
set(gcf, 'Renderer', 'painter');
saveas(filteredSignalFig, fullfile(paths.save, 'ExampleFilteredSignal.svg'));

% plot hilbert signal
hilbertFig          = figure;
hilbertPlotTime     = thisTime;
hilbertPlotReal     = real(thisHilbert);
hilbertPlotImag     = imag(thisHilbert);
hilbertPlotTime(params.interpolated)    = NaN;
hilbertPlotReal(params.interpolated)    = NaN;
hilbertPlotImag(params.interpolated)    = NaN;
plot3(hilbertPlotTime, hilbertPlotReal, hilbertPlotImag, 'Color', 'k', 'LineWidth', 2);
hold on;
plot3(thisTime(params.interpolated), real(thisHilbert(params.interpolated)), imag(thisHilbert(params.interpolated)), 'Color', 'b', 'LineWidth', 2);
set(gca, 'tickDir', 'out');
xlabel('Time (s)');
ylabel('Real part (µV)');
zlabel('Imaginary part (µV)');
ylim([-400, 400]);
yticks([-400, 0, 400]);
zlim([-400, 400]);
zticks([-400, 0, 400]);
view(45, 45);
set(gcf, 'Renderer', 'painter');
saveas(hilbertFig, fullfile(paths.save, 'ExampleAnalyticSignal.svg'));

% change perspective to complex plane
view(90, 0);
saveas(hilbertFig, fullfile(paths.save, 'ExampleAnalyticSignalComplexPlane.svg'));

% plot in polar histogram
polarPlotFig        = figure;
polarplot(angle(thisHilbert), abs(thisHilbert), 'Color', 'k', 'LineWidth', 2);
hold on;
polarplot(angle(thisHilbert(params.interpolated)), abs(thisHilbert(params.interpolated)), 'Color', 'b', 'LineWidth', 2);
set(gcf, 'Renderer', 'painter');
saveas(polarPlotFig, fullfile(paths.save, 'ExampleAnalyticSignalPhaseEstimates.svg'));

% plot hilbert phase estimate
hilbertPhaseFig   = figure;
plot(thisTime, thisFiltered, 'Color', 'k', 'LineWidth', 2);
hold on;
plot(thisTime(params.interpolated), thisFiltered(params.interpolated), 'Color', 'b', 'LineWidth', 2);
ylim([-400, 400]);
yticks([-400, 0, 400]);
yyaxis right;
plot(thisTime, rad2deg(angle(thisHilbert)), 'LineWidth', 2);
hold on;
plot(thisTime(params.interpolated), rad2deg(angle(thisHilbert(params.interpolated))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-');
set(gca, 'tickDir', 'out');
ylim([-350, 350]);
yticks([-180, 0, 180]);
box off;
xlabel('Time (s)');
ylabel('Angle (degree)');
set(gcf, 'Renderer', 'painter');
saveas(hilbertPhaseFig, fullfile(paths.save, 'ExampleHilbertAngle.svg'));

% plot generalized phase estimate
generalizedPhaseFig   = figure;
plot(thisTime, thisFiltered, 'Color', 'k', 'LineWidth', 2);
hold on;
plot(thisTime(params.interpolated), thisFiltered(params.interpolated), 'Color', 'b', 'LineWidth', 2);
ylim([-400, 400]);
yticks([-400, 0, 400]);
yyaxis right;
plot(thisTime, rad2deg(angle(thisPhase)), 'LineWidth', 2);
hold on;
plot(thisTime(params.interpolated), rad2deg(angle(thisPhase(params.interpolated))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-');
set(gca, 'tickDir', 'out');
ylim([-350, 350]);
yticks([-180, 0, 180]);
box off;
xlabel('Time (s)');
ylabel('Angle (degree)');
set(gcf, 'Renderer', 'painter');
saveas(generalizedPhaseFig, fullfile(paths.save, 'ExampleGeneralizedAngle.svg'));
