%==========================================================================
% This script creates example figures for cells with a firing rate increase
% during the encoding period (during object presentation).
%
% Tim Guth, 2025
%==========================================================================

% start
clear; close all; clc;

% add paths
paths                   = [];
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results from phase analysis
paths.behlog            = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.SUData            = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.save              = 'D:\TreasureHunt\ObjectCells_20231208';

% add functions
addpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions');
addpath(genpath('D:\External\Functions'));

%% settings

% parameters
param                   = [];
param.segBlDur          = 1500; % baseline duration in ms before segment
param.segEncDur         = 1500; % encoding duration in ms
param.buffer            = 2000; % time buffer in ms before and after segment
param.binWidth          = 50; % time bins for firing rate in ms
param.blWindow          = [-1500, 0]; % baseline window in ms
param.timeVector        = -param.segBlDur:param.binWidth:param.segEncDur + param.buffer;

% baseline index
[~, blIdxStart]         = min(abs(param.timeVector - param.blWindow(1)));
[~, blIdxEnd]           = min(abs(param.timeVector - param.blWindow(2)));
param.blIndices         = blIdxStart:blIdxEnd; % indices for baseline

%% example cell indices
exampleIdx = [6, 1, 56, 1; ...
    12, 0, 12, 2; ...
    14, 0, 33, 1; ...
    14, 0, 34, 1; ...
    17, 0, 4, 1; ...
    17, 0, 20, 4];

%% load subject information
phaseSettings = load(fullfile(paths.phaseRes, 'settings.mat'));
phaseRes      = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));
brainReg      = {phaseRes.phaseRes.brainRegion}';
allIdx        = cat(1, phaseRes.phaseRes.idx);
bObjCell      = cat(1, phaseRes.phaseRes.bObjectCell);

% loop through cells
for iCell = 1:size(exampleIdx, 1)
    
    %% load data

    % cell name
    unitName        = regexprep(num2str(exampleIdx(iCell, :)), '\s+', '_');

    % index in phase results
    [~, resIdx]     = ismember(exampleIdx(iCell, :), allIdx, 'rows');

    % subject ID
    subNum          = exampleIdx(iCell, 1);
    subjectID       = phaseSettings.subjects{subNum};

    % session name
    sessName        = strcat('session_', num2str(exampleIdx(iCell, 2)));

    % encoding segment info
    segmentInfo     = load(fullfile(paths.behlog, subjectID, sessName, 'segmentInfo_20230208.mat'));
    segments        = segmentInfo.encoding  / (segmentInfo.fsample / 1000); % convert to ms

    % get available single unit data
    mwChanDir       = TG_GetChanDir_20210812(paths.SUData, subjectID, sessName);

    % path for this channel's unit data
    mwChanPath      = fullfile(paths.SUData, subjectID, sessName, mwChanDir(exampleIdx(iCell, 3)).name);

    % load cluster
    t               = load(fullfile(mwChanPath, 'times_datacut.mat'));

    % get cluster
    thisClus        = t.cluster_class(t.cluster_class(:, 1) == exampleIdx(iCell, 4), 2); % cluster-number, time in ms
    
    %% get segment-wise spikes and create rasterplot

    % create figure
    rasFig          = figure('units', 'centimeters', 'position', [10, 10, 12, 12]);
    subplot(2, 1, 2);
    hold on;
    
    % loop through segments
    segSpikeTime            = cell(size(segments, 1), 1);
    for iSeg = 1:size(segments, 1)

        %% get segment information

        % get segment information
        segSpikeIdx             = find(thisClus >= (segments(iSeg, 1) - param.buffer) & thisClus <= (segments(iSeg, 2) + param.buffer));
        segNumSpikes            = size(segSpikeIdx, 1);
        segSpikeTimeAbs         = thisClus(segSpikeIdx);
        segSpikeTime{iSeg, 1}   = (segSpikeTimeAbs - segments(iSeg, 1));

        % plot horizontal lines for each element
        yVal                    = iSeg * ones(size(segSpikeTime{iSeg, 1}));
        for iSpike = 1:size(segSpikeTime{iSeg, 1}, 1)
            line([segSpikeTime{iSeg, 1}(iSpike), segSpikeTime{iSeg, 1}(iSpike)], [yVal(iSpike) - 0.4, yVal(iSpike) + 0.4], 'Color', 'k', 'LineWidth', 1);
        end
    end

    % adjust x-axis
    xlim([-1500, 3000]);
    xlabel('Time (ms)');
    xline([0, 1500], 'LineWidth', 1);
    ylim([0, size(segments, 1) + 1]);
    set(gca, 'YDir', 'reverse', 'tickDir', 'out');
    yticks([1, 50, 100]);
    ylabel('Segment number');
    set(gca, 'FontSize', 14);

    %% add mean firing rate

    % initialize firing rate vector
    spikeHistAll = zeros(size(segSpikeTime, 1), size(param.timeVector, 2) - 1);

    % compute firing rate for each segment
    for iSeg = 1:size(segSpikeTime, 1)

        % spike times
        spikeTimes              = segSpikeTime{iSeg};

        % compute histogram for the current cell
        spikeHist               = histcounts(spikeTimes, param.timeVector);

        % subtract baseline firing rate
        blRate                  = mean(spikeHist(param.blIndices));
        spikeHist               = spikeHist - blRate;

        % accumulate spike histogram
        spikeHistAll(iSeg, :)   = spikeHist;
    end

    % compute mean firing rate and SEM
    meanFiringRate      = mean(spikeHistAll, 1);
    semFiringRate       = std(spikeHistAll, 1) / sqrt(size(spikeHistAll, 1));

    % time axis
    timeAxis            = (param.timeVector(1 : end - 1) + param.timeVector(2 : end)) / 2;
    
    % plot firing rate with shaded area for SEM
    subplot(2, 1, 1);
    hold on;
    [lineOut, fillOut]  = TG_ShadeSEM_20210714(timeAxis, spikeHistAll, [0.4, 0.4, 0.4], 0.6);
    
    % firing rate
    xlim([-1500, 3000]);
    xline([0, 1500], 'LineWidth', 1);
    ylabel('Firing Rate (Hz)');
    ylim([-0.5, 1]);
    set(gca, 'tickDir', 'out');
    title(strsplit(brainReg{resIdx}, '_'), 'Interpreter', 'none');
    set(gca, 'FontSize', 14);

    % print
    set(rasFig, 'renderer', 'painters');
    saveas(rasFig, fullfile(paths.save, [unitName, '_Rasterplot.svg']), 'svg');

    %% spike density plot

    % get spike data
    thisClusSpike   = t.spikes(t.cluster_class(:, 1) == exampleIdx(iCell, 4), :);
    spkFig          = figure('units', 'centimeters', 'position', [10, 10, 3, 2.5]);
    hold on;

    % spike density plot
    spikeTime       = ((0:size(thisClusSpike, 2) - 1) ./ t.par.sr) .* 1000; % [msec]
    outPlot         = TG_DensityPlot(spikeTime, thisClusSpike);
    hold off;
    set(gca, ...
        'xlim', round([min(spikeTime), max(spikeTime)]), ...
        'xtick', round([min(spikeTime), max(spikeTime)]), ...
        'ylim', [outPlot.lbound, outPlot.ubound], ...
        'ytick', [outPlot.lbound, outPlot.ubound], ...
        'ticklength', [0, 0]);
    box off;
    colormap(spkFig, 'parula');
    % title('Spike density plot');

    % enhance axes
    xlabel('ms', ...
        'FontUnits', 'centimeters', 'units', 'normalized', 'position', [0.5, 0, 0]);
    ylabel('\muV', ...
        'FontUnits', 'centimeters', 'units', 'normalized', 'position', [-0.05, 0.5, 0]);

    % print
    saveas(spkFig, fullfile(paths.save, [unitName, '_SpikeDensityPlot']), 'png');
end
  
