%==========================================================================
% This script tests whether the MRVL, its rank in a distribution
% of surrogate MRVLs and a calculated z-value depend on the number of
% spikes
%
% Tim Guth, 2025
%==========================================================================

% start
clear; close all; clc;

% add paths
paths.save      = 'D:\TreasureHunt\Simulations_20231120';

% add functions
addpath(genpath('D:\External\Functions\'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% settings
spikeNum        = [10; 100; 1000];
cellNum         = 10000;
nSur            = 1001;
rng(444, 'twister'); % set random seed

% preallocate
subZ            = nan(1, cellNum);
zval            = nan(size(spikeNum, 1), cellNum);
PPC             = nan(size(spikeNum, 1), cellNum);
mrvl            = nan(size(spikeNum, 1), cellNum);
rank            = nan(size(spikeNum, 1), cellNum);
zscored         = nan(size(spikeNum, 1), cellNum);

% loop through example cells
for iC = 1:cellNum

    % display progress
    disp(iC);

    % random phase-locking neuron
    randPhasesPL    = cat(1, vmrand(0, 2, [250, 1]), 2 * pi * rand(750, 1) - pi);

    % subsampling 1000 spikes to 100 spikes
    cellSubZ = nan(nSur, 1);
    for iSp = 1:nSur

        % draw a 100 spikes subsample from the 1000 spikes
        subPhases               = datasample(randPhasesPL, spikeNum(2), 'Replace', false);
        [~, cellSubZ(iSp, 1)]   = circ_rtest(subPhases);
    end

    % calculate mean of the subsamples
    subZ(1, iC) = mean(cellSubZ);

    % loop through sample sizes
    for iSn = 1:size(spikeNum, 1)

        % draw n random spikes of neuron
        smpPhases            = datasample(randPhasesPL, spikeNum(iSn), 'Replace', false);

        % pairwise phase consistency
        PPC(iSn, iC)         = TG_PPC_20241128(smpPhases);

        % Rayleigh z-value
        [~, zval(iSn, iC)]   = circ_rtest(smpPhases);

        % mean resultant vector length
        mrvl(iSn, iC)        = round(circ_axialmean(smpPhases), 5);

        % surrogates
        surDistr             = 2 * pi * rand(spikeNum(iSn), nSur) - pi;
        surMrvl              = circ_axialmean(surDistr);

        % rank of mrvl
        rank(iSn, iC)        = sum(mrvl(iSn, iC) > surMrvl) / nSur;

        % z-score mrvl
        zscored(iSn, iC)     = (mrvl(iSn, iC) - mean(surMrvl)) / std(surMrvl);
    end
end

%% random phase-locking neurons

% 10 spikes
f1 = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
p1 = polarhistogram(randPhasesPL(1:spikeNum(1)), 24, 'FaceColor', [0.2, 0.2, 1]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
rlim([0, 2]);
print(f1, fullfile(paths.save, '10Spikes'), '-dsvg', '-r300');

% 100 spikes
f2 = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
p2 = polarhistogram(randPhasesPL(1:spikeNum(2)), 24, 'FaceColor', [0.4, 0.4, 0.4]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
rlim([0, 20]);
print(f2, fullfile(paths.save, '100Spikes'), '-dsvg', '-r300');

% 1000 spikes
f3 = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
p3 = polarhistogram(randPhasesPL(1:spikeNum(3)), 24, 'FaceColor', [0.4, 0.8, 0.2]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
rlim([0, 100]);
print(f3, fullfile(paths.save, '1000Spikes'), '-dsvg', '-r300');

%% histograms for PPC
mPPC = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% distributions
histogram(PPC(1, :), 'FaceColor', [0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(PPC(2, :), 'FaceColor', [0.4, 0.4, 0.4], 'EdgeColor', 'none');
histogram(PPC(3, :), 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'none');

% means
xline(mean(PPC(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 2);
xline(mean(PPC(2, :)), 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
xline(mean(PPC(3, :)), 'Color', [0.4, 0.8, 0.2], 'LineWidth', 2);

% settings
ylim([0, 1500]);
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of simulated units');
xlabel('Pairwise phase consistency');
print(mPPC, fullfile(paths.save, 'PPC'), '-dsvg', '-r300');

%% histograms for mrvl
mF = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% mrvl distributions
histogram(mrvl(1, :), 'FaceColor',[0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(mrvl(2, :), 'FaceColor', [0.4, 0.4, 0.4], 'EdgeColor', 'none');
histogram(mrvl(3, :), 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'none');

% means
xline(mean(mrvl(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 2);
xline(mean(mrvl(2, :)), 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
xline(mean(mrvl(3, :)), 'Color', [0.4, 0.8, 0.2], 'LineWidth', 2);

% settings
ylim([0, 1500]);
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of simulated units');
xlabel('MRVL');
print(mF, fullfile(paths.save, 'MRVL'), '-dsvg', '-r300');

%% histograms for Rayleigh's z
mRz = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% distributions
histogram(zval(1, :), 'FaceColor', [0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(zval(2, :), 'FaceColor', [0.4, 0.4, 0.4], 'EdgeColor', 'none');
histogram(zval(3, :), 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'none');

% means
xline(mean(zval(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 2);
xline(mean(zval(2, :)), 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
xline(mean(zval(3, :)), 'Color', [0.4, 0.8, 0.2], 'LineWidth', 2);

% settings
ylim([0, 1500]);
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of simulated units');
xlabel('Rayleigh z-value');
print(mRz, fullfile(paths.save, 'Rayleigh'), '-dsvg', '-r300');

%% histograms for rank
mRank = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% distributions
histogram(rank(1, :), 'FaceColor', [0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(rank(2, :), 'FaceColor', [0.4, 0.4, 0.4], 'EdgeColor', 'none');
% histogram(rank(3, :), 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'none');

% means
xline(mean(rank(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 2);
xline(mean(rank(2, :)), 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
xline(mean(rank(3, :)), 'Color', [0.4, 0.8, 0.2], 'LineWidth', 2);

% settings
ylim([0, 1500]);
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of simulated units');
xlabel('Rank in 1001 surrogates');
print(mRank, fullfile(paths.save, 'Rank'), '-dsvg', '-r300');

%% histograms for zval
mZval = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% distributions
histogram(zscored(1, :), 'FaceColor', [0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(zscored(2, :), 'FaceColor', [0.4, 0.4, 0.4], 'EdgeColor', 'none');
histogram(zscored(3, :), 'FaceColor', [0.4, 0.8, 0.2], 'EdgeColor', 'none');

% means
xline(mean(zscored(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 2);
xline(mean(zscored(2, :)), 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
xline(mean(zscored(3, :)), 'Color', [0.4, 0.8, 0.2], 'LineWidth', 2);

% settings
ylim([0, 1500]);
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of simulated units');
xlabel('Z-score (against 1001 surrogates)');
print(mZval, fullfile(paths.save, 'Zscore'), '-dsvg', '-r300');

%% pairwise phase consistency 100 vs 1000 spikes

% mean and SEM of Rayleigh z-value
meanPPC          = mean(PPC, 2);
semPPC           = std(PPC, [], 2) / sqrt(size(PPC, 2));

% plot PPC for 100 and 1000 spikes
ppcCompFig       = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
b2               = bar([meanPPC(2), meanPPC(3)], 'FaceColor', 'flat', 'FaceAlpha', 0.6);
set(gca, 'XTickLabel', {'100 spikes'; '1000 spikes'});
b2.CData(1, :)   = [0.4, 0.4, 0.4];
b2.CData(2, :)   = [0.4, 0.8, 0.2];
hold on;
s3                  = swarmchart(ones(size(PPC(2, :))), PPC(2, :), ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.9, ...
    'MarkerEdgeColor', [0.4, 0.4, 0.4], ...
    'MarkerEdgeAlpha', 0.5);
s3.SizeData         = 5;
s3.XJitter          = 'rand';
s3.XJitterWidth     = 0.4;
s4                  = swarmchart(ones(size(PPC(3, :))) * 2, PPC(3, :), ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.9, ...
    'MarkerEdgeColor', [0.4, 0.8, 0.2], ...
    'MarkerEdgeAlpha', 0.5);
s4.SizeData         = 5;
s4.XJitter          = 'rand';
s4.XJitterWidth     = 0.4;
ylabel('PPC');
ylim([0, 0.2]);
box off;
set(gca, 'tickdir', 'out');
print(ppcCompFig, fullfile(paths.save, 'PPC100vs1000'), '-dsvg', '-r300');

%% Rayleigh z-value

% mean and SEM of Rayleigh z-value
meanZ           = mean(zval, 2);
semZ            = std(zval, [], 2) / sqrt(size(zval, 2));

% plot Rayleigh z-values for 100 and 1000 spikes
rayFig          = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
b1               = bar([meanZ(2), meanZ(3)], 'FaceColor', 'flat', 'FaceAlpha', 0.6);
set(gca, 'XTickLabel', {'100 spikes'; '1000 spikes'});
b1.CData(1, :)   = [0.4, 0.4, 0.4];
b1.CData(2, :)   = [0.4, 0.8, 0.2];
hold on;
s1                  = swarmchart(ones(size(zval(2, :))), zval(2, :), ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.9, ...
    'MarkerEdgeColor', [0.4, 0.4, 0.4], ...
    'MarkerEdgeAlpha', 0.5);
s1.SizeData         = 5;
s1.XJitter          = 'rand';
s1.XJitterWidth     = 0.4;
s2                  = swarmchart(ones(size(zval(3, :))) * 2, zval(3, :), ...
    'k', 'filled', ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.9, ...
    'MarkerEdgeColor', [0.4, 0.8, 0.2], ...
    'MarkerEdgeAlpha', 0.5);
s2.SizeData         = 5;
s2.XJitter          = 'rand';
s2.XJitterWidth     = 0.4;
ylabel('Rayleigh z-value');
ylim([0, 80]);
box off;
set(gca, 'tickdir', 'out');
print(rayFig, fullfile(paths.save, 'Rayleigh100vs1000'), '-dsvg', '-r300');

%% correlation Rayleigh z-value and PPC (1000 spikes)

% correlation
[rhoCorr1000, pvalCorr1000] = corr(PPC(3, :)', zval(3, :)', 'type', 'Spearman');

% plot
corrFig                     = figure('units', 'centimeters', 'position', [15, 15, 8, 6]);
plot(PPC(3, :), zval(3, :), 'x', 'Color', [0.64, 0.88, 0.52]);
set(gca, 'tickdir', 'out', 'box', 'off');
xlabel('PPC');
ylabel('Rayleigh z-value');
xlim([0, 0.08]);
ylim([0, 80]);
print(corrFig, fullfile(paths.save, 'PPCvsRayleigh'), '-dsvg', '-r300');


