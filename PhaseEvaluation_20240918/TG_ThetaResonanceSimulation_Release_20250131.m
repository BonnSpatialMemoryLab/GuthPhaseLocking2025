%==========================================================================
% This script checks, whether theta resonance of individual neurons leads
% to theta-phase locking to aperiodic (1/f²) signals.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths                   = struct();
paths.save              = 'D:\TreasureHunt\PhaseEvaluation_20240918';

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% parameters
params.duration         = 3600;     % duration of the signal in seconds
params.fs               = 2000;     % sampling rate in Hz
params.alpha            = 2;        % slope of 1/f aperiodic signal
params.freqRange        = [1, 10];  % frequency range
params.nSim             = 100;      % number of simulated signals
params.thetaFreqs       = 1:0.1:10; % frequencies of 'resonant' neurons
params.exampleDur       = 10;       % in seconds
params.nSur             = 101;      % number of surrogates  

% set random seed
randSeed    = 444;
rng(randSeed, 'twister');
randSeedNum = randi(100000, params.nSim, 1);

% loop through simulations
allPval     = cell(params.nSim, 1);
for iSim = 1:params.nSim

    % display progress
    disp(iSim);

    % set random seed
    rng(randSeedNum(iSim));

    % number of samples
    nSamples                        = params.duration * params.fs;

    % time vector
    timeVector                      = (0:nSamples - 1) / params.fs;

    % frequency vector
    nPos                            = (nSamples / 2) + 1;
    frequencyVector                 = (0:(nPos - 1)) * (params.fs / nSamples);

    % define power for each frequency
    powerSpectrum                   = zeros(1, nPos);
    validFreqs                      = frequencyVector >= params.freqRange(1, 1) & frequencyVector <= params.freqRange(1, 2);
    powerSpectrum(validFreqs)       = 1 ./ (frequencyVector(validFreqs) .^ params.alpha); % apply 1/f

    % random phases for each positive-frequency bin
    randomPhases                    = exp(2 * pi * 1i * rand(size(frequencyVector)));

    % create Fourier coefficients with random phases
    fourierCoeffsPos                = powerSpectrum .* randomPhases;

    % create negative side by conjugating indices
    fourierCoeffsNeg                = conj(fourierCoeffsPos(nPos-1:-1:2));

    % create a spectrum
    fullSpectrum                    = [fourierCoeffsPos, fourierCoeffsNeg];

    % inverse FFT to get the time-domain signal
    signal                          = real(ifft(fullSpectrum));

    % normalize the signal
    signal                          = (signal - mean(signal)) / std(signal);

    % compute phase with the generalized phase approach
    signalPhase                     = TG_generalized_phase_vector(signal', params.fs, 1);
    signalAngles                    = angle(signalPhase);

    % simulate spike trains with theta resonance
    thisPval                        = nan(size(params.thetaFreqs'));
    for iTrain = 1:size(params.thetaFreqs, 2)

        % simulate spike train
        spikeTrain                              = false(1, size(signal, 2));
        stepSize                                = params.fs / params.thetaFreqs(iTrain);
        currentIndex                            = 1; % start at first sample
        while currentIndex <= size(signal, 2)
            spikeTrain(1, round(currentIndex))  = true; % place a spike at the rounded index
            currentIndex                        = currentIndex + stepSize; % Increment by exact step size
        end

        % get spike phases
        spikePhases                             = signalAngles(spikeTrain, 1);

        % get PPC
        thisPPC                                 = TG_PPC_20241128(spikePhases);

        % create surrogate spike phases
        surSpikePhases                          = 2 * pi * rand(size(spikePhases, 1), params.nSur) - pi;

        % compute surrogate PPCs
        surPPC      = nan(params.nSur, 1);
        parfor iSur = 1:params.nSur
            surPPC(iSur, 1)                              = TG_PPC_20241128(surSpikePhases(:, iSur));
        end

        % p-value
        thisPval(iTrain, 1)                     = 1 - (sum(thisPPC > surPPC) / params.nSur);

        % plot example
        if iSim == 1 && iTrain == 11

            % plot the time-domain signal and spike train
            timeFig     = figure;
            plot(timeVector(timeVector < (params.exampleDur + 0.1)), signal(timeVector < (params.exampleDur + 0.1)), 'Color', 'k', 'LineWidth', 2);
            hold on;
            xlim([0, params.exampleDur]);
            xlabel('Time (s)');
            ylabel('Amplitude');
            ylim([-5, 5]);
            yticks([-5, 0, 5]);
            set(gca, 'TickDir', 'out', 'box', 'off');
            title('Simulated aperiodic Signal');
            set(gcf, 'Renderer', 'painter');

            % plot spike train
            spikeTimes              = find(spikeTrain) / params.fs;
            spikeTimesToPlot        = spikeTimes(spikeTimes <= (params.exampleDur + 0.1));
            xData                   = [spikeTimesToPlot; spikeTimesToPlot];
            yData                   = repmat([-4, -3], length(spikeTimesToPlot), 1)';
            plot(xData, yData, 'Color', 'k', 'LineWidth', 2);

            % save figure
            saveas(timeFig, fullfile(paths.save, 'SimulatedThetaResonance_AperiodicSignal.svg'));

            % plot the frequency-domain signal
            fourierSignal           = fft(signal);
            powerSpectrum           = (abs(fourierSignal) .^ 2) / nSamples;
            fourierFreqs            = (0:(nSamples - 1)) * (params.fs / nSamples);
            freqFig                 = figure;
            plot(fourierFreqs, powerSpectrum, 'Color', 'k', 'LineWidth', 2);
            set(gca, 'TickDir', 'out', 'box', 'off', 'XScale', 'log', 'YScale', 'log');
            xlim([0.9, 11]);
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title('Power of simulated aperiodic LFP Signal');
            set(gcf, 'Renderer', 'painter');
            saveas(freqFig, fullfile(paths.save, 'SimulatedThetaResonance_AperiodicSignalPowerSpectrum.svg'));

            % plot in polar histogram
            polarFig    = figure;
            polarhistogram(spikePhases, 24, 'FaceColor', 'k');
            set(gca, 'ThetaTickLabel', {'0'; ''; ''; '90'; ''; ''; '\pm180'; ''; ''; '-90'; ''; ''});
            rlim([0, 400]);
            set(gcf, 'Renderer', 'painter');
            title(['Simulated phases P = ', num2str(thisPval(iTrain, 1))]);
            saveas(polarFig, fullfile(paths.save, 'SimulatedThetaResonance_Phases.svg'));
        end
    end

    % collect p-values
    allPval{iSim, 1} = thisPval;
end

% concatenate p-values
allPval     = cat(1, allPval{:});

% find significant part
sigRatio    = sum(allPval < 0.05) / size(allPval, 1);

% plot histogram with P-values
pvalFig     = figure;
histogram(allPval, 20, 'FaceColor', 'k', 'Normalization', 'probability');
xline(0.05, 'Color', 'k', 'LineStyle', ':');
xlabel('P-value');
ylabel('Probability');
xticks([0, 0.5, 1]);
yticks([0, 0.1, 0.2, 0.3]);
ylim([0, 0.3]);
title(['P-values of simulated cells - significant: ', num2str(sigRatio)]);
set(gca, 'TickDir', 'out', 'box', 'off');
set(gcf, 'Renderer', 'painter');
saveas(pvalFig, fullfile(paths.save, 'SimulatedThetaResonance_Pvalues.svg'));
