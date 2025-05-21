%==========================================================================
% This script finds oscillations in the downsampled LFP using Bycycle
% (Cole & Voytek, Journal of Neurophysiology, 2019) on the peaks and
% troughs identified with the generalized phase approach.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.lfpData   = 'D:\TreasureHunt\MicroFiltered_20210930'; % data path
paths.phaseData = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.save      = 'D:\TreasureHunt\MicroBycycle_20241216'; % save folder

% add functions
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% parameters
params                                  = struct();
params.filterBand                       = [1, 10]; % frequency band
params.filterBandStr                    = strjoin(strsplit(num2str(params.filterBand)), '_');
params.phaseDataName                    = strcat('datacutFt2000HzBP_generalized_', params.filterBandStr, '.mat');
params.method                           = 'bycycle'; % phase extraction method
params.amplitudeFractionThreshold       = 0.2; % normalized amplitude threshold
params.amplitudeConsistencyThreshold    = 0.3; % difference in the decay and rise voltage within a cycle
params.periodConsistencyThreshold       = 0.5; % difference between a cycle's period and the period of the adjacent cycles
params.monotonicityThreshold            = 0.6; % fraction of monotonic voltage changes in decay and rise phases
params.minCyclesThreshold               = 2; % minimum number of cycles

% load and save name
loadName            = strcat('datacutFt2000HzBP_', regexprep(num2str(params.filterBand), '\s+', '_'), '.mat');
saveName            = strcat('datacutFt2000HzBP_', params.method, '_', regexprep(num2str(params.filterBand), '\s+', '_'), '.mat');

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

%% save settings
save(fullfile(paths.save, 'settings'));

%% loop through subjects
for iSub = 1:length(subjects)

    % get sessions
    sessions = dir(fullfile(paths.lfpData, subjects{iSub}, 'session*'));

    % loop through sessions
    for iSess = 1:size(sessions, 1)

        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);

        % get available microwires
        chanDir = TG_GetChanDir_20210812(paths.lfpData, subjects{iSub}, sessions(iSess).name);

        % loop through channels
        for iWire = 1:size(chanDir, 1)

            % print channel name
            fprintf('\tChannel name: %s.\n', chanDir(iWire).name);

            % channel directory
            thisChanDir                     = fullfile(paths.lfpData, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name);

            %% load data

            % load LFP data
            data                            = load(fullfile(thisChanDir, loadName));
            lfpData                         = data.trial{1, 1};

            % sampling rate
            fsample                         = data.fsample;

            % load phase data
            try
                phaseData                   = load(fullfile(paths.phaseData, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, params.phaseDataName));
                thisPhaseData               = phaseData.trial{1, 1};
                if size(unique(thisPhaseData), 2) <= 1
                    continue;
                end
            catch
                continue;
            end

            % identify troughs
            bTrough                         = diff(angle(thisPhaseData)) < -6; % find transitions from pi to -pi (troughs)
            troughIdx                       = [1, bwlabel(~bTrough)];
            troughIdx(troughIdx == 0)       = troughIdx(find(troughIdx == 0) - 1);
            [~, troughSmp, newTroughIdx]    = unique(troughIdx);
            troughSmp                       = troughSmp(2:end);

            % identify peaks
            shiftedPhases                   = circ_dist(angle(thisPhaseData), pi); % shift signal by pi
            bPeak                           = diff(shiftedPhases) < -6; % find transitions from pi to -pi in shifted signal (peaks)
            peakIdx                         = [1, bwlabel(~bPeak)];
            peakIdx(peakIdx == 0)           = peakIdx(find(peakIdx == 0) - 1);
            [~, peakSmp, newPeakIdx]        = unique(peakIdx);
            peakSmp                         = peakSmp(2:end);

            % samples for peaks and troughs
            peakTroughSmp                   = [peakSmp; troughSmp];
            [peakTroughSmp, sortIdx]        = sort(peakTroughSmp);

            % index for peaks and troughs
            peakTroughIdx                   = [ones(size(peakSmp)); 2 * ones(size(troughSmp))];
            peakTroughIdx                   = peakTroughIdx(sortIdx);

            % remove trough or peak repeats
            labeledRepeats                  = [bwlabel(diff(peakTroughIdx) == 0); 0];
            for iRepeat = 1:max(labeledRepeats)

                % find repeat
                thisRepeatIdx       = [find(labeledRepeats == iRepeat); find(labeledRepeats == iRepeat, 1, 'last') + 1];

                % get repeat samples
                thisRepeatSamples   = peakTroughSmp(thisRepeatIdx);

                % get repeat values
                thisRepeatValues    = lfpData(thisRepeatSamples);

                % get index of included peak/trough
                if peakTroughIdx(thisRepeatIdx(1)) == 1
                    [~, keepThis]       = max(thisRepeatValues);
                else
                    [~, keepThis]       = min(thisRepeatValues);
                end

                % get samples to exclude
                thisRepeatIdx(keepThis)         = [];

                % remove samples
                peakTroughIdx(thisRepeatIdx)    = [];
                peakTroughSmp(thisRepeatIdx)    = [];
                labeledRepeats(thisRepeatIdx)   = [];
            end

            % get "cleaned" troughs and peaks
            peakSmp                     = peakTroughSmp(peakTroughIdx == 1);
            troughSmp                   = peakTroughSmp(peakTroughIdx == 2);

            % get values
            peakValues                  = lfpData(peakSmp);
            troughValues                = lfpData(troughSmp);

            % cycle-by-cycle frequencies
            cycleSamples                = diff(peakSmp);
            cycleFreq                   = phaseData.fsample ./ cycleSamples;

            %% run Bycycle (on flipped data)

            % add python environment
            pyenv('Version', 'C:/Users/Tim Guth/anaconda3/envs/matlab_py38/python.exe');

            % import bycycle functions
            py.importlib.import_module('bycycle');
            py.importlib.import_module('pandas');

            % flip data (so that Bycycle detects cycles from peak to peak)
            flippedLfpData      = -lfpData;
            flippedPeaks        = troughSmp;
            flippedTroughs      = peakSmp;

            % make sure to start with a peak and end with a trough (for Bycycle)
            flippedTroughs      = flippedTroughs(flippedTroughs > flippedPeaks(1));
            flippedPeaks        = flippedPeaks(flippedPeaks < flippedTroughs(end));

            % convert data to python arrays
            lfpFilteredPy       = py.numpy.array(flippedLfpData);
            fsamplePy           = py.float(fsample);
            filterBandPy        = py.list(params.filterBand);

            % find rise and decay points using Bycycle
            outPy               = py.bycycle.cyclepoints.find_zerox(lfpFilteredPy, py.numpy.array(int64(flippedPeaks - 1)), py.numpy.array(int64(flippedTroughs - 1)));
            risesPy             = outPy{1};
            decayPy             = outPy{2};

            % extract rises and decays as Matlab array
            rises               = int64(py.array.array('d', py.numpy.nditer(risesPy)))';
            decays              = int64(py.array.array('d', py.numpy.nditer(decayPy)))';

            % construct Python dataframe with peaks, rises, decays, last troughs and next troughs
            dfSamplesPy         = py.pandas.DataFrame(pyargs('data', py.dict(pyargs( ...
                'sample_peak', py.numpy.array(int64(flippedPeaks(2:end) - 1)), ... % flipped data
                'sample_last_zerox_decay', py.numpy.array(int64(decays(1:end-1))), ...
                'sample_zerox_decay', py.numpy.array(int64(decays(2:end))), ...
                'sample_zerox_rise', risesPy, ...
                'sample_last_trough', py.numpy.array(int64(flippedTroughs(1:end-1) - 1)), ... % flipped data
                'sample_next_trough', py.numpy.array(int64(flippedTroughs(2:end) - 1)) ... % flipped data
                ))));

            % compute shape features
            outPy               = py.bycycle.features.shape.compute_durations(dfSamplesPy);
            periodPy            = outPy{1};
            timePeakPy          = outPy{2};
            timeTroughPy        = outPy{3};
            outPy               = py.bycycle.features.shape.compute_extrema_voltage(dfSamplesPy, lfpFilteredPy);
            voltPeakPy          = outPy{1};
            voltTroughPy        = outPy{2};
            symFeaturesPy       = py.bycycle.features.shape.compute_symmetry(dfSamplesPy, lfpFilteredPy, periodPy, timePeakPy, timeTroughPy);
            bandAmpPy           = py.bycycle.features.shape.compute_band_amp(dfSamplesPy, lfpFilteredPy, fsamplePy, filterBandPy);

            % shape features
            dfShapeFeaturesPy   = py.pandas.concat(py.tuple({dfSamplesPy, ...
                py.pandas.DataFrame(py.dict(pyargs( ...
                'period', periodPy, ...
                'time_peak', timePeakPy, ...
                'time_trough', timeTroughPy, ...
                'volt_peak', voltPeakPy, ...
                'volt_trough', voltTroughPy, ...
                'time_decay', symFeaturesPy{'time_decay'}, ...
                'time_rise', symFeaturesPy{'time_rise'}, ...
                'volt_decay', symFeaturesPy{'volt_decay'}, ...
                'volt_rise', symFeaturesPy{'volt_rise'}, ...
                'volt_amp', symFeaturesPy{'volt_amp'}, ...
                'time_rdsym', symFeaturesPy{'time_rdsym'}, ...
                'time_ptsym', symFeaturesPy{'time_ptsym'}, ...
                'band_amp', bandAmpPy ...
                )))}), pyargs('axis', int64(1)));

            % burst features
            dfBurstFeaturesPy   = py.bycycle.features.burst.compute_burst_features(dfShapeFeaturesPy, lfpFilteredPy, py.str('cycles'));

            % burst detection
            dfBurstFeaturesPy   = py.bycycle.burst.detect_bursts_cycles(dfBurstFeaturesPy, ...
                pyargs('amp_fraction_threshold', params.amplitudeFractionThreshold, ...
                'amp_consistency_threshold', params.amplitudeConsistencyThreshold, ...
                'period_consistency_threshold', params.periodConsistencyThreshold, ...
                'monotonicity_threshold', params.monotonicityThreshold, ...
                'min_n_cycles', params.minCyclesThreshold));

            % concatenate all results
            dfFeatures          = py.pandas.concat(py.tuple({dfShapeFeaturesPy, dfBurstFeaturesPy}), pyargs('axis', int64(1)));

            % loop through channel results
            df                  = dfFeatures.to_dict;
            resultsArray        = nan(double(py.len(df{'sample_peak'})), 22);
            for iRow = 1:(double(py.len(df{'sample_peak'})))
                resultsArray(iRow, 1)   = df{'sample_peak'}{iRow - 1};
                resultsArray(iRow, 2)   = df{'sample_zerox_decay'}{iRow - 1};
                resultsArray(iRow, 3)   = df{'sample_zerox_rise'}{iRow - 1};
                resultsArray(iRow, 4)   = df{'sample_last_trough'}{iRow - 1};
                resultsArray(iRow, 5)   = df{'sample_next_trough'}{iRow - 1};
                resultsArray(iRow, 6)   = df{'period'}{iRow - 1};
                resultsArray(iRow, 7)   = df{'time_peak'}{iRow - 1};
                resultsArray(iRow, 8)   = df{'time_trough'}{iRow - 1};
                resultsArray(iRow, 9)   = df{'volt_peak'}{iRow - 1};
                resultsArray(iRow, 10)  = df{'volt_trough'}{iRow - 1};
                resultsArray(iRow, 11)  = df{'time_rise'}{iRow - 1};
                resultsArray(iRow, 12)  = df{'volt_decay'}{iRow - 1};
                resultsArray(iRow, 13)  = df{'volt_rise'}{iRow - 1};
                resultsArray(iRow, 14)  = df{'volt_amp'}{iRow - 1};
                resultsArray(iRow, 15)  = df{'time_rdsym'}{iRow - 1};
                resultsArray(iRow, 16)  = df{'time_ptsym'}{iRow - 1};
                resultsArray(iRow, 17)  = df{'band_amp'}{iRow - 1};
                resultsArray(iRow, 18)  = df{'amp_fraction'}{iRow - 1};
                resultsArray(iRow, 19)  = df{'amp_consistency'}{iRow - 1};
                resultsArray(iRow, 20)  = df{'period_consistency'}{iRow - 1};
                resultsArray(iRow, 21)  = df{'monotonicity'}{iRow - 1};
                resultsArray(iRow, 22)  = df{'is_burst'}{iRow - 1};
            end

            % results table (adjusted labels to the flipped data)
            tableColumns                = {'sample_trough', 'sample_zerox_rise', 'sample_zerox_decay', ...
                'sample_last_peak', 'sample_next_peak', 'period','time_trough', ...
                'time_peak', 'volt_trough', 'volt_peak', 'time_decay', 'volt_rise', ...
                'volt_decay', 'volt_amp', 'time_drsym', 'time_tpsym', 'band_amp', ...
                'amp_fraction', 'amp_consistency', 'period_consistency', ...
                'monotonicity', 'is_burst'};
            bycycleTable                = array2table(resultsArray, 'VariableNames', tableColumns);

            % flip back peak and trough voltages
            bycycleTable.volt_peak      = -bycycleTable.volt_peak;
            bycycleTable.volt_trough    = -bycycleTable.volt_trough;

            % cycle frequency
            bycycleFreq                 = fsample ./ (bycycleTable.sample_next_peak - bycycleTable.sample_last_peak);

            % burst starts and ends
            burstStartIdx               = find(diff([0; bycycleTable.is_burst]) > 0);
            burstEndIdx                 = find(diff([bycycleTable.is_burst; 0]) < 0);

            % burst start and end samples
            burstStartSample            = bycycleTable.sample_last_peak(burstStartIdx);
            burstEndSample              = bycycleTable.sample_next_peak(burstEndIdx);

            % burst frequency
            burstCycles                 = burstEndIdx - burstStartIdx + 1;
            burstSamples                = burstEndSample - burstStartSample + 1;
            burstDuration               = burstSamples / fsample; % in seconds
            burstFrequency              = burstCycles ./ burstDuration;

            % create logical index for bursts
            bBurst                      = false(size(lfpData));
            dataBurstFreq               = nan(size(lfpData));
            for iBurst = 1:size(burstStartSample, 1)
                bBurst(burstStartSample(iBurst):burstEndSample(iBurst))         = true;
                dataBurstFreq(burstStartSample(iBurst):burstEndSample(iBurst))  = burstFrequency(iBurst);
            end

            % save results
            saveDir                     = fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name);
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end
            saveData                    = struct();
            saveData.bycycleTable       = bycycleTable;
            saveData.params             = params;
            saveData.fsample            = fsample;
            saveData.bBurst             = bBurst;
            saveData.dataBurstFreq      = dataBurstFreq;
            saveData.burstCycles        = burstCycles;
            saveData.burstDuration      = burstDuration;
            saveData.burstFrequency     = burstFrequency;
            TG_save(fullfile(saveDir, saveName), saveData);
        end
    end
end
