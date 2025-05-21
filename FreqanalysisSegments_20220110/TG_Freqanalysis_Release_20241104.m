%==========================================================================
% This script performs power analysis for LFP data. It compares power
% during successful and unsuccessful encoding and recall segments
% using cluster-based permutation testing.
%
% Tim Guth, 2025
%==========================================================================

%% settings
clc; close all; clear;
rng(444);

% paths
paths           = struct();
paths.behlog    = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.SUData    = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.LFPData   = 'D:\TreasureHunt\MicroDownsampled_20210910'; % LFP data
paths.artifact  = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.subInfo   = 'D:\TreasureHunt\SubjectInformation_20210111'; % subject information
paths.phaseRes  = 'D:\TreasureHunt\PhaseAnalysis_20230921\1_10_Hz_Results'; % phase analysis folder
paths.save      = 'D:\TreasureHunt\FreqanalysisSegments_20220110\Freqanalysis_20241104'; % save folder
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
params               = struct();
params.freqMin       = 1; % minimum frequency to analyze
params.freqMax       = 40; % maximum frequency to analyze
params.freqNum       = 40; % number of frequencies of interest
params.buffer        = 5; % borders before and after segment in seconds
params.segDur        = 4; % duration of each segment to cut in seconds
params.timeRes       = 0.01; % length of time bins for TFR analysis
params.condition     = {'encoding'; 'objectRecall'; 'locationRecall'};
params.c4aName       = {'Cluster4Analysis_LK_20200416_113912', 'Cluster4Analysis_TG_20210315_020120'};

% for random locations in the Treasure Hunt arena to normalize drop error
RandLoccfg           = struct();
RandLoccfg.maxR      = 50;
RandLoccfg.minR      = 0;
RandLoccfg.N         = 1000000;
RandLoccfg.centerX   = 370;
RandLoccfg.centerY   = 360;
randLocs             = TG_RandomPointsInCircle(RandLoccfg);

% for ft_freqanalysis
TFRcfg               = [];
TFRcfg.method        = 'wavelet';
TFRcfg.pad           = 'nextpow2';
TFRcfg.output        = 'pow';
TFRcfg.keeptrials    = 'yes';
TFRcfg.width         = 5; % number of wavelet cycles
TFRcfg.foi           = round(logspace(log10(params.freqMin), log10(params.freqMax), params.freqNum), 2); % see Miller et al. 2018 Nature Communications
TFRcfg.toi           = -params.buffer:params.timeRes:params.segDur + params.buffer;

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

%% load phase results to get list of included wires
phaseRes            = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

%% wire index
unitIdx             = cat(1, phaseRes.phaseRes.idx);
[wireIdx, selIdx]   = unique(unitIdx(:, 1:3), 'rows'); % all included wires

%% preallocations

% main results
allRes = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % get sessions
    sessions         = dir(fullfile(paths.LFPData, subjects{iSub}, 'session*'));
    
    %% get each channel's brain region

    % load information about brain region
    micro2regionFile = dir(fullfile(paths.subInfo, subjects{iSub}, 'Microwires', 'Micro2Region.txt'));
    fileID           = fopen(fullfile(micro2regionFile.folder, micro2regionFile.name));
    micro2Region     = textscan(fileID, '%s %s');
    fclose(fileID);
    
    % load information about left or right hemisphere
    micro2macroFile  = dir(fullfile(paths.subInfo, subjects{iSub}, 'Microwires', 'Micro2Macro.txt'));
    fileID           = fopen(fullfile(micro2macroFile.folder, micro2macroFile.name));
    micro2macro      = textscan(fileID, '%*s %s');
    fclose(fileID);
    
    % channel list
    chanNumbers      = cellfun(@str2num, micro2Region{:, 1}, 'UniformOutput', false);
    chanList         = cat(2, chanNumbers{:, 1});
    
    % information about hemisphere
    hemisphereShort  = cellfun(@(x) x(end), micro2macro{:}, 'UniformOutput', false);
    hemisphere       = replace(hemisphereShort, {'R'; 'L'; 'd'}, {'right'; 'left'; 'NotImplanted'});
    
    % create list of each channel's brain region
    brainRegionsList       = cell(size(chanList, 2), 2);
    brainRegionsList(:, 2) = num2cell(chanList');
    regionNames            = [];
    for iMacro = 1:size(chanNumbers, 1)
        regionNames = cat(1, regionNames, repmat(strcat(micro2Region{1, 2}(iMacro, 1), '_', hemisphere{iMacro, 1}), numel(chanNumbers{iMacro, 1}), 1));
    end
    brainRegionsList(:, 1) = regionNames;
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % original session index
        sessIdx = split(sessions(iSess).name, '_');
        sessIdx = str2double(sessIdx{2});
        
        % display session information
        fprintf('\n=============================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % get available data
        LFPChanDir = TG_GetChanDir_20210812(paths.LFPData, subjects{iSub}, sessions(iSess).name); % LFP data
        
        % brain region sanity check
        if ~isequal(size(LFPChanDir, 1), size(brainRegionsList, 1))
            error('different number of microwire channels than channels in the brain region file');
        end

        %% chest-wise memory performance

        % load behavioral logfile and segment info file
        trialInfo     = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo     = trialInfo.trialInfo;
        segmentInfo   = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'segmentInfo_20230208.mat'));
        
        % each chest's Treasure Hunt trial index
        encTrialIdx   = transpose(rude(trialInfo.numChestsPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numChestsPerTrial))));
        objTrialIdx   = transpose(rude(trialInfo.numObjRecallPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numObjRecallPerTrial))));
        locTrialIdx   = transpose(rude(trialInfo.numLocRecallPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numLocRecallPerTrial))));
        
        % index for encoding segments
        encSegmentIdx = trialInfo.chestIdx;
        
        % prelloacte - object recall
        encObjRec     = cell(size(encTrialIdx, 1), 1);
        objSegmentIdx = nan(size(trialInfo.audioResponse, 1), 1);
        
        % prelloacte - location recall
        encLocRec     = nan(size(encTrialIdx, 1), 1);
        locRec        = nan(size(locTrialIdx, 1), 1);
        locSegmentIdx = nan(size(locTrialIdx, 1), 1);
        
        % loop through chests
        for iChest = 1:size(encTrialIdx, 1)
            
            % object recall
            if strcmp(trialInfo.recallType{encTrialIdx(iChest)}, 'OBJECT')
                thisChestLoc = trialInfo.chestLoc(iChest, :);
                idxObjRec    = abs((trialInfo.cueingLoc(:, 1) - thisChestLoc(1))) < 0.0001 ...
                    & abs(trialInfo.cueingLoc(:, 3) - thisChestLoc(3)) < 0.0001 ...
                    & encTrialIdx(iChest) == objTrialIdx;
                encObjRec(iChest, 1)        = trialInfo.audioResponse(idxObjRec);
                objSegmentIdx(idxObjRec, 1) = iChest;
            
            % location recall
            elseif strcmp(trialInfo.recallType{encTrialIdx(iChest)}, 'LOCATION')
                thisChestLabel = trialInfo.TREASURE_LABEL_EN(iChest, 1);
                idxLocRec      = strcmp(trialInfo.cueingObject, thisChestLabel) & encTrialIdx(iChest) == locTrialIdx;
                if sum(idxLocRec == 1)
                    D                           = pdist2(trialInfo.CORRECT_TEST_POSITION(idxLocRec, [1, 3]), trialInfo.CHOSEN_TEST_POSITION(idxLocRec, [1, 3]));
                    DSurro                      = pdist2(trialInfo.CORRECT_TEST_POSITION(idxLocRec, [1, 3]), randLocs);
                    encLocRec(iChest, 1)        = sum(D < DSurro) / numel(DSurro);
                    locRec(idxLocRec, 1)        = sum(D < DSurro) / numel(DSurro);
                    locSegmentIdx(idxLocRec, 1) = iChest;
                else
                    warning('No location recall for %s, session %d, chest %d', subjects{iSub}, iSess, iChest);
                end
            end
        end
        
        % memory performance
        bGoodMemEnc     = strcmp(encObjRec, 'CORRECT') | encLocRec > median(encLocRec, 'omitnan'); % or: > 0.9; % logical index for good memory performance during encoding
        bGoodMemEncAbs  = strcmp(encObjRec, 'CORRECT') | encLocRec > 0.9;
        bGoodMemObj     = strcmp(trialInfo.audioResponse, 'CORRECT'); % logical index for successful object recall
        bGoodMemLoc     = locRec > median(locRec, 'omitnan'); % or '> 0.9'; % logical index for good location recall performance
        bGoodMemLocAbs  = locRec > 0.9;

        % concatenate object and location recall
        bGoodMemRec     = cat(1, bGoodMemObj, bGoodMemLoc);
        bGoodMemRecAbs  = cat(1, bGoodMemObj, bGoodMemLocAbs);
        recTrialIdx     = cat(1, objTrialIdx, locTrialIdx);
        recSegmentIdx   = cat(1, objSegmentIdx, locSegmentIdx);
        
        % classify object (1) and location recall (2) segments
        encObjOrLoc                = zeros(size(encSegmentIdx, 1), 1);
        encObjOrLoc(objSegmentIdx(~isnan(objSegmentIdx))) = 1;
        encObjOrLoc(locSegmentIdx) = 2;
        recObjOrLoc                = [ones(size(bGoodMemObj, 1), 1); ones(size(bGoodMemLoc, 1), 1) * 2];
        
        % preallocate
        allChanRes = cell(size(LFPChanDir, 1), 1);
        
        %% loop through channels
        parfor iWire = 1:size(LFPChanDir, 1)

            % brain region sanity check
            if ~isequal(str2double(LFPChanDir(iWire).name(5:end)), brainRegionsList{iWire, 2})
                error('mismatch between microwire channel and channel in the brain region file');
            end
            
            % continue if there were no analyzed units from this wire
            thisWireIdx = cat(2, iSub, sessIdx, iWire);

            if ~ismember(thisWireIdx, wireIdx, 'rows')
                fprintf('Excluded channel.\n');
                continue;
            end
            
            %% load data and perform time-frequency analysis

            % load data
            LFPData         = load(fullfile(paths.LFPData, subjects{iSub}, sessions(iSess).name, LFPChanDir(iWire).name, 'datacutFt2000Hz.mat'));
            
            % load artifact data
            artifactData    = load(fullfile(paths.artifact, subjects{iSub}, sessions(iSess).name, LFPChanDir(iWire).name, 'detectedArtifacts.mat'));
            bArtifact       = artifactData.bArtifact;

            % loop through encoding and recall
            segmentsOrig    = [];
            chanRes         = [];
            for iCond = 1:size(params.condition, 1)

                % include only encoding or recall periods
                if contains(params.condition{iCond}, 'encoding')
                    segmentsOrig    = segmentInfo.encoding;
                elseif contains(params.condition{iCond}, 'objectRecall')
                    segmentsOrig    = segmentInfo.objRecall;
                elseif contains(params.condition{iCond}, 'locationRecall')
                    segmentsOrig    = segmentInfo.locRecall;
                end

                % define data to cut
                cfg                     = [];
                cfg.trl                 = round(segmentsOrig / (segmentInfo.fsample / LFPData.fsample)); % adjust sampling rate of segmentInfo to data
                if contains(params.condition{iCond}, 'locationRecall')
                    cfg.trl(:, 1)            = cfg.trl(:, 2) - (params.segDur * LFPData.fsample) - (params.buffer * LFPData.fsample);
                else
                    cfg.trl(:, 1)            = cfg.trl(:, 1) - (params.buffer * LFPData.fsample);
                end
                cfg.trl(:, 2)            = cfg.trl(:, 2) + (params.buffer * LFPData.fsample);
                cfg.trl(:, 3)            = -(params.buffer * LFPData.fsample);

                % cut data
                allSegmentData           = ft_redefinetrial(cfg, LFPData);

                % loop through trials and set artifact data to NaN
                for iSeg = 1:size(cfg.trl, 1)
                    bArtSeg                                 = bArtifact(cfg.trl(iSeg, 1):cfg.trl(iSeg, 2));
                    allSegmentData.trial{1, 1}(1, bArtSeg)  = NaN;
                end
                
                % time-frequency analysis
                TFR_allData              = ft_freqanalysis(TFRcfg, allSegmentData);
                
                % log-transformation
                TFR_allData.powspctrm    = log10(TFR_allData.powspctrm);
                
                % reshape and normalize TFR
                squeezedData             = permute(squeeze(TFR_allData.powspctrm), [3, 1, 2]);
                reshData1                = reshape(squeezedData, size(squeezedData, 1) * size(squeezedData, 2), size(squeezedData, 3));
                normData                 = normalize(reshData1, 1);
                reshData2                = reshape(normData, size(squeezedData, 1), size(squeezedData, 2), size(squeezedData, 3));
                permData                 = permute(reshData2, [2, 4, 3, 1]);
                TFR_allData.powspctrm    = permData; % z-transformed powerspectrum
                
                % change axis values to 'trick' fieldtrip to plot all frequencies with same height
                TFR_allData.freq         = 1:numel(TFR_allData.freq); % linearly increasing vector to 'trick' fieldtrip to plot all frequencies with same height
                
                % collect results across channels
                chanRes.idx              = [iSub, sessIdx, iWire];
                chanRes.subject          = subjects{iSub};
                chanRes.session          = sessions(iSess).name;
                chanRes.channel          = LFPChanDir(iWire).name;
                chanRes.brainRegion      = brainRegionsList{iWire, 1};
                
                % collect all results for each condition
                if strcmp(params.condition{iCond}, 'encoding')
                    chanRes.encBothPowspctrm = squeeze(mean(TFR_allData.powspctrm, 'omitnan'));
                    chanRes.encSuccPowspctrm = squeeze(mean(TFR_allData.powspctrm(bGoodMemEnc, :, :, :), 'omitnan'));
                    chanRes.encFailPowspctrm = squeeze(mean(TFR_allData.powspctrm(~bGoodMemEnc, :, :, :), 'omitnan'));
                elseif strcmp(params.condition{iCond}, 'objectRecall')
                    chanRes.objBothPowspctrm = squeeze(mean(TFR_allData.powspctrm, 'omitnan'));
                    chanRes.objSuccPowspctrm = squeeze(mean(TFR_allData.powspctrm(bGoodMemObj, :, :, :), 'omitnan'));
                    chanRes.objFailPowspctrm = squeeze(mean(TFR_allData.powspctrm(~bGoodMemObj, :, :, :), 'omitnan'));
                elseif strcmp(params.condition{iCond}, 'locationRecall')
                    chanRes.locBothPowspctrm = squeeze(mean(TFR_allData.powspctrm, 'omitnan'));
                    chanRes.locSuccPowspctrm = squeeze(mean(TFR_allData.powspctrm(bGoodMemLoc, :, :, :), 'omitnan'));
                    chanRes.locFailPowspctrm = squeeze(mean(TFR_allData.powspctrm(~bGoodMemLoc, :, :, :), 'omitnan'));
                end
            end
            
            % collect results across channels
            allChanRes{iWire, 1} = chanRes;
            
        end
        
        % collect all results across subjects
        allRes = cat(2, allRes, [allChanRes{:}]);
    end
end

% shut down parallel pool
delete(gcp);

% save average results for each channel
save(fullfile(paths.save, 'results'), '-v7.3');

%% additional analysis and figures

%==========================================================================
% in this part of the script power effects during successful and
% failed encoding and object/location recall are compared
%==========================================================================

% load previously saved files
r = load(fullfile(paths.save, 'results'));

% configuration for random number generator
rng(444);

% configurations for ft_freqstatistics
FScfg                     = [];
FScfg.channel             = {'merged_channel'};
FScfg.parameter           = 'powspctrm';
FScfg.method              = 'montecarlo';
FScfg.statistic           = 'ft_statfun_depsamplesT';
FScfg.correctm            = 'cluster';
FScfg.clusterstatistic    = 'maxsum';
FScfg.correcttail         = 'alpha';
FScfg.neighbours          = [];
FScfg.tail                = 0;
FScfg.clustertail         = 0;
FScfg.clusteralpha        = 0.05;
FScfg.alpha               = 0.05;
FScfg.numrandomization    = 1001;
FScfg.ivar                = 1;
FScfg.uvar                = 2;

% configurations for ft_singleplot
Plotcfg                   = [];
Plotcfg.figure            = 'gcf';
Plotcfg.interactive       = 'no';
Plotcfg.xlim              = [-1, 5];
Plotcfg.zlim              = [-0.05, 0.05];
newYLabels                = 2 .^ (ceil((log10(min(r.TFRcfg.foi)) / log10(2))):floor((log10(max(r.TFRcfg.foi)) / log10(2))));  % create new y tick values
Plotcfg.newYLabels        = num2cell(newYLabels);
Plotcfg.yticks            = (log10(newYLabels) - log10(r.params.freqMin) ...
                            + ((log10(r.params.freqMax) - log10(r.params.freqMin)) / (r.params.freqNum - 1))) ...
                            * ((r.params.freqNum - 1) / (log10(r.params.freqMax) - log10(r.params.freqMin))); % create new y ticks positions
Plotcfg.parameter         = 'powspctrm';

% preallocate
biggestClusStat           = cell(size(r.params.condition, 1), 2);

%% loop through encoding, object recall and location recall
for iCond = 1:size(r.params.condition, 1)
    
    % condition-specific settings/variables
    if strcmp(r.params.condition{iCond}, 'encoding')
        succPow             = cat(3, r.allRes.encSuccPowspctrm);
        failPow             = cat(3, r.allRes.encFailPowspctrm);
        FScfg.latency       = [0, 1.5];
    elseif strcmp(r.params.condition{iCond}, 'objectRecall')
        succPow             = cat(3, r.allRes.objSuccPowspctrm);
        failPow             = cat(3, r.allRes.objFailPowspctrm);
        FScfg.latency       = [0, 4];
    elseif strcmp(r.params.condition{iCond}, 'locationRecall')
        succPow             = cat(3, r.allRes.locSuccPowspctrm);
        failPow             = cat(3, r.allRes.locFailPowspctrm);
        FScfg.latency       = [0, 4];
    end

    % print number of channels
    fprintf('Total number of channels: %d.\n', size(r.allRes, 2));

    % latency start and end time index
    [~, latStartIdx]            = min(abs(r.TFRcfg.toi - FScfg.latency(1, 1)));
    [~, latEndIdx]              = min(abs(r.TFRcfg.toi - FScfg.latency(1, 2)));

    % define variables
    allSessSucc         = {};
    allSessFail         = {};
    numSess             = 0; % counter variable
    
    % loop through subjects
    for iSub = 1:size(r.subjects, 1)
        
        % get sessions
        sessions = dir(fullfile(r.paths.LFPData, r.subjects{iSub}, 'session*'));
        
        % loop through sessions
        for iSess = 1:size(sessions, 1)
            
            % get index for this session
            oneSessIdx               = strcmp({r.allRes.subject}, r.subjects{iSub}) & strcmp({r.allRes.session}, sessions(iSess).name);
            
            % remove session if there are no channels
            if oneSessIdx == 0
                continue;
            end
            
            % number of sessions
            numSess                  = numSess + 1;
            
            % get powerspectrum for successful and failed trials
            oneSessSucc              = mean(succPow(:, :, oneSessIdx), 3, 'omitnan');
            oneSessFail              = mean(failPow(:, :, oneSessIdx), 3, 'omitnan');
            
            % create fieldtrip structure for this session - successful
            oneSessSuccTFR           = [];
            oneSessSuccTFR.label     = {'merged_channel'};
            oneSessSuccTFR.dimord    = 'chan_freq_time';
            oneSessSuccTFR.freq      = 1:size(oneSessSucc, 1);
            oneSessSuccTFR.time      = r.TFRcfg.toi;
            oneSessSuccTFR.powspctrm = permute(oneSessSucc, [3, 1, 2]);
            
            % create fieldtrip structure for this session - fail
            oneSessFailTFR           = [];
            oneSessFailTFR.label     = {'merged_channel'};
            oneSessFailTFR.dimord    = 'chan_freq_time';
            oneSessFailTFR.freq      = 1:size(oneSessFail, 1);
            oneSessFailTFR.time      = r.TFRcfg.toi;
            oneSessFailTFR.powspctrm = permute(oneSessFail, [3, 1, 2]);
            
            % collect results across sessions
            allSessSucc{1, numSess}  = oneSessSuccTFR;
            allSessFail{1, numSess}  = oneSessFailTFR;
        end
    end
    
    %% grand average over sessions
    cfg         = [];
    cfg.nanmean = 'yes';
    grAvgSucc   = ft_freqgrandaverage(cfg, allSessSucc{:});
    grAvgFail   = ft_freqgrandaverage(cfg, allSessFail{:});
    
    %% ft_freqstatistics
    
    % design matrix for ft_freqstatistics
    FScfg.design             = [ones(1, size(allSessSucc, 2)), ones(1, size(allSessFail, 2)) * 2; ...
        1:size(allSessSucc, 2), 1:size(allSessFail, 2)];

    % apply ft_freqstatistics (permutation test)
    statSuccVsFail                  = ft_freqstatistics(FScfg, allSessSucc{:}, allSessFail{:});
    
    % collect time-frequency clusters with smallest p-values for each region
    biggestClusStat{iCond, 1}       = [statSuccVsFail.posclusters(1).prob, statSuccVsFail.posclusters(1).clusterstat, statSuccVsFail.posclusters(1).stddev, statSuccVsFail.posclusters(1).cirange];
    biggestClusStat{iCond, 2}       = [statSuccVsFail.negclusters(1).prob, statSuccVsFail.negclusters(1).clusterstat, statSuccVsFail.negclusters(1).stddev, statSuccVsFail.negclusters(1).cirange];
    
    % subtract powerspectra
    grAvgSuccVsFail                 = grAvgSucc;
    grAvgSuccVsFail.powspctrm       = grAvgSucc.powspctrm - grAvgFail.powspctrm;

    % add mask
    origMask                        = statSuccVsFail.mask;
    paddedMask                      = cat(3, false(1, 40, sum(grAvgSuccVsFail.time < FScfg.latency(1))), ...
        statSuccVsFail.mask, false(1, 40, sum(grAvgSuccVsFail.time > FScfg.latency(2))));
    grAvgSuccVsFail.mask            = paddedMask;

    % create TFR figure - difference
    TFR_Diff                = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
    Plotcfg.maskparameter   = 'mask';
    Plotcfg.maskstyle       = 'outline';
    ft_singleplotTFR(Plotcfg, grAvgSuccVsFail);
    yticks(Plotcfg.yticks); % change positions of y ticks
    set(gca, 'YTickLabel', Plotcfg.newYLabels); % change values of y ticks
    title(strcat('(', num2str(size(allSessSucc, 2)), {' '}, 'sessions)'), 'FontSize', 15);
    set(gca, 'fontsize', 15);
    if strcmp(r.params.condition{iCond}, 'encoding')
        xline([0, 1.5], ':', 'LineWidth', 2);
        xlim([-1, 2.5]);
    elseif strcmp(r.params.condition{iCond}, 'objectRecall')
        xline([0, 4], ':', 'LineWidth', 2);
    elseif strcmp(r.params.condition{iCond}, 'locationRecall')
        xline([0, 4], ':', 'LineWidth', 2);

        % change x-axis labels
        xticks                  = get(gca, 'XTick');
        xticklabels             = get(gca, 'XTickLabel');
        xticklabels_numeric     = str2double(xticklabels);
        new_xticklabels_numeric = xticklabels_numeric - 4;
        new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
        set(gca, 'XTickLabel', new_xticklabels);
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    % save figure
    set(gca, 'TickDir', 'out', 'box', 'off');
    set(gcf, 'Renderer', 'painters');
    saveas(TFR_Diff, fullfile(r.paths.save, strcat('TFR_diff_', r.params.condition{iCond}, '.svg')));

    % remove mask plot fields from Plotcfg
    Plotcfg = rmfield(Plotcfg, 'maskparameter');
    Plotcfg = rmfield(Plotcfg, 'maskstyle');
    
    % create TFR figure - successful
    TFR_Succ     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
    ft_singleplotTFR(Plotcfg, grAvgSucc);
    yticks(Plotcfg.yticks); % change positions of y ticks
    set(gca, 'YTickLabel', Plotcfg.newYLabels); % change values of y ticks
    title(strcat('(', num2str(size(allSessSucc, 2)), {' '}, 'sessions)'), 'FontSize', 15);
    set(gca, 'fontsize', 15);
    if strcmp(r.params.condition{iCond}, 'encoding')
        xline([0, 1.5], ':', 'LineWidth', 2);
        xlim([-1, 2.5]);
    elseif strcmp(r.params.condition{iCond}, 'objectRecall')
        xline([0, 4], ':', 'LineWidth', 2);
    elseif strcmp(r.params.condition{iCond}, 'locationRecall')
        xline([0, 4], ':', 'LineWidth', 2);

        % change x-axis labels
        xticks                  = get(gca, 'XTick');
        xticklabels             = get(gca, 'XTickLabel');
        xticklabels_numeric     = str2double(xticklabels);
        new_xticklabels_numeric = xticklabels_numeric - 4;
        new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
        set(gca, 'XTickLabel', new_xticklabels);
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % save figure
    set(gca, 'TickDir', 'out', 'box', 'off');
    set(gcf, 'Renderer', 'painters');
    saveas(TFR_Succ, fullfile(r.paths.save, strcat('TFR_succ_', r.params.condition{iCond}, '.svg')));
    
    % create TFR figure - fail
    TFR_Fail     = figure('units', 'centimeters', 'Position', [10, 10, 7, 8]);
    ft_singleplotTFR(Plotcfg, grAvgFail);
    yticks(Plotcfg.yticks); % change positions of y ticks
    set(gca, 'YTickLabel', Plotcfg.newYLabels); % change values of y ticks
    title(strcat('(', num2str(size(allSessFail, 2)), {' '}, 'sessions)'), 'FontSize', 15);
    set(gca, 'fontsize', 15);
    if strcmp(r.params.condition{iCond}, 'encoding')
        xline([0, 1.5], ':', 'LineWidth', 2);
        xlim([-1, 2.5]);
    elseif strcmp(r.params.condition{iCond}, 'objectRecall')
        xline([0, 4], ':', 'LineWidth', 2);
    elseif strcmp(r.params.condition{iCond}, 'locationRecall')
        xline([0, 4], ':', 'LineWidth', 2);

        % change x-axis labels
        xticks                  = get(gca, 'XTick');
        xticklabels             = get(gca, 'XTickLabel');
        xticklabels_numeric     = str2double(xticklabels);
        new_xticklabels_numeric = xticklabels_numeric - 4;
        new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
        set(gca, 'XTickLabel', new_xticklabels);
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % save figure
    set(gca, 'TickDir', 'out', 'box', 'off');
    set(gcf, 'Renderer', 'painters');
    saveas(TFR_Fail, fullfile(r.paths.save, strcat('TFR_fail_', r.params.condition{iCond}, '.svg')));
end
