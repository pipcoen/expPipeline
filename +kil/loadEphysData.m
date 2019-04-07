function [eph] = loadEphysData(x, kilosortOutput)
if ~exist('kilosortOutput', 'var'); kilosortOutput = x.kilosortOutput; end
%% Load timeline and associated inputs
if ~exist(x.rawTimeline, 'file'); error('No timeline file exists for requested ephys session'); end
fprintf('Loading timeline... \n');
timeline = load(x.rawTimeline); timeline = timeline.Timeline;

inputNames = {timeline.hw.inputs.name}';
% Get acqLive signal
acqLiveTrace = timeline.rawDAQData(:,strcmp(inputNames, 'acqLive'));
acqLiveTrace = acqLiveTrace>(max(acqLiveTrace)/2);
% acqLiveSrtEnd = timeline.rawDAQTimestamps(find(diff(acqLiveTrace))+1);

% Get wheel position %%%%NEVER USED?
% wheelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
% wheelPosition(wheelPosition > 2^31) = wheelPosition(wheelPosition > 2^31) - 2^32;

%Get whether stim was flickering %%%%NEVER USED?
% if contains('stimScreen', inputNames); stimScreenFlicker = range(timeline.rawDAQData(:,strcmp('stimScreen', inputNames))) > 2; end

% Get flipper flips
flipperTrace = timeline.rawDAQData(:,strcmp(inputNames, 'flipper')) > 2;
flipperFlips = sort([strfind(flipperTrace', [0 1]), strfind(flipperTrace', [1 0])])'+1;
flipperFlipTimesTimeline = timeline.rawDAQTimestamps(flipperFlips)';

%% Load task/behavior
fprintf('Loading block file... \n');
block  = load(x.processedData, 'blk'); block  = block.blk;

%% Load ephys data (single long recording)
fprintf('Loading ephys... \n');
% These are the digital channels going into the FPGA
acqLiveSyncIdx = 2;
flipperSyncIdx = 4;
%%
% Load clusters and sync/photodiode
clusterGroups = tdfread([kilosortOutput '\cluster_group.tsv']);
sync = load([kilosortOutput '\sync.mat']); sync = sync.sync;

%%
% Read header information
headerFID = fopen([kilosortOutput '\dat_params.txt']);
headerInfo = textscan(headerFID,'%s %s', 'delimiter',{' = '});
fclose(headerFID);
headerInfo = [headerInfo{1}'; headerInfo{2}'];
header = struct(headerInfo{:});
%%
% Load spike data
ephysSampleRate = str2double(header.apSampleRate);
spikeTimes = double(readNPY([kilosortOutput '\spike_times.npy']))./ephysSampleRate;
spikeTemplates = readNPY([kilosortOutput '\spike_templates.npy']);
templates = readNPY([kilosortOutput '\templates.npy']);
channelPositions = readNPY([kilosortOutput '\channel_positions.npy']);
channelMap = readNPY([kilosortOutput '\channel_map.npy']); %#ok<NASGU>
winv = readNPY([kilosortOutput '\whitening_mat_inv.npy']);
templateAmplitudes = readNPY([kilosortOutput '\amplitudes.npy']);

% Default channel map/positions are from end: make from surface
channelPositions(:,2) = max(channelPositions(:,2)) - channelPositions(:,2);


if exist('flipperFlipTimesTimeline','var') && length(sync) >= 4
    % Get flipper experiment differences by long delays
    flipThresh = 1; % time between flips to define experiment gap (s)
    flipTimes = sync(flipperSyncIdx).timestamps;
    flipperStEnIdx = [[1;find(diff(flipTimes) > flipThresh)+1], [find(diff(flipTimes) > flipThresh); length(flipTimes)]];
    experimentDurations = diff(flipTimes(flipperStEnIdx),[],2);
    [~, currExpIdx] = min(abs(experimentDurations-timeline.rawDAQTimestamps(end)));
    flipperFlipTimesFPGA = flipTimes(flipperStEnIdx(currExpIdx,1):flipperStEnIdx(currExpIdx,2));
   
    % Check that number of flipper flips in timeline matches ephys
    if length(flipperFlipTimesFPGA) ~= length(flipperFlipTimesTimeline)
        fprintf([x.subject ' ' x.expDate ': WARNING = Flipper flip times different in timeline/ephys \n']);
    else, spikeTimesTimeline = interp1(flipperFlipTimesFPGA,flipperFlipTimesTimeline,spikeTimes,'linear','extrap');
    end
    
    % Get the spike/lfp times in timeline time (accounts for clock drifts)
end
%%
if ~exist('spikeTimesTimeline', 'var')
    warning([x.subject ' ' x.expDate ': WARNING = Could not sync flipper. Using AcqLive... \n']);
    % Get acqLive times for current experiment
    acqLiveTimes = sync(acqLiveSyncIdx).timestamps;
    acqLiveStEnIdx = [find(sync(acqLiveSyncIdx).values == 1) find(sync(acqLiveSyncIdx).values == 0)];
    experimentDurations = diff(acqLiveTimes(acqLiveStEnIdx),[],2);
    [~, currExpIdx] = min(abs(experimentDurations-timeline.rawDAQTimestamps(end)));
    acqLiveTimesFPGA = acqLiveTimes(acqLiveStEnIdx(currExpIdx,1):acqLiveStEnIdx(currExpIdx,2));
    acqLiveTimesTimeline = timeline.rawDAQTimestamps([find(acqLiveTrace,1),find(acqLiveTrace,1,'last')+1]);
    
    % Check that the experiment time is the same within threshold
    % (it should be almost exactly the same)
    if abs(diff(acqLiveTimesFPGA) - diff(acqLiveTimesTimeline)) > 1
        error([x.subject ' ' x.expDate ': acqLive duration different in timeline and ephys']);
    else, spikeTimesTimeline = interp1(acqLiveTimesFPGA,acqLiveTimesTimeline,spikeTimes,'linear','extrap');
    end
end

%%

% Get the depths of each template
[spikeAmps, ~, templateDepths, ~, ~, templateDuration, waveforms] = ...
    kil.templatePositionsAmplitudes(templates,winv,channelPositions(:,2),spikeTemplates,templateAmplitudes);
% Eliminate spikes that were classified as not "good"    
fprintf('Removing noise and MUA templates... \n');

goodTemplatesList = clusterGroups.cluster_id(contains(num2cell(clusterGroups.group,2), 'good '));
goodTemplatesIdx = ismember(0:size(templates,1)-1,goodTemplatesList);


% Throw out all non-good template data
% templates = templates(goodTemplatesIdx,:,:);
templateDepths = templateDepths(goodTemplatesIdx);
waveforms = waveforms(goodTemplatesIdx,:);
templateDuration = templateDuration(goodTemplatesIdx);
%%
% Throw out all non-good spike data
goodSpikeIdx = ismember(spikeTemplates,goodTemplatesList);
goodSpikeIdx(spikeTimesTimeline<-10 | spikeTimesTimeline>timeline.rawDAQTimestamps(end)+10) = 0;
spikeTemplates = spikeTemplates(goodSpikeIdx);
spikeTimesTimeline = spikeTimesTimeline(goodSpikeIdx);
spikeAmps = spikeAmps(goodSpikeIdx);

% Re-name the spike templates according to the remaining templates
% (and make 1-indexed from 0-indexed)
newSpikeIdx = nan(max(spikeTemplates)+1,1);
newSpikeIdx(goodTemplatesList+1) = 1:length(goodTemplatesList);
spikeTemplates = newSpikeIdx(spikeTemplates+1);

%%
fields2copy = {'subject'; 'expDate'; 'expNum'; 'expDef'; 'kilosortOutput'};
for i = 1:length(fields2copy); eph.(fields2copy{i}) = x.(fields2copy{i}); end
eph.recordingSiteIdx = single(spikeTimesTimeline);
eph.spikeTimes = single(spikeTimesTimeline);
eph.spikeAmps = single(spikeAmps);
eph.clusterID = uint16(spikeTemplates);
eph.clusterDepths = templateDepths;
eph.clusterDuration = templateDuration;
eph.waveforms = waveforms;
eph = prc.catStructs(eph, prc.filtStruct(x.aligned, x.validTrials));
%%
fprintf('Finished loading experiment... \n');








