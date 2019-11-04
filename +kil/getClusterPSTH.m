function [clusterPSTH, timeValues] = getClusterPSTH(blk,eventTimes,tBins)
%Assigning default values
experimentIdxs = blk.pen.expRef;
if ~exist('eventTimes', 'var') || isempty(eventTimes)
    switch blk.expDef
        case 'multiSpaceWorldPassive'
            eventTimes = [nanmean([blk.eph_visStimPeriodOnOff(:,1) blk.eph_audStimPeriodOnOff(:,1)],2); blk.ephRewardTimes];
            eventTimes = [eventTimes repmat(blk.trialExperimentIdx,2,1)];
            eventTimes = eventTimes(~isnan(eventTimes(:,1)),:);
        case 'multiSpaceWorld'
            eventTimes = [nanmean([blk.eph_visStimPeriodOnOff(:,1) blk.eph_audStimPeriodOnOff(:,1)],2) blk.trialExperimentIdx];
            eventTimes = eventTimes(~isnan(eventTimes(:,1)),:);
    end
elseif size(eventTimes,2) == 1 && size(eventTimes,1) == blk.tot.trials; eventTimes = [eventTimes blk.tri.expRef];
elseif size(eventTimes,2) == 1 && length(unique(experimentIdxs)) == 1; eventTimes = [eventTimes eventTimes*0+experimentIdxs];
elseif size(eventTimes,2) == 1, error('Could not figure out experiment info for event times');
end

%% Note: this looping is faster because it means smaller subsets are indexed when searching.
clusterPSTH = cell(blk.tot.clusters,1);
for i = 1:length(blk.pen.ephysRecordIdx)
    pentrationRef = blk.pen.ephysRecordIdx==blk.pen.ephysRecordIdx(i);
    tBlk = prc.filtBlock(blk, pentrationRef, 'pen');
    if ~exist('tWin', 'var') || isempty(tBins)
        tBins = -0.5:0.01:0.5;
    end
    timeValues = tBins(1:end-1) + diff(tBins)./2;
    
    selEvents = sort(eventTimes(eventTimes(:,2)==experimentIdxs(i)));
    cluIdx = (1:length(tBlk.clu.depths))';
    spkTimes = cell2mat(tBlk.clu.spkTimes);
    spkCluster = cell2mat(cellfun(@(x,y) x*0+y, tBlk.clu.spkTimes, num2cell(cluIdx), 'uni', 0));
    
    if ~isempty(selEvents)
        eventWindows = selEvents+tBins;
        spikeCounts = histcounts2(spkTimes, spkCluster, [sort(eventWindows(:));1e6], 1:(cluIdx(end)+1));
        spikeCounts = permute(reshape(spikeCounts, [size(eventWindows') length(cluIdx)]),[2 1 3]);
        clusterPSTH(blk.clu.penetrationRef==find(pentrationRef)) = num2cell(spikeCounts(:,1:end-1,:),[1 2]);
    end
end
end



