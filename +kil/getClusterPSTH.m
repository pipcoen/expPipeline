function [clusterPSTH, timeValues] = getClusterPSTH(blk,eventTimes,tBins)
%Assigning default values
experimentIdxs = cell2mat(blk.subExpPenLink(:,2));
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
elseif size(eventTimes,2) == 1 && size(eventTimes,1) == length(blk.trialStartEnd); eventTimes = [eventTimes blk.trialExperimentIdx];
elseif size(eventTimes,2) == 1 && length(unique(experimentIdxs)) == 1; eventTimes = [eventTimes eventTimes*0+experimentIdxs];
elseif size(eventTimes,2) == 1, error('Could not figure out experiment info for event times');
end

%% Note: this looping is faster because it means smaller subsets are indexed when searching.
clusterPSTH = cell(length(blk.eph_clusterAmps),1);
for i = cell2mat(blk.subExpPenLink(:,3))'
    tDat = prc.filtStruct(blk, i, 'penetration');
    if ~exist('tWin', 'var') || isempty(tBins)
        tBins = -0.5:0.01:0.5;
    end
    timeValues = tBins(1:end-1) + diff(tBins)./2;
    
    selEvents = sort(eventTimes(eventTimes(:,2)==tDat.subExpPenLink{:,2},1));
    clusters = 1:length(tDat.eph_clusterAmps);
    
    if ~isempty(selEvents)
        eventWindows = selEvents+tBins;
        spikeCounts = histcounts2(tDat.eph_spikeTimes, tDat.eph_spikeCluster, [sort(eventWindows(:));1e6], 1:(clusters(end)+1));
        spikeCounts = permute(reshape(spikeCounts, [size(eventWindows') length(clusters)]),[2 1 3]);
        clusterPSTH(blk.eph_clusterPenetrationIdx==i) = num2cell(spikeCounts(:,1:end-1,:),[1 2]);
    end
end
end



