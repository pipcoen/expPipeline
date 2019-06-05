function [clusterSigLevel, clusterTTestData] = findResponsiveCells(blk,eventTimes,tWin)
%Assigning default values
sessionIDs = double(unique(blk.ephClusterSession));
if ~exist('eventTimes', 'var') || isempty(eventTimes)
    switch blk.expDef
        case 'multiSpaceWorldPassive'
            eventTimes = [nanmean([blk.ephVisStimPeriodOnOff(:,1) blk.ephAudStimPeriodOnOff(:,1)],2); blk.ephRewardTimes];
            eventTimes = [eventTimes repmat(blk.sessionIdx,2,1)];
            eventTimes = eventTimes(~isnan(eventTimes(:,1)),:);
        case 'multiSpaceWorld'
            eventTimes = [nanmean([blk.ephVisStimPeriodOnOff(:,1) blk.ephAudStimPeriodOnOff(:,1)],2) blk.sessionIdx];
            eventTimes = eventTimes(~isnan(eventTimes(:,1)),:);
    end
elseif size(eventTimes,2) == 1 && size(eventTimes,1) == length(blk.sessionIdx); eventTimes = [eventTimes blk.sessionIdx];
elseif size(eventTimes,2) == 1 && length(unique(sessionIDs)) == 1; eventTimes = [eventTimes eventTimes*0+sessionIDs];
else, error('Could not figure out session info for event times');
end

if ~exist('tWin', 'var') || isempty(tWin)
    switch blk.expDef
        case 'multiSpaceWorldPassive'; tWin = [-0.5 -0.1 0.05 0.25];
        case 'multiSpaceWorld'; tWin = [-0.5 -0.1 0.05 0.25];
    end
elseif numel(tWin)~=4; error('tWin should be 1x4 vector');
elseif ~all([(tWin(1:2)<=0) (tWin(3:4)>=0)]); error('pre/post windows should be negative/positive (or zero)');
end

%% Note: this looping is faster because it means smaller subsets are indexed when searching.
clusterSigLevel = nan*ones(length(blk.ephClusterAmps),1);
clusterTTestData = cell(length(blk.ephClusterAmps),1);
for i = sessionIDs'
    selEvents = sort(eventTimes(eventTimes(:,2)==i,1));
    currSessionIdx = i == blk.ephSpikeSession;
    sessionSpikeTimes = blk.ephSpikeTimes(currSessionIdx);
    sessionSpikeCluster = blk.ephSpikeCluster(currSessionIdx);
    sessClusters  = unique(sessionSpikeCluster);
    
    if ~isempty(selEvents)
        eventWindows = selEvents+tWin;
        spikeCounts = histcounts2(sessionSpikeTimes, sessionSpikeCluster, sort(eventWindows(:)), 1:max(sessClusters)+1);
        spikeCountsPre = spikeCounts(1:4:end,sessClusters)./range(tWin(1:2));
        spikeCountsPost = spikeCounts(3:4:end,sessClusters)./range(tWin(3:4));
        ttestData = spikeCountsPost - spikeCountsPre;
        
        [~, pVal] = ttest(ttestData);
        clusterTTestData(sessClusters) = num2cell(ttestData,1)';
        clusterSigLevel(sessClusters) = pVal';
    end
end
end



