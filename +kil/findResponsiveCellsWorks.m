function [clusterSigLevel, clusterTTestData] = findResponsiveCells(blk,eventTimes,tWin)
%Assigning default values
sessionIDs = double(unique(blk.ephClusterSession));
if ~exist('eventTimes', 'var')
    switch blk.expDef
        case 'multiSpaceWorldPassive'
            eventTimes = [nanmean([blk.ephVisStimPeriodOnOff(:,1) blk.ephAudStimPeriodOnOff(:,1)],2); blk.ephRewardTimes];
            eventTimes = [eventTimes repmat(blk.sessionIdx,2,1)];
            eventTimes = eventTimes(~isnan(eventTimes(:,1)),:);
            if ~exist('tWin', 'var'); tWin = [-0.5 -0.1 0 0.5]; end
        case 'multiSpaceWorld'
            eventTimes = [nanmean([blk.ephVisStimPeriodOnOff(:,1) blk.ephAudStimPeriodOnOff(:,1)],2) blk.sessionIdx];
            eventTimes = eventTimes(~isnan(eventTimes(:,1)),:);
            if ~exist('tWin', 'var'); tWin = [-1 -0.1 0 1]; end
    end
elseif size(eventTimes,2) == 1 && size(eventTimes,1) == length(blk.sessionIdx); eventTimes = [eventTimes blk.sessionIdx];
elseif size(eventTimes,2) == 1 && length(unique(sessionIDs)) == 1; eventTimes = [eventTimes eventTimes*0+1];
else error('Could not figure out session info for event times');
end

if ~exist('tWin', 'var'); tWin = [-0.5 -0.1 0 0.5];
elseif numel(tWin)~=4; error('tWin should be 1x4 vector');
elseif ~all([(tWin(1:2)<=0) (tWin(3:4)>=0)]); error('pre/post windows should be negative/positive (or zero)');
end

%% Note: this looping is faster because it means smaller subsets are indexed when searching.
clusterSigLevel = nan*ones(length(blk.ephClusterAmps),1);
clusterTTestData = cell(length(blk.ephClusterAmps),1);
for i = sessionIDs'
    selEvents = eventTimes(i == eventTimes(:,2));
    sessClust = unique(blk.ephSpikeCluster(i == blk.ephSpikeSession));
    ttestData = zeros(length(sessClust), length(selEvents));
    for j = 1:length(selEvents)
        idx2Take = i == blk.ephSpikeSession & blk.ephSpikeTimes<(selEvents(j)+tWin(4)) & blk.ephSpikeTimes>(selEvents(j)+tWin(1));
        spikeTimes = blk.ephSpikeTimes(idx2Take)-selEvents(j);
        spikeCluster = blk.ephSpikeCluster(idx2Take);
        spikeCounts = cell2mat(arrayfun(@(x) [sum(x==spikeCluster(spikeTimes<tWin(2))) sum(x==spikeCluster(spikeTimes>tWin(3)))], sessClust, 'uni', 0));
        diffSpikeRates = diff(spikeCounts./[range(tWin(1:2)) range(tWin(3:4))], [], 2);
        ttestData(:,j) = diffSpikeRates;
    end
    [~, pVal] = ttest(ttestData');
    clusterTTestData(sessClust) = num2cell(ttestData',1)';
    clusterSigLevel(sessClust) = pVal;
end
end



