function clusterSigLevel = findResponsiveCells(blk,eventTimes,tWin)
%Assigning default values
if ~exist('tWin', 'var'); tWin = [-0.5 -0.1 0 0.5];
elseif numel(tWin)~=4; error('tWin should be 1x4 vector');
elseif ~all([(tWin(1:2)<=0) (tWin(3:4)>=0)]); error('pre/post windows should be negative/positive (or zero)');
end

if ~exist('eventTimes', 'var')
    switch blk.expDef
        case 'multiSpaceWorldPassive'
            eventTimes = [nanmean([blk.ephVisStimPeriodOnOff(:,1) blk.ephAudStimPeriodOnOff(:,1)],2); blk.ephRewardTimes];
            eventTimes = [eventTimes repmat(blk.sessionIdx,2,1)];
            eventTimes = eventTimes(~isnan(eventTimes(:,1)),:);
    end
elseif size(eventTimes,1) == 1; eventTimes = [eventTimes eventTimes*0+1];
end

%% Note: this looping is faster because it means smaller subsets are indexed when searching.
clusterSigLevel = nan*ones(length(blk.ephClusterAmps),1);
for i = unique(blk.ephSessionIdx)'
    selEvents = eventTimes(i == eventTimes(:,2));
    sessClust = unique(blk.ephClusterID(i == blk.ephSessionIdx));
    ttestData = zeros(length(sessClust), length(selEvents));
    for j = 1:length(selEvents)
        idx2Take = i == blk.ephSessionIdx & blk.ephSpikeTimes<(selEvents(j)+tWin(4)) & blk.ephSpikeTimes>(selEvents(j)+tWin(1));
        spikeTimes = blk.ephSpikeTimes(idx2Take)-selEvents(j);
        clusterID = blk.ephClusterID(idx2Take);
        spikeCounts = cell2mat(arrayfun(@(x) [sum(x==clusterID(spikeTimes<tWin(2))) sum(x==clusterID(spikeTimes>tWin(3)))], sessClust, 'uni', 0));
        diffSpikeRates = diff(spikeCounts./[range(tWin(1:2)) range(tWin(3:4))], [], 2);
        ttestData(:,j) = diffSpikeRates;
    end
    [~, pVal] = ttest(ttestData');
    clusterSigLevel(sessClust) = pVal;
end
end



