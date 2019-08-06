function [clusterSigLevel, clusterTTestData] = findResponsiveCells(blk,eventTimes,tWin)
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
clusterSigLevel = nan*ones(length(blk.eph_clusterAmps),1);
clusterTTestData = cell(length(blk.eph_clusterAmps),1);

for i = cell2mat(blk.subExpPenLink(:,3))'
    tDat = prc.filtStruct(blk, i, 'penetration');
    if ~exist('tWin', 'var') || isempty(tWin)
        switch blk.expDef{1}
            case 'multiSpaceWorldPassive'; tWin = [-0.5 -0.1 0.05 0.25];
            case 'multiSpaceWorld'; tWin = [-0.5 -0.1 0.05 0.25];
        end
    elseif numel(tWin)~=4; error('tWin should be 1x4 vector');
    elseif ~all([(tWin(1:2)<=0) (tWin(3:4)>=0)]); error('pre/post windows should be negative/positive (or zero)');
    end
    
    selEvents = sort(eventTimes(eventTimes(:,2)==tDat.subExpPenLink{:,2},1));
    clusters = 1:length(tDat.eph_clusterAmps);
    
    if ~isempty(selEvents)
        eventWindows = selEvents+tWin;
        spikeCounts = histcounts2(tDat.eph_spikeTimes, tDat.eph_spikeCluster, sort(eventWindows(:)), 1:(clusters(end)+1));
        spikeCountsPre = spikeCounts(1:4:end,clusters)./range(tWin(1:2));
        spikeCountsPost = spikeCounts(3:4:end,clusters)./range(tWin(3:4));
        ttestData = spikeCountsPost - spikeCountsPre;
        
        [~, pVal] = ttest(ttestData);
        clusterTTestData(blk.eph_clusterPenetrationIdx==i) = num2cell(ttestData,1)';
        clusterSigLevel(blk.eph_clusterPenetrationIdx==i) = pVal';
    end
end
end



