function [clusterSigLevel, clusterTTestBlka] = findResponsiveCells(blk,eventTimes,tWin)
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
clusterSigLevel = nan*ones(blk.tot.clusters,1);
clusterTTestBlka = cell(blk.tot.clusters,1);

for i = 1:length(blk.pen.ephysRecordIdx)
    pentrationRef = blk.pen.ephysRecordIdx==blk.pen.ephysRecordIdx(i);
    tBlk = prc.filtBlock(blk, pentrationRef, 'pen');
    if ~exist('tWin', 'var') || isempty(tWin)
        switch blk.exp.expDef{1}
            case 'multiSpaceWorldPassive'; tWin = [-0.5 -0.01 0.01 0.5];
            case 'multiSpaceWorld'; tWin = [-0.5 -0.1 0.05 0.25];
        end
    elseif numel(tWin)~=4; error('tWin should be 1x4 vector');
    elseif ~all([(tWin(1:2)<=0) (tWin(3:4)>=0)]); error('pre/post windows should be negative/positive (or zero)');
    end
    
    selEvents = sort(eventTimes(eventTimes(:,2)==experimentIdxs(i)));
    cluIdx = (1:length(tBlk.clu.depths))';
    spkTimes = cell2mat(tBlk.clu.spkTimes);
    spkCluster = cell2mat(cellfun(@(x,y) x*0+y, tBlk.clu.spkTimes, num2cell(cluIdx), 'uni', 0));
    if ~isempty(selEvents)
        eventWindows = selEvents+tWin;
        spikeCounts = histcounts2(spkTimes, spkCluster, sort(eventWindows(:)), 1:(cluIdx(end)+1));
        spikeCountsPre = spikeCounts(1:4:end,cluIdx)./range(tWin(1:2));
        spikeCountsPost = spikeCounts(3:4:end,cluIdx)./range(tWin(3:4));
        ttestBlka = spikeCountsPost - spikeCountsPre;
        
        %         [~, pVal] = ttest(ttestBlka);
        pVal = zeros(1, size(ttestBlka,2));
        for k = 1:size(ttestBlka,2); pVal(1,k) = signrank(ttestBlka(:,k)); end
        clusterTTestBlka(blk.clu.penetrationRef==find(pentrationRef)) = num2cell(ttestBlka,1)';
        clusterSigLevel(blk.clu.penetrationRef==find(pentrationRef)) = pVal';
    end
end
end



