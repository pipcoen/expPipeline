function sortedByTrial = indexByTrial(blk, prmTimes, prmValues, substractStartTime)
if ~exist('substractStartTime', 'var'); substractStartTime =0*prmValues(1,:); end
trialStart = blk.trialStart;
trialEnd = blk.trialEnd;

[eventCount, ~, trialIdx] = histcounts(prmTimes, sort([trialStart;trialEnd+realmin]));
outOfBounds = trialIdx==0;
prmValues(outOfBounds,:) = []; trialIdx(outOfBounds) = [];

idxBounds = [find(diff([-10;trialIdx])>0) find(diff([trialIdx;1e6])>0)];
idxBounds(mod(unique(trialIdx),2)==0,:) = [];
eventCount(2:2:end) = [];

uniqueIdx = unique(trialIdx(mod(trialIdx,2)>0));
sortedByTrial = cell(length(uniqueIdx),1);
for i = 1:length(sortedByTrial)
    subtractValues = repmat(substractStartTime, eventCount(i), 1);
    if any(substractStartTime); subtractValues = blk.stimPeriodStart(i)*subtractValues; end

    sortedByTrial{i} = prmValues(idxBounds(i,1):idxBounds(i,2),:);
    if isempty(sortedByTrial{i}); continue; end
    sortedByTrial{i} = sortedByTrial{i} - subtractValues;
end
end