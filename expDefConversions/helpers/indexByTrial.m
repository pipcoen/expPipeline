function sortedByTrial = indexByTrial(blk, prmTimes, prmValues, substractStartTime)
if ~exist('substractStartTime', 'var'); substractStartTime =0*prmValues(1,:); end
trialStart = blk.trialStart;
trialEnd = blk.trialEnd;

trialIdx = arrayfun(@(x) x+0*prmTimes(prmTimes>=trialStart(x) & prmTimes<=trialEnd(x)), 1:length(trialEnd), 'uni', 0)';
trialIdx = cell2mat(trialIdx);

sortedByTrial = cell(max(trialIdx),1);
for i = 1:length(sortedByTrial)
    subtractValues = repmat(substractStartTime, sum(trialIdx==i), 1);
    if any(substractStartTime); subtractValues = blk.stimPeriodStart(i)*subtractValues; end

    sortedByTrial{i} = prmValues(trialIdx==i,:);
    if isempty(sortedByTrial{i}); continue; end
    sortedByTrial{i} = sortedByTrial{i} - subtractValues;
end
end