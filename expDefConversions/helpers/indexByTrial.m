function sortedByTrial = indexByTrial(blk, prmTimes, prmValues, substractStartTime)
if ~exist('substractStartTime', 'var'); substractStartTime =0*prmValues(1,:); end
trialStart = blk.trialStart;
trialEnd = blk.trialEnd;

trialIdx = arrayfun(@(x) x+0*prmTimes(prmTimes>=trialStart(x) & prmTimes<=trialEnd(x))', 1:length(trialEnd), 'uni', 0)';
trialIdx = cell2mat(trialIdx);

sortedByTrial = cell(max(trialIdx),1);
for i = 1:length(sortedByTrial)
    tSub = (substractStartTime==1)*stim.sSrt(i);
    tVal = prmValues(tIdx==i,:);
    if isempty(tVal); continue; end
    tSub(substractStartTime==2) = tVal(1,substractStartTime==2);
    oPut{i} = single(bsxfun(@minus, tVal, tSub));
end
end