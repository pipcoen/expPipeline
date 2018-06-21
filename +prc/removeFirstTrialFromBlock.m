function x = removeFirstTrialFromBlock(x)
fieldList = fieldnames(x.oldBlock.events);
trial2Idx = max(find(x.oldBlock.events.repeatNumValues==1,2));
trial2Start = x.oldBlock.events.endTrialTimes(trial2Idx-1)+0.0001;
for i = 1:2:length(fieldList)
    if contains(fieldList{i}, {'rTotValues'; 'totalRewardValues'}); continue; end 
    if isempty(x.oldBlock.events.(fieldList{i})); continue; end 
    if length(x.oldBlock.events.(fieldList{i+1})) == 1; continue; end 
    cutOff = find(x.oldBlock.events.(fieldList{i+1})<trial2Start, 1, 'last');
    eventRatio = size(x.oldBlock.events.(fieldList{i}),2)/size(x.oldBlock.events.(fieldList{i+1}),2);
    if round(eventRatio) ~= eventRatio; keyboard; end
    x.oldBlock.events.(fieldList{i})(:,1:(cutOff*eventRatio)) = [];
    x.oldBlock.events.(fieldList{i+1})(:,1:cutOff) = [];
end
x.oldBlock.events.trialNumValues = x.oldBlock.events.trialNumValues-trial2Idx+1;
x.oldBlock.paramsTimes(1:(trial2Idx-1)) = [];
x.oldBlock.paramsValues(1:(trial2Idx-1)) = [];
if isstruct(x.galvoLog) && length(fieldnames(x.galvoLog))>1
    galvoIdx = find(x.galvoLog.trialNum>=trial2Idx, 1);
    fieldList = fieldnames(x.galvoLog);
    for i = 1:length(fieldList)
        if strcmp(fieldList{i}, 'stereotaxCalib'); continue; end
        x.galvoLog.(fieldList{i})(1:(galvoIdx-1),:) = [];
    end
    x.galvoLog.trialNum = x.galvoLog.trialNum-trial2Idx+1;
end
end


