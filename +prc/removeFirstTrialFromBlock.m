function x = removeFirstTrialFromBlock(x)
fieldList = fieldnames(x.oldBlock.events);
trial2Idx = max(find(x.oldBlock.events.repeatNumValues==1,2));
trial2Start = x.oldBlock.events.endTrialTimes(trial2Idx-1)+0.001;
for i = 3:2:length(fieldList)
    x.oldBlock.events.(fieldList{i})(:,x.oldBlock.events.(fieldList{i+1})<trial2Start) = [];
    x.oldBlock.events.(fieldList{i+1})(:,x.oldBlock.events.(fieldList{i+1})<trial2Start) = [];
end
x.oldBlock.events.trialNumValues = x.oldBlock.events.trialNumValues-trial2Idx+1;
x.oldBlock.paramsTimes(1:(trial2Idx-1)) = [];
x.oldBlock.paramsValues(1:(trial2Idx-1)) = [];
if isstruct(x.galvoLog)
    galvoIdx = find(x.galvoLog.trialNum>=trial2Idx, 1);
    fieldList = fieldnames(x.galvoLog);
    for i = 1:length(fieldList)-1
        x.galvoLog.(fieldList{i})(1:(galvoIdx-1),:) = [];
    end
    x.galvoLog.trialNum = x.galvoLog.trialNum-galvoIdx+1;
end
end


