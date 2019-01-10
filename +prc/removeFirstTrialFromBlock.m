function [truncatedBlock, truncatedGalvo] = removeFirstTrialFromBlock(block, galvoLog)
fieldList = fieldnames(block.events);
trial2Idx = max(find(block.events.repeatNumValues==1,2));
trial2Start = block.events.endTrialTimes(trial2Idx-1)+0.0001;
for i = 1:2:length(fieldList)
    if contains(fieldList{i}, {'rTotValues'; 'totalRewardValues'}); continue; end 
    if isempty(block.events.(fieldList{i})); continue; end 
    if length(block.events.(fieldList{i+1})) == 1; continue; end 
    cutOff = find(block.events.(fieldList{i+1})<trial2Start, 1, 'last');
    eventRatio = size(block.events.(fieldList{i}),2)/size(block.events.(fieldList{i+1}),2);
    if round(eventRatio) ~= eventRatio; keyboard; end
    block.events.(fieldList{i})(:,1:(cutOff*eventRatio)) = [];
    block.events.(fieldList{i+1})(:,1:cutOff) = [];
end
%%
fieldList = fieldnames(block.outputs);
for i = 1:2:length(fieldList)
    if isempty(block.outputs.(fieldList{i})); continue; end 
    if length(block.outputs.(fieldList{i+1})) == 1; continue; end 
    cutOff = find(block.outputs.(fieldList{i+1})<trial2Start, 1, 'last');
    eventRatio = size(block.outputs.(fieldList{i}),2)/size(block.outputs.(fieldList{i+1}),2);
    if round(eventRatio) ~= eventRatio; keyboard; end
    block.outputs.(fieldList{i})(:,1:(cutOff*eventRatio)) = [];
    block.outputs.(fieldList{i+1})(:,1:cutOff) = [];
end
%%

block.events.trialNumValues = block.events.trialNumValues-trial2Idx+1;
block.paramsTimes(1:(trial2Idx-1)) = [];
block.paramsValues(1:(trial2Idx-1)) = [];
if isstruct(galvoLog) && length(fieldnames(galvoLog))>1
    galvoIdx = find(galvoLog.trialNum>=trial2Idx, 1);
    fieldList = fieldnames(galvoLog);
    for i = 1:length(fieldList)
        if strcmp(fieldList{i}, 'stereotaxCalib'); continue; end
        galvoLog.(fieldList{i})(1:(galvoIdx-1),:) = [];
    end
    galvoLog.trialNum = galvoLog.trialNum-trial2Idx+1;
end
truncatedGalvo = galvoLog;
truncatedBlock = block;
end


