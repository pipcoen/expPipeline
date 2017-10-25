function [newBlock, newParams] = orientTuning_Blk_Proc(x)
v = x.standardizedBlock.paramsValues;
e = x.standardizedBlock.events;
n = x.newBlock;
p = x.standardizedParams;

n.gratingSpatialFreq = v(1).gratingSF;
n.gratingTemporalFreq = v(1).gratingTF;
n.stimulusDuration = v(1).stimulusDuration;
n.gratingOrientation = [v.gratingOrient]';
n.visContrast = v(1).visContrast;

n.stimPeriodOnOff = indexByTrial(n, e.stimPeriodOnOffTimes', [e.stimPeriodOnOffTimes' e.stimPeriodOnOffValues']);
n.wheelTimeValue = indexByTrial(n, n.rawWheelTimeValue(:,1), n.rawWheelTimeValue);

newParams = p;
newBlock = n;
end