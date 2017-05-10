function [blk, prm] = orientTuning_Blk_Proc(~, b, blk, prm)
pVals = b.paramsValues;
gSTD = [pVals(1).gratingSF pVals(1).gratingTF pVals(1).stimulusDuration];
gOri = [pVals.gratingOrient]';
vCon = pVals(1).visualContrast;

stmV = sigOnOffTimes(b.events.sPreValues, b.events.sPreTimes);
blk.stmV = indexByTrial(blk, stmV(:,1), stmV, [1,1]);

whTV = indexByTrial(blk, blk.wrTV(:,1), blk.wrTV, [1,2]);

blk.gSTD = gSTD;
blk.gOri = gOri;
blk.vCon = vCon;
blk.whTV = whTV;
end