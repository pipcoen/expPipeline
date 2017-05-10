function [blk, prm] = multiTemporalPassive_Blk_Proc(~, b, blk, prm)
pVals = b.paramsValues;
ckAD = pVals(1).clickAmpDur';     %Click Duration and Rate
AVIx = [pVals.audVisIndex]';     %Click Duration and Rate
vCon = pVals(1).visualContrast; 

stmV = sigOnOffTimes(b.events.stmVValues, b.events.stmVTimes);
blk.stmV = indexByTrial(blk, stmV(:,1), stmV, [1,1]);
stmA = sigOnOffTimes(b.events.stmAValues, b.events.stmATimes);
blk.stmA = indexByTrial(blk, stmA(:,1), stmA(:,1), 1);

whTV = indexByTrial(blk, blk.wrTV(:,1), blk.wrTV, [1,2]);

blk.ckAD = ckAD;
blk.AVIx = AVIx;
blk.vCon = vCon;
blk.whTV = whTV;
end