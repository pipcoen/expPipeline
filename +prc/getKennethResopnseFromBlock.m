function updatedBlk = getKennethResopnseFromBlock(origBlk)
if any(origBlk.tri.outcome.responseMade==0)
    error('Remove timeout trials before getting Kenneth times');
end
threshVal = 0.5;

blk = origBlk;
idxL = blk.tri.outcome.responseMade==1;
idxR = blk.tri.outcome.responseMade==2;

evalPnts = 0:0.01:0.5;
rawWheelTV = blk.tri.raw.wheelTimeValue;
rawVisAzi = cellfun(@(x) x(x(:,1)>0.5,:), blk.tri.raw.visAzimuthTimeValue, 'uni', 0);
wheelAtVizAzi = cellfun(@(x,y) interp1(x(:,1), x(:,2), y(:,1), 'nearest', 'extrap'), rawWheelTV, rawVisAzi, 'uni', 0);


allAzi = cell2mat(cellfun(@(x) diff(x(:,2)), rawVisAzi, 'uni', 0));
allWheel = cell2mat(cellfun(@diff, wheelAtVizAzi, 'uni', 0));

alignTimes = repmat({evalPnts}, blk.tot.trials,1);
wheelPos = cellfun(@(x,y) interp1(x(:,1), x(:,2), [y y(end)], 'nearest', 'extrap')*-1, rawWheelTV, alignTimes, 'uni', 0);
wheelPos = cell2mat(cellfun(@(x) x(1:end-1), wheelPos, 'uni', 0));
wheelPosL = normalize(wheelPos(idxL,:), 'scale');
wheelPosR = normalize(wheelPos(idxR,:), 'scale');

wheelLMean = nanmean(wheelPosL);
wheelRMean = nanmean(wheelPosR);

kMoveL = bsxfun(@minus,  wheelPos, wheelLMean).^2;
kMoveR = bsxfun(@minus,  wheelPos, wheelRMean).^2;
moveTime = cell2mat(cellfun(@(x) max([find(abs(x)>threshVal,1),nan]), num2cell(kMoveL-kMoveR,2), 'uni', 0));
moveTime(~isnan(moveTime)) = evalPnts(moveTime(~isnan(moveTime)));
origBlk.tri.outcome.timeToKenMove = moveTime;
updatedBlk = origBlk;
end