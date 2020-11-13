function grids = getGridsFromBlock(blk)
grids = prc.makeGrid(blk);

grids.conditions = prc.makeGrid(blk, blk.tri.stim.conditionLabel, @mean);
grids.performance = prc.makeGrid(blk, blk.tri.outcome.feedbackGiven==1, @mean,'abscondition');
grids.numTrials = prc.makeGrid(blk, blk.tri.outcome.responseMade==1, @length);
grids.numRightTurns = prc.makeGrid(blk, blk.tri.outcome.timeDirFirstMove(:,2)==2, @sum);
grids.fracRightTurns = prc.makeGrid(blk, blk.tri.outcome.timeDirFirstMove(:,2)==2, @nanmean);

[~,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), grids.numRightTurns, grids.numTrials, 'uni', 0);
grids.fracRightTurnsLowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
grids.fracRightTurnsHighBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));
grids.fracRightTurnsLowBound(isnan(grids.fracRightTurns)) = nan;
grids.fracRightTurnsHighBound(isnan(grids.fracRightTurns)) = nan;

grids.timeToFirstMove = prc.makeGrid(blk, blk.tri.outcome.timeDirFirstMove(:,1), @nanmedian, [], 1);
grids.timeToFirstMoveSE = nanstd(grids.timeToFirstMove,[],3)./sqrt(size(grids.timeToFirstMove,3));
grids.timeToFirstMove = nanmean(grids.timeToFirstMove, 3);

grids.timeChoiceCross = prc.makeGrid(blk, blk.tri.outcome.timeDirChoiceCross(:,1), @nanmedian, [], 1);
grids.timeChoiceCrossSE = nanstd(grids.timeChoiceCross,[],3)./sqrt(size(grids.timeChoiceCross,3));
grids.timeChoiceCross = nanmean(grids.timeChoiceCross, 3);

grids.timeToChoiceInit = prc.makeGrid(blk, blk.tri.outcome.timeDirChoiceInit(:,1), @nanmedian, [], 1);
grids.timeToChoiceInitSE = nanstd(grids.timeToChoiceInit,[],3)./sqrt(size(grids.timeToChoiceInit,3));
grids.timeToChoiceInit = nanmean(grids.timeToChoiceInit, 3);
end