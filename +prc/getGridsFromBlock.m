function grids = getGridsFromBlock(blk)
grids = prc.makeGrid(blk);

grids.conditions = prc.makeGrid(blk, blk.tri.stim.conditionLabel, @mean);
grids.performance = prc.makeGrid(blk, blk.tri.outcome.feedbackGiven==1, @mean,'abscondition');
grids.numTrials = prc.makeGrid(blk, blk.tri.outcome.responseMade==1, @length);
grids.numRightTurns = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @sum);
grids.fracRightTurns = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @mean);

[~,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), grids.numRightTurns, grids.numTrials, 'uni', 0);
grids.fracRightTurnsLowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
grids.fracRightTurnsHighBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));
grids.fracRightTurnsLowBound(isnan(grids.fracRightTurns)) = nan;
grids.fracRightTurnsHighBound(isnan(grids.fracRightTurns)) = nan;


meadianMad = @(x) mad(x,1);

grids.timeToFirstMove = prc.makeGrid(blk, blk.tri.outcome.timeToFirstMove, @nanmedian);
grids.timeToFirstMoveMAD = prc.makeGrid(blk, blk.tri.outcome.timeToFirstMove, meadianMad);
grids.timeToThreshMove = prc.makeGrid(blk, blk.tri.outcome.threshMoveTime, @nanmedian);
grids.timeToThreshMoveMAD = prc.makeGrid(blk, blk.tri.outcome.threshMoveTime, meadianMad);
end