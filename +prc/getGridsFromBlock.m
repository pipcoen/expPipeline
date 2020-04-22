function grids = getGridsFromBlock(blk)
grids = prc.makeGrid(blk);
grids.numTrials = prc.makeGrid(blk, blk.tri.outcome.responseMade==1, @length, 1);
grids.numRightTurns = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @sum, 1);
grids.fracRightTurns = prc.makeGrid(blk, blk.tri.outcome.responseMade==2, @mean, 1);

[~,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), grids.numRightTurns, grids.numTrials, 'uni', 0);
grids.fracRightTurnsLowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
grids.fracRightTurnsHighBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));
end