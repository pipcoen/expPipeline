function grids = getGridsFromBlock(blk)
%% Function to extract useful "grids" from a block, arranged with aud/vis values in rows/columns

outC = blk.tri.outcome;
responseCorrect = blk.tri.stim.correctResponse==outC.responseCalc;

grids = prc.makeGrid(blk);                                                       %Grids of aud/vis values and the condition labels
grids.performance = prc.makeGrid(blk, responseCorrect==1, @mean,'abscondition'); %Grid of performance
grids.numTrials = prc.makeGrid(blk, outC.responseCalc==1, @length);              %Grid of trial numbers

%Get grid of fraction of right turns for each session, then calculate the SE across sessions. Mean of these sessions is the overall fration
grids.fracRightTurns = prc.makeGrid(blk, outC.responseCalc==2, @nanmean, [], 1);
grids.fracRightTurnsSE = nanstd(grids.fracRightTurns,[],3)./sqrt(size(grids.fracRightTurns,3));
grids.fracRightTurns = nanmean(grids.fracRightTurns, 3);

%Get confidence binomial intervals for fraction of right turns (more appropriate for combined mice)
numRightTurns = prc.makeGrid(blk, blk.tri.outcome.responseCalc==2, @sum);
[~,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), numRightTurns, grids.numTrials, 'uni', 0);
grids.fracRightTurnsLowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0)).*(grids.fracRightTurns*0+1);
grids.fracRightTurnsHighBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0)).*(grids.fracRightTurns*0+1);


%Get grid of meadian reaction times for each session, then calculate the SE across sessions. Mean of these sessions is the overall reaction time
grids.reactionTime = prc.makeGrid(blk, outC.reactionTime, @nanmedian, [], 1);
grids.reactionTimeSE = nanstd(grids.reactionTime,[],3)./sqrt(size(grids.reactionTime,3));
grids.reactionTime = nanmean(grids.reactionTime, 3);
grids.reactionTimeComb = prc.makeGrid(blk, outC.reactionTime, @nanmedian);
grids.timeToResponseThreshComb = prc.makeGrid(blk, outC.timeToResponseThresh, @nanmedian);
end