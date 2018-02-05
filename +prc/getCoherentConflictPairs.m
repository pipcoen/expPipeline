function paired = getCoherentConflictPairs(blk)
%% Function to combine and filter block files.
%%
trialTypes = prc.makeGrid(blk, blk.trialType, @mean, 'abscondition');
pairedTrials = ((trialTypes.*flip(trialTypes,2))>0) & trialTypes == 3;
blk.timeToWheelMove(isnan(blk.timeToWheelMove)) = nanmean(blk.timeToWheelMove);
reactionTimes = nanmedian(prc.makeGrid(blk, blk.timeToWheelMove, @median, 'abscondition', 1),3);
paired.yData = [{reactionTimes(pairedTrials.*trialTypes==3)*1000}; {reactionTimes(flip((pairedTrials.*trialTypes==3),2))*1000}];
paired.xTickLocations = [0.5, 1.5];
paired.linkedGroups = {[1 2]};
paired.pairs2test = paired.linkedGroups;