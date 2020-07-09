function viewConflictChoiceVsRT(obj)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices
%	'rea'--------------------timeToFirstMove vs fration of rightward choices

% xDat = 0:50:500;
xDat = 0:20:250;
[visDat, audDat, conDat] = deal(nan*ones(length(obj.blks), length(xDat)-1));
% [visDat, audDat, conDat] = deal(nan*ones(length(obj.blks), length(xDat)));
for i  = 1:length(obj.blks)
    bigBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    bigBlk = prc.filtBlock(bigBlk, ~isnan(bigBlk.tri.outcome.timeToFirstMove) & bigBlk.tri.outcome.timeToFirstMove<max(xDat));   
%     bigBlk = prc.filtBlock(bigBlk, abs(bigBlk.tri.stim.visContrast) == max(abs(bigBlk.tri.stim.visContrast)) | bigBlk.tri.trialClass.auditory);

%     bigBlk = prc.filtBlock(bigBlk, bigBlk.tri.outcome.firstMoveReliable==1);
    minTrials = 10;
    perRel(i,1) = mean(bigBlk.tri.outcome.firstMoveReliable(bigBlk.tri.outcome.timeToFirstMove<0.1));
    
    blk = prc.filtBlock(bigBlk, bigBlk.tri.trialClass.visual);
    moveDir = blk.tri.outcome.firstMoveDirection;
    time2Move = blk.tri.outcome.timeToFirstMove;
    visChosen = (moveDir == 2) == (blk.tri.stim.visInitialAzimuth>0);
    [binCounts, ~, binIdx] = histcounts(time2Move*1000, xDat);
    visChosenForEachRT = arrayfun(@(x) mean(visChosen(binIdx==x)), 1:length(xDat)-1);
    visDat(i,binCounts>minTrials) = visChosenForEachRT(binCounts>minTrials);
    visReliable(i,1) = mean(blk.tri.outcome.firstMoveReliable(time2Move<0.1));
    
    
    
    blk = prc.filtBlock(bigBlk, bigBlk.tri.trialClass.auditory);
    moveDir = blk.tri.outcome.firstMoveDirection;
    time2Move = blk.tri.outcome.timeToFirstMove;
    audChosen = (moveDir == 2) == (blk.tri.stim.audInitialAzimuth>0);
    [binCounts, ~, binIdx] = histcounts(time2Move*1000, xDat);
    audChosenForEachRT = arrayfun(@(x) mean(audChosen(binIdx==x)), 1:length(xDat)-1);
    audDat(i,binCounts>minTrials) = audChosenForEachRT(binCounts>minTrials);   
    audReliable(i,1) = mean(blk.tri.outcome.firstMoveReliable(time2Move<0.1));
    
    blk = prc.filtBlock(bigBlk, bigBlk.tri.trialClass.conflict);
    moveDir = blk.tri.outcome.firstMoveDirection;
    time2Move = blk.tri.outcome.timeToFirstMove;
    audChosen = (moveDir == 2) == (blk.tri.stim.audInitialAzimuth>0);
    [binCounts, ~, binIdx] = histcounts(time2Move*1000, xDat);
    audChosenForEachRT = arrayfun(@(x) mean(audChosen(binIdx==x)), 1:length(xDat)-1);
    conDat(i,binCounts>minTrials) = audChosenForEachRT(binCounts>minTrials);
    conReliable(i,1) = mean(blk.tri.outcome.firstMoveReliable(time2Move<0.1));

end
xTicks = xDat(2:end)-(diff(xDat(1:2)))/2;

% audDat = audDat(perRel>0.50, :);
% visDat = visDat(perRel>0.50, :);
% conDat = conDat(perRel>0.50, :);

audDat = abs(audDat-0.5);
visDat = abs(visDat-0.5);
conDat = abs(conDat-0.5);

audPlot = zeros(size(audDat,2),3);
audPlot(:,2) = nanmean(audDat,1);
audPlot(:,1) = nanmean(audDat,1)-(nanstd(audDat,1));
audPlot(:,3) = nanmean(audDat,1)+(nanstd(audDat,1));

visPlot = zeros(size(visDat,2),3);
visPlot(:,2) = nanmean(visDat,1);
visPlot(:,1) = nanmean(visDat,1)-(nanstd(visDat,1));
visPlot(:,3) = nanmean(visDat,1)+(nanstd(visDat,1));

conPlot = zeros(size(conDat,2),3);
conPlot(:,2) = nanmean(conDat,1);
conPlot(:,1) = nanmean(conDat,1)-(nanstd(conDat,1));
conPlot(:,3) = nanmean(conDat,1)+(nanstd(conDat,1));

figureSize = get(gcf, 'position');
plt.lineWithPatch(xTicks, audPlot, 'c');
plt.lineWithPatch(xTicks, visPlot, 'm');
plt.lineWithPatch(xTicks, conPlot, 'k');
xlabel('Reaction time (ms)');
ylabel('Fraction of "correct" initial movments');


    
%     blk = prc.filtBlock(bigBlk, bigBlk.tri.trialClass.visual);
%     blk = prc.filtBlock(blk, abs(blk.tri.stim.visContrast) == max(abs(blk.tri.stim.visContrast)));
%     wheelPosVis = cellfun(@(x) interp1(x(:,1), x(:,2), xDat, 'nearest', 'extrap'), blk.tri.raw.wheelTimeValue, 'uni', 0);
%     wheelPosVis = cell2mat(arrayfun(@(x,y) [x{1}(xDat<=y) nan*x{1}(xDat>y)], wheelPosVis, blk.tri.outcome.threshMoveTime, 'uni', 0));
%     wheelPosVis = bsxfun(@times, wheelPosVis, ((blk.tri.outcome.threshMoveDirection-2)*2)+1)*-1;
%     wheelPosVis = bsxfun(@rdivide, wheelPosVis, nanmax(wheelPosVis));
%     visDat(i,:) = nanmean(wheelPosVis);
% 
%     blk = prc.filtBlock(bigBlk, bigBlk.tri.trialClass.auditory);
%     wheelPosAud = cellfun(@(x) interp1(x(:,1), x(:,2), xDat, 'nearest', 'extrap'), blk.tri.raw.wheelTimeValue, 'uni', 0);
%     wheelPosAud = cell2mat(arrayfun(@(x,y) [x{1}(xDat<=y) nan*x{1}(xDat>y)], wheelPosAud, blk.tri.outcome.threshMoveTime, 'uni', 0));
%     wheelPosAud = bsxfun(@times, wheelPosAud, ((blk.tri.outcome.threshMoveDirection-2)*2)+1)*-1;
%     wheelPosAud = bsxfun(@rdivide, wheelPosAud, nanmax(wheelPosAud));
%     audDat(i,:) = nanmean(wheelPosAud);
% 
%     blk = prc.filtBlock(bigBlk, bigBlk.tri.trialClass.conflict);
%     wheelPosCon = cellfun(@(x) interp1(x(:,1), x(:,2), xDat, 'nearest', 'extrap'), blk.tri.raw.wheelTimeValue, 'uni', 0);
%     wheelPosCon = cell2mat(arrayfun(@(x,y) [x{1}(xDat<=y) nan*x{1}(xDat>y)], wheelPosCon, blk.tri.outcome.threshMoveTime, 'uni', 0));
%     wheelPosCon = bsxfun(@times, wheelPosCon, sign(blk.tri.stim.audInitialAzimuth)*-1);
%     wheelPosCon = bsxfun(@rdivide, wheelPosCon, nanmax(wheelPosCon));
%     conDat(i,:) = nanmean(wheelPosCon);
    

% plot(xDat(2:end)-(diff(xDat(1:2))), nanmean(yDat));
end