function viewRightLeftWheelSeparationOverTimeAllRL(obj)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices
%	'rea'--------------------timeToFirstMove vs fration of rightward choices

% xDat = 0:50:500;
xDat = 0:0.005:0.25;
moveDatAV = cell(2,1);
for i  = 1:length(obj.blks)
    bigBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    bigBlk = prc.filtBlock(bigBlk, ~isnan(bigBlk.tri.outcome.timeToFirstMove) & bigBlk.tri.outcome.timeToFirstMove<max(xDat));
    bigBlk = prc.filtBlock(bigBlk, ~bigBlk.tri.trialClass.blank);
    bigBlk.tri.wheelPos = cellfun(@(x) interp1(x(:,1), x(:,2), xDat, 'nearest', 'extrap'), bigBlk.tri.raw.wheelTimeValue, 'uni', 0);
    filtIdx = {~bigBlk.tri.trialClass.visual; ~bigBlk.tri.trialClass.auditory; bigBlk.tri.trialClass.conflict | bigBlk.tri.trialClass.coherent};
    
    
    for j = 1:2
        blk = prc.filtBlock(bigBlk, filtIdx{j});
        if j == 1 || j==3; rightIdx = blk.tri.stim.audInitialAzimuth>0; else, rightIdx = blk.tri.stim.visInitialAzimuth>0;  end
%         wheelPos = cell2mat(blk.tri.wheelPos);
        wheelPos = cell2mat(arrayfun(@(x,y) [x{1}(xDat<=y) nan*x{1}(xDat>y)], blk.tri.wheelPos, blk.tri.outcome.threshMoveTime, 'uni', 0));
        moveDatAV{j}(i,:) = abs(nanmean(wheelPos(rightIdx, :),1) - nanmean(wheelPos(~rightIdx, :),1));
        if j == 1; normVal = moveDatAV{j}(i,end); end
        moveDatAV{j}(i,:) = moveDatAV{j}(i,:)./normVal;
    end
end
%%
cla
meanDataAVC = cell2mat(cellfun(@(x) mean(x), moveDatAV, 'uni', 0));
stdDataAVC = cell2mat(cellfun(@(x) std(x), moveDatAV, 'uni', 0));
lowBoundAVC = meanDataAVC-stdDataAVC./sqrt(length(obj.blks));
upBoundAVC = meanDataAVC+stdDataAVC./sqrt(length(obj.blks));

plotData = cat(3, meanDataAVC, lowBoundAVC, upBoundAVC);
plotOpt.Marker = 'none';
plt.rowsOfGrid(xDat, plotData, [1 0 1; 0 1 1; 0 0 0], plotOpt);
end