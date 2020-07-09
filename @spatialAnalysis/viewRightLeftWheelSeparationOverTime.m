function viewRightLeftWheelSeparationOverTime(obj)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices
%	'rea'--------------------timeToFirstMove vs fration of rightward choices

% xDat = 0:50:500;
xDat = 0:0.005:0.25;
moveDatAVC = cell(4,1);
for i  = 1:length(obj.blks)
    bigBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    bigBlk = prc.filtBlock(bigBlk, ~isnan(bigBlk.tri.outcome.timeToFirstMove) & bigBlk.tri.outcome.timeToFirstMove<max(xDat));
%     contast2Use = max(abs(bigBlk.tri.stim.visContrast(bigBlk.tri.trialClass.visual)));
%     bigBlk = prc.filtBlock(bigBlk, (abs(bigBlk.tri.stim.visContrast) == contast2Use) | bigBlk.tri.trialClass.auditory);
    
    bigBlk.tri.wheelPos = cellfun(@(x) interp1(x(:,1), x(:,2), xDat, 'nearest', 'extrap'), bigBlk.tri.raw.wheelTimeValue, 'uni', 0);
    bigBlk.tri.right = (bigBlk.tri.stim.visInitialAzimuth>0 & bigBlk.tri.trialClass.visual) | bigBlk.tri.stim.audInitialAzimuth > 0;
    bigBlk.tri.left = (bigBlk.tri.stim.visInitialAzimuth<0 & bigBlk.tri.trialClass.visual) | bigBlk.tri.stim.audInitialAzimuth < 0;
    filtIdx = {bigBlk.tri.trialClass.auditory; bigBlk.tri.trialClass.visual; bigBlk.tri.trialClass.conflict};
    
    for j = 1:3
        blk = prc.filtBlock(bigBlk, filtIdx{j});
%         wheelPos = cell2mat(blk.tri.wheelPos);
        wheelPos = cell2mat(arrayfun(@(x,y) [x{1}(xDat<=y) nan*x{1}(xDat>y)], blk.tri.wheelPos, blk.tri.outcome.threshMoveTime, 'uni', 0));
        moveDatAVC{j}(i,:) = abs(nanmean(wheelPos(blk.tri.right, :),1) - nanmean(wheelPos(blk.tri.left, :),1));
        if j == 1; normVal = moveDatAVC{j}(i,end); end
        moveDatAVC{j}(i,:) = moveDatAVC{j}(i,:)./normVal;
        
        shuffIdx = randperm(blk.tot.trials);
        moveDatAVC{j+3}(i,:) = abs(nanmean(wheelPos(blk.tri.right(shuffIdx), :),1) - nanmean(wheelPos(blk.tri.left(shuffIdx), :),1));
        moveDatAVC{j+3}(i,:) = moveDatAVC{j+3}(i,:)./normVal;
    end
end
%%
cla
testDat = moveDatAVC(1:length(moveDatAVC)/2);
controlDat = moveDatAVC(length(moveDatAVC)/2+1:length(moveDatAVC));
sigTest = zeros(length(testDat),length(xDat));
for i = 1:length(testDat)
    for j = 1:length(xDat)
        sigTest(i,j) = ttest(testDat{i}(:,j), controlDat{i}(:,j), 'Alpha',0.01);
    end
end

meanDataAVC = cell2mat(cellfun(@(x) mean(x), testDat, 'uni', 0));
stdDataAVC = cell2mat(cellfun(@(x) std(x), testDat, 'uni', 0));
lowBoundAVC = meanDataAVC-stdDataAVC./sqrt(length(obj.blks));
upBoundAVC = meanDataAVC+stdDataAVC./sqrt(length(obj.blks));

plotData = cat(3, meanDataAVC, lowBoundAVC, upBoundAVC);
plotOpt.Marker = 'none';
colors = [1 0 1; 0 1 1; 0 0 0];
plt.rowsOfGrid(xDat, plotData(:,:,:), colors, plotOpt);
for i = 1:length(testDat)
    plot(xDat(sigTest(i,:)>0), sigTest(i,sigTest(i,:)>0)+0.1*i, '.', 'color', colors(i,:))
end

end