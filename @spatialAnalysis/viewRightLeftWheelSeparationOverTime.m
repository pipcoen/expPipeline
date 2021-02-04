function viewRightLeftWheelSeparationOverTime(obj, contasts2Use)
%% A method for the spatialAnalysis class to plot data without any fit for a all the blocks.
% INPUTS(default values)
% plotType('res')--------String indicating the type of data to plot. Options are
%	'res'--------------------contrast vs fration of rightward choices

% xDat = 0:50:500;
xDat = 0.01:0.005:0.3;
moveDatAVC = cell(3,1);
if ~exist('contasts2Use', 'var'); contasts2Use = 0; end
for i  = 1:length(obj.blks)
    bigBlk = spatialAnalysis.getBlockType(obj.blks(i),'norm',1);
    
    if contasts2Use
        bigBlk = prc.filtBlock(bigBlk, abs(bigBlk.tri.stim.visContrast) == contasts2Use | bigBlk.tri.trialType.auditory);
    end
    
    bigBlk.tri.wheelPos = cellfun(@(x) interp1(x(:,1), x(:,2), xDat, 'nearest', 'extrap'), bigBlk.tri.raw.wheelTimeValue, 'uni', 0);
    bigBlk.tri.right = (bigBlk.tri.stim.visInitialAzimuth>0 & bigBlk.tri.trialType.visual) | bigBlk.tri.stim.audInitialAzimuth > 0;
    bigBlk.tri.left = (bigBlk.tri.stim.visInitialAzimuth<0 & bigBlk.tri.trialType.visual) | bigBlk.tri.stim.audInitialAzimuth < 0;
    filtIdx = {bigBlk.tri.trialType.auditory; bigBlk.tri.trialType.visual; bigBlk.tri.trialType.conflict};
    
    if ~any(bigBlk.tri.trialType.visual); continue; end
    for shuffRepeat = 1:100
        for j = 1:3
            blk = prc.filtBlock(bigBlk, filtIdx{j});
            blk = prc.filtBlock(blk, prc.makeFreqUniform((blk.tri.right)*2+(blk.tri.left)));
            wheelPos = cell2mat(blk.tri.wheelPos);
            wheelPos = bsxfun(@rdivide, bsxfun(@minus, wheelPos, nanmean(wheelPos)), nanstd(wheelPos)+0.01);
            controlDat = abs(nanmean(wheelPos(blk.tri.right, :),1) - nanmean(wheelPos(blk.tri.left, :),1));
            
            shuffIdx = randperm(blk.tot.trials);
            shuffDat = abs(nanmean(wheelPos(blk.tri.right(shuffIdx), :),1) - nanmean(wheelPos(blk.tri.left(shuffIdx), :),1));
            moveDatAVC{j}(i,:,shuffRepeat) = double(controlDat - shuffDat);
        end
    end
end
%%
cla
testDat = cellfun(@(x) nanmean(x,3), moveDatAVC, 'uni', 0);
sigTest = zeros(length(testDat),length(xDat));
for i = 1:length(testDat)
    for j = 1:length(xDat)
        [~, sigTest(i,j)] = ttest(testDat{i}(:,j));
        sigTest(i,j) = sigTest(i,j)<0.01;
    end
end
meanDataAVC = cell2mat(cellfun(@(x) nanmean(x), testDat, 'uni', 0));
stdDataAVC = cell2mat(cellfun(@(x) nanstd(x), testDat, 'uni', 0));
lowBoundAVC = meanDataAVC-stdDataAVC./sqrt(length(obj.blks));
upBoundAVC = meanDataAVC+stdDataAVC./sqrt(length(obj.blks));

plotData = cat(3, meanDataAVC, lowBoundAVC, upBoundAVC);
plotOpt.Marker = 'none';
colors = [1 0 1; 0 1 1; 0 0 0];
plt.rowsOfGrid(xDat, plotData(:,:,:), colors, plotOpt);
for i = 1:length(testDat)
    plot(xDat(sigTest(i,:)>0), sigTest(i,sigTest(i,:)>0)-0.05*i, '.', 'color', colors(i,:))
end

end