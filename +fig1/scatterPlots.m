function scatterPlots
%%
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1); end
[perAud, perVis80, perVis40, perMul40] = deal(nan*ones(length(behBlks.blks), 1));
[threshTimeAV, threshTimeCohCon] = deal(nan*ones(length(behBlks.blks), 2));
[timeToThreshMove, timeToFirstMove] = deal(cell(length(behBlks.blks), 1));

vis2Use = 0.4;
exampIdx = zeros(length(behBlks.blks), 1);
for i = 1:length(behBlks.blks)
    if strcmp(behBlks.blks(i).exp.subject{1}, 'PC022'); exampIdx(i)  =1; end
    nBlk = spatialAnalysis.getBlockType(behBlks.blks(i), 'norm');
    grds = prc.getGridsFromBlock(nBlk);
    
    if ~any(grds.visValues(:) == vis2Use); continue; end
    perAud(i) = grds.performance(grds.visValues == 0 & grds.audValues > 0);
    perVis40(i) = grds.performance(grds.visValues == vis2Use & grds.audValues == 0);
    if any(grds.visValues(:) == 0.8)
        perVis80(i) = grds.performance(grds.visValues == 0.8 & grds.audValues == 0);
    end
    perMul40(i) = grds.performance(grds.visValues == vis2Use & grds.audValues > 0);
       
    threshTimeAV(i,1) = mean(grds.timeToThreshMove(grds.visValues == 0 & grds.audValues ~= 0));
    threshTimeAV(i,2) = mean(grds.timeToThreshMove(abs(grds.visValues) == vis2Use & grds.audValues == 0));
    timeToThreshMove{i,1} = nBlk.tri.outcome.threshMoveTime;
    timeToFirstMove{i,1} = nBlk.tri.outcome.timeToFirstMove;
    
    cohIdx = grds.visValues.*grds.audValues > 0 &~isnan(grds.timeToThreshMove) & abs(grds.visValues) == vis2Use;
    conIdx = grds.visValues.*grds.audValues < 0 &~isnan(grds.timeToThreshMove) & flipud(cohIdx) & abs(grds.visValues) == vis2Use;
    threshTimeCohCon(i,:) = [mean(grds.timeToThreshMove(cohIdx)) mean(grds.timeToThreshMove(conIdx))];
end

%%
popIdx = ~any(isnan([perVis40, perAud, perMul40 threshTimeCohCon threshTimeAV]),2) & ~exampIdx;
exampIdx = exampIdx>0;
allIdx = (popIdx + exampIdx)>0;
%%
figure;
axHeight = 250;
axWidth = 250;
nCols = 3;
nRows = 2;
figHeight = nRows*axHeight;
figWidth = nCols*axWidth;

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
hold on
scatter(axesHandle, perAud(popIdx), perVis40(popIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, perAud(exampIdx), perVis40(exampIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, mean(perAud(allIdx)), mean(perVis40(allIdx)), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(perAud, perVis40);
pVal = round(pVal, 2, 'significant');
title(['n=' num2str(sum(allIdx)) '   P<' num2str(pVal)]);
xlim([0.5 1]);
ylim([0.5 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); hold on;
scatter(axesHandle, perMul40(popIdx), perVis40(popIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, perMul40(exampIdx), perVis40(exampIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, mean(perMul40(allIdx)), mean(perVis40(allIdx)), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(perMul40, perVis40);
pVal = round(pVal, 2, 'significant');
title(['n=' num2str(sum(allIdx)) '   P<' num2str(pVal)]);
xlim([0.5 1]);
ylim([0.5 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); hold on;
histogram(cell2mat(timeToThreshMove),100, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 1);
xlim([0 1.5]);
axis square; 
box off;

axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); hold on;
scatter(axesHandle, threshTimeAV(popIdx,1), threshTimeAV(popIdx,2), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, threshTimeAV(exampIdx, 1), threshTimeAV(exampIdx, 2), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, mean(threshTimeAV(allIdx, 1)), mean(threshTimeAV(allIdx, 2)), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(threshTimeAV(allIdx,1), threshTimeAV(allIdx,2));
pVal = round(pVal, 2, 'significant');
title(['n=' num2str(sum(allIdx)) '   P<' num2str(pVal)]);
axis square; 
xlim([0.2 0.4])
ylim([0.2 0.4])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); hold on;
scatter(axesHandle, threshTimeCohCon(popIdx,1), threshTimeCohCon(popIdx,2), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, threshTimeCohCon(exampIdx, 1), threshTimeCohCon(exampIdx, 2), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, mean(threshTimeCohCon(allIdx, 1)), mean(threshTimeCohCon(allIdx, 2)), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(threshTimeCohCon(allIdx,1), threshTimeCohCon(allIdx,2));
pVal = round(pVal, 2, 'significant');
title(['n=' num2str(sum(allIdx)) '   P<' num2str(pVal)]);
axis square; 
hold on
xlim([0.2 0.4])
ylim([0.2 0.4])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg); hold on;
cellfun(@(x,y) plot(x,[y, y], 'k'), num2cell(threshTimeCohCon(popIdx,:),2), num2cell(threshTimeAV(popIdx,2)));
scatter(axesHandle, threshTimeCohCon(popIdx,1),  threshTimeAV(popIdx,2),25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, threshTimeCohCon(popIdx,2),threshTimeAV(popIdx,2), 25, 'w', 'filled', 'MarkerEdgeColor', 'k');

plot(threshTimeCohCon(exampIdx, :), threshTimeAV(exampIdx, 2)*[1 1], 'k')
plot(mean(threshTimeCohCon(allIdx, :)), mean(threshTimeAV(allIdx, 2))*[1 1], 'k')
scatter(axesHandle, threshTimeCohCon(exampIdx, 1), threshTimeAV(exampIdx, 2), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, mean(threshTimeCohCon(allIdx, 1)), mean(threshTimeAV(allIdx, 2)), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, threshTimeCohCon(exampIdx, 2), threshTimeAV(exampIdx, 2), 50, 'w', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axesHandle, mean(threshTimeCohCon(allIdx, 2)), mean(threshTimeAV(allIdx, 2)), 25, 'w', '^', 'filled', 'MarkerEdgeColor', 'k');



[~, pVal] = ttest(threshTimeCohCon(allIdx,1), threshTimeCohCon(allIdx,2));
pVal = round(pVal, 2, 'significant');
title(['n=' num2str(sum(allIdx)) '   P<' num2str(pVal)]);
axis square; 
hold on
xlim([0.2 0.4])
ylim([0.2 0.4])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;

%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\1_scatterPlots', '-pdf', '-painters');
end