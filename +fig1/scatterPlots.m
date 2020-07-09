function scatterPlots
%%
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1); end

%%
[perAud, perHighVis, perLowVis, perHighMul, perLowMul, perMaxVis, threshTimeAV, conPer] = deal(nan*ones(length(behBlks.blks), 1));
[timeToThreshMove, timeToFirstMove] = deal(cell(length(behBlks.blks), 1));

for i = 1:length(behBlks.blks)
    nBlk = spatialAnalysis.getBlockType(behBlks.blks(i), 'norm');
    grds = prc.getGridsFromBlock(nBlk);
    
    avPairs = nBlk.exp.conditionParametersAV{1};
    visValues = unique(abs(avPairs(~any(avPairs==0,2) & sum(sign(avPairs),2)~=0,2)));
    perAud(i) = grds.performance(grds.visValues == 0 & grds.audValues > 0);
    perLowVis(i) = grds.performance(grds.visValues == min(visValues) & grds.audValues == 0);
    perHighVis(i) = grds.performance(grds.visValues == max(visValues) & grds.audValues == 0);
    perLowMul(i) = grds.performance(grds.visValues == min(visValues) & grds.audValues > 0);
    perHighMul(i) = grds.performance(grds.visValues == max(visValues) & grds.audValues > 0);
    perMaxVis(i) = max(grds.performance(grds.audValues == 0)); 
    
    avDiffPer = grds.performance(grds.visValues > 0 & grds.audValues == 0)-perAud(i);
    selectedVis = grds.visValues(grds.visValues > 0 & grds.audValues == 0);
    closestVis = selectedVis(abs(avDiffPer) == min(abs(avDiffPer)));
    
    if min(abs(avDiffPer))<0.02
        conPer(i) = grds.performance(grds.visValues == closestVis*-1 & grds.audValues >0);
    end
    
    threshTimeAV(i,1) = mean(grds.timeToThreshMove(grds.visValues == 0 & grds.audValues ~= 0));
    threshTimeAV(i,2) = mean(grds.timeToThreshMove(abs(grds.visValues) == closestVis & grds.audValues == 0));
    timeToFirstMove{i,1} = nBlk.tri.outcome.timeToFirstMove;
    timeToThreshMove{i,1} = nBlk.tri.outcome.threshMoveTime;
end

%%
figure;
numRow = 1; numCol = 10;
axesHandle = plt.tightSubplot(numRow,numCol,1:2,0.05,[0.15 0.15],[0.05 0.05]);
set(gcf, 'position', get(gcf, 'position').*[2 1 0 0] + [0 0 numCol*150 300]);

scatter(axesHandle, perAud, perMaxVis, 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
xlim([0.7 1]);
ylim([0.7 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(numRow,numCol,3:4,0.05,[0.15 0.15],[0.05 0.05]);
scatter(axesHandle, max([perAud perHighVis], [], 2), perHighMul, 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
xlim([0.7 1]);
ylim([0.7 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(numRow,numCol,5:6,0.05,[0.15 0.15],[0.05 0.05]);
histogram(cell2mat(timeToThreshMove),100, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 1)
xlim([0 1.5]);
axis square; 
box off;

axesHandle = plt.tightSubplot(numRow,numCol,7:8,0.05,[0.15 0.15],[0.05 0.05]);
scatter(axesHandle, threshTimeAV(:,1), threshTimeAV(:,2), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
axis square; 
hold on
xlim([0.25 0.45])
ylim([0.25 0.45])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;
end