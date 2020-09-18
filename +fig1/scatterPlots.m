function scatterPlots
%%
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1, 'raw'); end
[perAud, perHighVis, perLowVis, perHighMul, perLowMul, perMaxVis] = deal(nan*ones(length(behBlks.blks), 1));
[threshTimeAV, threshTimeCohCon] = deal(nan*ones(length(behBlks.blks), 2));
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
    
  
    threshTimeAV(i,1) = mean(grds.timeToThreshMove(grds.visValues == 0 & grds.audValues ~= 0));
    threshTimeAV(i,2) = mean(grds.timeToThreshMove(abs(grds.visValues) == closestVis & grds.audValues == 0));
    timeToFirstMove{i,1} = nBlk.tri.outcome.timeToFirstMove;
    timeToThreshMove{i,1} = nBlk.tri.outcome.threshMoveTime;
    
    cohIdx = grds.visValues.*grds.audValues > 0 &~isnan(grds.timeToThreshMove);
    conIdx = grds.visValues.*grds.audValues < 0 &~isnan(grds.timeToThreshMove) & flipud(cohIdx);
    threshTimeCohCon(i,:) = [mean(grds.timeToThreshMove(cohIdx)) mean(grds.timeToThreshMove(conIdx))];
        
end

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
scatter(axesHandle, perAud, perMaxVis, 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
xlim([0.7 1]);
ylim([0.7 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
scatter(axesHandle, perHighMul, perHighVis, 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
xlim([0.7 1]);
ylim([0.7 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);
histogram(cell2mat(timeToThreshMove),100, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 1)
xlim([0 1.5]);
axis square; 
box off;

axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg);
scatter(axesHandle, threshTimeAV(:,1), threshTimeAV(:,2), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
axis square; 
hold on
xlim([0.25 0.45])
ylim([0.25 0.45])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;


axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg);
scatter(axesHandle, threshTimeCohCon(:,1), threshTimeCohCon(:,2), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
axis square; 
hold on
xlim([0.25 0.45])
ylim([0.25 0.45])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;

axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg);
behBlks.viewRightLeftWheelSeparationOverTime;
axis square; 
xlim([0 0.25]);

export_fig('D:\Dropbox (Personal)\TalksAndApps\Papers\Coen_2020\1_scatterPlots', '-pdf', '-painters');
end