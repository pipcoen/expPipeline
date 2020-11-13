function modelFitAnalysis
%% This function plots the data panels for figure one of the ms
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
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 100 figWidth, figHeight]);

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
plotOpt = plt.compareModels({'audOnly_Cross5';'visOnly_Cross5';'simpLogSplitVSplitA_Cross5';'fullEmp_Cross5'}); close;

scatter(axesHandle, plotOpt.yData{4}-plotOpt.yData{1}, plotOpt.yData{5}-plotOpt.yData{1}, 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
axis square; 
hold on
xlim([0.2 0.5])
ylim([0.2 0.5])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off

axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
logLikDiff = cellfun(@(x) x-plotOpt.yData{5}, plotOpt.yData(1:4), 'uni', 0);
plt.jitter(logLikDiff);

export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\2_modelFitAnalysis', '-pdf', '-painters');
end