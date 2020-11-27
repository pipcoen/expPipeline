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
eIdx = double(contains(plotOpt.subjects, {'PC022'}));
eIdx(contains(plotOpt.subjects, {'PC051'})) = 2;
eIdx(end) = 3;

modelVbias =  plotOpt.yData{4}-plotOpt.yData{1};
fullEmpVbias =  plotOpt.yData{5}-plotOpt.yData{1};

%%
cla
hold on
scatter(axesHandle, modelVbias(~eIdx), fullEmpVbias(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, modelVbias(eIdx==1), fullEmpVbias(eIdx==1), 50, 'd', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, modelVbias(eIdx==2), fullEmpVbias(eIdx==2), 50, 's', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, modelVbias(eIdx==3), fullEmpVbias(eIdx==3), 25, '^', 'filled', 'MarkerEdgeColor', 'none');
axis square; 
xlim([0.2 0.5])
ylim([0.2 0.5])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);

[~, pVal] = ttest(modelVbias(1:end-1), fullEmpVbias(1:end-1));
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(length(eIdx)-1) '   P<' num2str(pVal)]);
box off
%%
axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
cla;
logLikDiff = cellfun(@(x) x-plotOpt.yData{5}, plotOpt.yData(1:4), 'uni', 0);
opt.faceColors = repmat([0 0 0], length(modelVbias),1);
opt.faceColors(eIdx>0,:) = [1 0 0; 0 0 1; 0 1 0];
opt.faceColors = repmat({opt.faceColors}, length(logLikDiff),1);
plt.jitter(logLikDiff, opt);
%%
axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);
cla;
logLikDiff = cellfun(@(x) x(1:end-1)-plotOpt.yData{5}(1:end-1), plotOpt.yData(1:4), 'uni', 0);
opt.faceColors = repmat([0 0 0], length(modelVbias)-1,1);
opt.faceColors(eIdx(1:end-1)>0,:) = [1 0 0; 0 0 1];
opt.faceColors = repmat({opt.faceColors}, length(logLikDiff),1);
opt.pairs2test = 'all';
xlim([0 5]);
plt.jitter(logLikDiff, opt);
%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\2_modelFitAnalysis', '-pdf', '-painters');
end