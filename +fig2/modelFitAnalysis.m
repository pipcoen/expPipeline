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
models2fit = {'audOnly_Cross5';'visOnly_Cross5';'simpLogSplitVSplitA_Cross5'; 'simpLogSplitVSplitAUnisensory_Cross5'; 'fullEmp_Cross5'; 'SimpEmp_Cross5'};
plotOpt = plt.compareModels(models2fit); close;
eIdx = double(contains(plotOpt.subjects, {'PC022'}));
eIdx(contains(plotOpt.subjects, {'PC051'})) = 2;
eIdx(end) = 3;

audModVsBias =  plotOpt.yData{2}-plotOpt.yData{1};
visModVsBias =  plotOpt.yData{3}-plotOpt.yData{1};
addModVBias =  plotOpt.yData{4}-plotOpt.yData{1};
addModUniVBias =  plotOpt.yData{5}-plotOpt.yData{1};
fullEmpVbias =  plotOpt.yData{6}-plotOpt.yData{1};

%%
cla
hold on
scatter(axesHandle, addModVBias(~eIdx), fullEmpVbias(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModVBias(eIdx==1), fullEmpVbias(eIdx==1), 50, 'd', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModVBias(eIdx==2), fullEmpVbias(eIdx==2), 50, 's', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModVBias(eIdx==3), fullEmpVbias(eIdx==3), 25, '^', 'filled', 'MarkerEdgeColor', 'none');
axis square; 
xlim([0.2 0.5])
ylim([0.2 0.5])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);

[~, pVal] = ttest(addModVBias(1:end-1), fullEmpVbias(1:end-1));
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(length(eIdx)-1) '   P<' num2str(pVal)]);
box off

%%
axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg);
cla
hold on
scatter(axesHandle, addModVBias(~eIdx), fullEmpVbias(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModVBias(eIdx==1), fullEmpVbias(eIdx==1), 50, 'dk', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModVBias(eIdx==2), fullEmpVbias(eIdx==2), 50, 'sk', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModVBias(eIdx==3), fullEmpVbias(eIdx==3), 25, '^k', 'filled', 'MarkerEdgeColor', 'none');

scatter(axesHandle, audModVsBias(~eIdx), fullEmpVbias(~eIdx), 25, 'r', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, audModVsBias(eIdx==1), fullEmpVbias(eIdx==1), 50, 'dr', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, audModVsBias(eIdx==2), fullEmpVbias(eIdx==2), 50, 'sr', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, audModVsBias(eIdx==3), fullEmpVbias(eIdx==3), 25, '^r', 'filled', 'MarkerEdgeColor', 'none');

scatter(axesHandle, visModVsBias(~eIdx), fullEmpVbias(~eIdx), 25, 'b', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, visModVsBias(eIdx==1), fullEmpVbias(eIdx==1), 50, 'db', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, visModVsBias(eIdx==2), fullEmpVbias(eIdx==2), 50, 'sb', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, visModVsBias(eIdx==3), fullEmpVbias(eIdx==3), 25, '^b', 'filled', 'MarkerEdgeColor', 'none');

axis square; 
xlim([0 0.5])
ylim([0 0.5])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);

[~, pVal] = ttest(addModVBias(1:end-1), fullEmpVbias(1:end-1));
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(length(eIdx)-1) '   P<' num2str(pVal)]);
box off

%%
axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
cla;
logLikDiff = cellfun(@(x) x-plotOpt.yData{6}, plotOpt.yData([1:3 5]), 'uni', 0);
opt.faceColors = repmat([0 0 0], length(addModVBias),1);
opt.faceColors(eIdx>0,:) = [1 0 0; 0 0 1; 0 1 0];
opt.faceColors = repmat({opt.faceColors}, length(logLikDiff),1);
ylim([-0.6 0.05])
plt.jitter(logLikDiff, opt);

%%
axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);
cla;
logLikDiff = cellfun(@(x) x(1:end-1)-plotOpt.yData{6}(1:end-1), plotOpt.yData([1:3 5]), 'uni', 0);
opt.faceColors = repmat([0 0 0], length(addModVBias)-1,1);
opt.faceColors(eIdx(1:end-1)>0,:) = [1 0 0; 0 0 1];
opt.faceColors = repmat({opt.faceColors}, length(logLikDiff),1);
opt.pairs2test = 'all';
xlim([0 5]);
plt.jitter(logLikDiff, opt);

% %%
% axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg);
% cla
% hold on
% scatter(axesHandle, addModVBias(~eIdx), simpEmpVbias(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
% scatter(axesHandle, addModVBias(eIdx==1), simpEmpVbias(eIdx==1), 50, 'd', 'filled', 'MarkerEdgeColor', 'none');
% scatter(axesHandle, addModVBias(eIdx==2), simpEmpVbias(eIdx==2), 50, 's', 'filled', 'MarkerEdgeColor', 'none');
% scatter(axesHandle, addModVBias(eIdx==3), simpEmpVbias(eIdx==3), 25, '^', 'filled', 'MarkerEdgeColor', 'none');
% axis square; 
% xlim([0.2 0.5])
% ylim([0.2 0.5])
% plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
% 
% [~, pVal] = ttest(addModVBias(1:end-1), simpEmpVbias(1:end-1));
% pVal = round(pVal, 4, 'significant');
% title(['n=' num2str(length(eIdx)-1) '   P<' num2str(pVal)]);
% box off

%%
axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg);
cla;
logLikDiff = cellfun(@(x) x-plotOpt.yData{end}, plotOpt.yData([4 7 6]), 'uni', 0);
opt.faceColors = repmat([0 0 0], length(addModVBias),1);
opt.faceColors(eIdx>0,:) = [1 0 0; 0 0 1; 0 1 0];
opt.faceColors = repmat({opt.faceColors}, length(logLikDiff),1);
ylim([-0.03 0.01])
plt.jitter(logLikDiff, opt);


%%
axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg);
cla
hold on
scatter(axesHandle, addModUniVBias(~eIdx), fullEmpVbias(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModUniVBias(eIdx==1), fullEmpVbias(eIdx==1), 50, 'd', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModUniVBias(eIdx==2), fullEmpVbias(eIdx==2), 50, 's', 'filled', 'MarkerEdgeColor', 'none');
scatter(axesHandle, addModUniVBias(eIdx==3), fullEmpVbias(eIdx==3), 25, '^', 'filled', 'MarkerEdgeColor', 'none');
axis square; 
xlim([0.2 0.5])
ylim([0.2 0.5])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);

[~, pVal] = ttest(addModUniVBias(1:end-1), fullEmpVbias(1:end-1));
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(length(eIdx)-1) '   P<' num2str(pVal)]);
box off
%%
% export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\ModelFitAnalysis', '-pdf', '-painters');
end