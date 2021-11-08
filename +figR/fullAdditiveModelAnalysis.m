function fullAdditiveModelAnalysis
%% This function plots the data panels for figure one of the ms
figure;
gcaeight = 250;
axWidth = 250;
nCols = 3;
nRows = 2;
figHeight = nRows*gcaeight;
figWidth = nCols*axWidth;

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 200 figWidth, figHeight]);

%%
models2fit = {'simpLogSplitVSplitA_Cross5'; 'SimpEmp_Cross5'; 'fullEmp_Cross5'};
plotOpt = plt.compareModels(models2fit); close;
eIdx = double(contains(plotOpt.subjects, {'PC022'}));
eIdx(contains(plotOpt.subjects, {'PC051'})) = 2;
eIdx(end) = 3;
%%
plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); cla;
logLikDiff = cellfun(@(x) x-plotOpt.yData{end}, plotOpt.yData(2:3), 'uni', 0);
ylim([-0.015 0.0])
nXPnts = 2;
yDat = cell2mat(arrayfun(@(x) [logLikDiff{x}(~eIdx); logLikDiff{x}(eIdx>0); mean(logLikDiff{x})], 1:nXPnts, 'uni', 0));
xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));

set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
hold on
for i = 1:nXPnts-1
    cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(:,i:i+1),2), num2cell(yDat(:,i:i+1),2));
end

plot(gca, xDat(end-1,:), yDat(end-1,:),'sc', 'MarkerEdgeColor', 'c','MarkerFaceColor', 'c', 'MarkerSize',6);
plot(gca, xDat(end-2,:), yDat(end-2,:),'sc', 'MarkerEdgeColor', 'c','MarkerFaceColor', 'c', 'MarkerSize',6);
plot(gca, xDat(end,:), yDat(end,:),'^c', 'MarkerEdgeColor', 'c','MarkerFaceColor', 'c', 'MarkerSize',6);
[~, pVal] = ttest(yDat(1:end-1,1), yDat(1:end-1,2))
xlim([0 nXPnts]);
%%
behBlks = spatialAnalysis('all', 'behavior', 0, 1, '');
allParametersAV = cell2mat(arrayfun(@(x) x.exp.conditionParametersAV{1}, behBlks.blks, 'uni', 0));
[uniParametersAV, ~, rowIdx] = unique(allParametersAV, 'rows');
condFreq = histcounts(rowIdx,1:max(rowIdx)+1)';
cond2Use = uniParametersAV(condFreq>=15,:);

glmBlk = spatialAnalysis('all', 'behavior', 1, 0, '');
glmBlk.blks = spatialAnalysis.getBlockType(glmBlk.blks, 'norm');
audDiff = glmBlk.blks.tri.stim.audDiff;
visDiff = glmBlk.blks.tri.stim.visDiff;
glmBlk.blks = prc.filtBlock(glmBlk.blks, ismember([audDiff visDiff], cond2Use, 'rows'));

%%
numExp = glmBlk.blks.tot.experiments;
glmBlk.blks.exp.conditionParametersAV = repmat(glmBlk.blks.exp.conditionParametersAV(end),numExp,1);
axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('simpEmp', [], [], 1)
axis square;
%%
axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('simpEmp', [],'log', 1)
axis square;
glmBlk.viewGLMFits('simpLogSplitVSplitA', [],'log', 1);

%%
axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk = spatialAnalysis('PC022', 'behavior', 1, 0, '');
glmBlk.blks = spatialAnalysis.getBlockType(glmBlk.blks, 'norm');
audDiff = glmBlk.blks.tri.stim.audDiff;
visDiff = glmBlk.blks.tri.stim.visDiff;
glmBlk.blks = prc.filtBlock(glmBlk.blks, ismember([audDiff visDiff], cond2Use, 'rows'));
glmBlk.viewGLMFits('simpEmp', [],'log', 1)
axis square;
glmBlk.viewGLMFits('simpLogSplitVSplitA', [],'log', 1);

%%
axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk = spatialAnalysis('PC051', 'behavior', 1, 0, '');
glmBlk.blks = spatialAnalysis.getBlockType(glmBlk.blks, 'norm');
audDiff = glmBlk.blks.tri.stim.audDiff;
visDiff = glmBlk.blks.tri.stim.visDiff;
glmBlk.blks = prc.filtBlock(glmBlk.blks, ismember([audDiff visDiff], cond2Use, 'rows'));
glmBlk.viewGLMFits('simpEmp', [],'log', 1)
axis square;
glmBlk.viewGLMFits('simpLogSplitVSplitA', [],'log', 1);


%%
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\FullAddModelFitAnalysis', '-pdf', '-painters');
end