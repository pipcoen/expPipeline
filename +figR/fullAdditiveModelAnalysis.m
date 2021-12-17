function fullAdditiveModelAnalysis
models2fit = {'audOnly_Cross5'; 'visOnly_Cross5'; 'simpLogSplitVSplitAAudDom_Cross5'; 'simpLogSplitVSplitAAudExtraDom_Cross5';  ...
    'simpLogSplitVSplitA_Cross5'; 'simpLogSplitVSplitASplitT_Cross5'; 'SimpEmp_Cross5'; 'fullEmp_Cross5'};
plotOpt = plt.compareModels(models2fit); close;

tVals = {'audOnly'; 'visOnly'; 'AudDom'; 'exAudDom'; 'AddP'; 'SenBias'; 'AddNP'; 'Full'};
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

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('simpLogSplitVSplitA', [],'log', 1); cla
pow2use = glmBlk.glmFit{1}.prmFits(4);

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('simpLogSplitVSplitAAudDom', [],'log', 1,pow2use)
glmBlk.viewGLMFits('simpLogSplitVSplitA', [],'log', 1, pow2use); 

axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('simpLogSplitVSplitAAudExtraDom', [],'log', 1,pow2use)

axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('simpLogSplitVSplitASplitT', [],'log', 1,pow2use)

axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('simpEmp', [],'log', 1,pow2use)

axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('visOnly', [],'log', 1,pow2use)

axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg); cla
glmBlk.viewGLMFits('audOnly', [],'log', 1,pow2use)

export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\FullAddModelFitAnalysis_GLMs', '-pdf', '-painters');
%%
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

logLikDiff = cellfun(@(x) x-plotOpt.yData{1}, plotOpt.yData([6 2 3 4 5 7 8]), 'uni', 0);
tIdx = [2 3 4 5 7 8];
for i = 1:6   
    plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg); cla;
    plot(logLikDiff{1}, logLikDiff{i+1}, '.k'); box off;
    hold on;
    ylim([0 0.5])
    xlim([0 0.5])
    plot(xlim, xlim, '--k')
    title(tVals(tIdx(i)-1))
end
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\FullAddModelFitAnalysis_Scatter', '-pdf', '-painters');

%%
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

plt.tightSubplot(nRows,nCols,1:2,axesGap,botTopMarg,lftRgtMarg); cla;
logLikDiff = cellfun(@(x) x-plotOpt.yData{end}, plotOpt.yData([2 3 4 5 6 7 8 9]), 'uni', 0);
ylim([-0.5 0.0])
nXPnts = length(logLikDiff);
yDat = cell2mat(arrayfun(@(x) [logLikDiff{x}; mean(logLikDiff{x})], 1:nXPnts, 'uni', 0));
xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));

set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
hold on
for i = 1:nXPnts-1
    cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(1:end-1,i:i+1),2), num2cell(yDat(1:end-1,i:i+1),2));
    cellfun(@(x,y) plot(x,y, 'b','HandleVisibility','off'), num2cell(xDat(end,i:i+1),2), num2cell(yDat(end,i:i+1),2));
end
set(gca, 'XTick', xDat(1,:), 'XTickLabel', tVals)
for i = 1:nXPnts
   [~, p] = ttest(yDat(1:end-1,5), yDat(1:end-1,i))
end
xlim([0 nXPnts]);
ylabel('LLR Relative to Optimality')
set(gca, 'position', get(gca, 'position').*[1.3 1 1 1]);

%%

plt.tightSubplot(nRows,nCols,4:5,axesGap,botTopMarg,lftRgtMarg); cla;
logLikDiff = cellfun(@(x) x-plotOpt.yData{1}, plotOpt.yData([2 3 4 5 6 7 8 9]), 'uni', 0);
ylim([0 0.6])
nXPnts = length(logLikDiff);
yDat = cell2mat(arrayfun(@(x) [logLikDiff{x}; mean(logLikDiff{x})], 1:nXPnts, 'uni', 0));
xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));

set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
hold on
for i = 1:nXPnts-1
    cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(1:end-1,i:i+1),2), num2cell(yDat(1:end-1,i:i+1),2));
    cellfun(@(x,y) plot(x,y, 'b','HandleVisibility','off'), num2cell(xDat(end,i:i+1),2), num2cell(yDat(end,i:i+1),2));
end
set(gca, 'XTick', xDat(1,:), 'XTickLabel', tVals)
for i = 1:nXPnts
   [~, p] = ttest(yDat(1:end-1,3), yDat(1:end-1,i))
end
xlim([0 nXPnts]);
ylabel('LLR Relative to bias')
set(gca, 'position', get(gca, 'position').*[1.3 1 1 1]);


%%
plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg); cla;
opt.faceColors = repmat({[0 0 0]}, length(logLikDiff{1}),1);
opt.linkedGroups = [1 2];
plt.jitter(logLikDiff, opt);
end