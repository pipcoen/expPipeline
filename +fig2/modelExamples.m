function modelExamples
%% This function plots the data panels for figure one of the ms
mice2Plot = {'PC051'; 'PC022'};

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
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 200 figWidth, figHeight]);

for i = 1:length(mice2Plot)
    behBlks = spatialAnalysis(mice2Plot(i), 'behavior', 0, 1);
    behBlks.blks = prc.filtBlock(behBlks.blks, behBlks.blks.tri.stim.visContrast ~= 0.06);
    axesHandle = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
    behBlks.viewGLMFits('simpLogSplitVSplitA', [], [], 1)
    axis square;

    axesHandle = plt.tightSubplot(nRows,nCols,i+3,axesGap,botTopMarg,lftRgtMarg);
    behBlks.viewGLMFits('simpLogSplitVSplitA', [],'log', 1)
    axis square;
end


%%
behBlks = spatialAnalysis('all', 'behavior', 0, 1);
allParametersAV = cell2mat(arrayfun(@(x) x.exp.conditionParametersAV{1}, behBlks.blks, 'uni', 0));
[uniParametersAV, ~, rowIdx] = unique(allParametersAV, 'rows');
condFreq = histcounts(rowIdx,1:max(rowIdx)+1)';
cond2Use = uniParametersAV(condFreq>=15,:);

glmBlk = spatialAnalysis('all', 'behavior', 1, 0);
glmBlk.blks = spatialAnalysis.getBlockType(glmBlk.blks, 'norm');
audDiff = glmBlk.blks.tri.stim.audDiff;
visDiff = glmBlk.blks.tri.stim.visDiff;
glmBlk.blks = prc.filtBlock(glmBlk.blks, ismember([audDiff visDiff], cond2Use, 'rows'));

%%
numExp = glmBlk.blks.tot.experiments;
glmBlk.blks.exp.conditionParametersAV = repmat(glmBlk.blks.exp.conditionParametersAV(end),numExp,1);
axesHandle = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg);
glmBlk.viewGLMFits('simpLogSplitVSplitA', [], [], 1)
axis square;

axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg);
glmBlk.viewGLMFits('simpLogSplitVSplitA', [],'log', 1)
axis square;

%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\2_modelExamples', '-pdf', '-painters');
end