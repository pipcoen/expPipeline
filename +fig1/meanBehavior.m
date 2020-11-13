function meanBehavior(behBlks, axHand)
%% This function plots the data panels for figure one of the ms
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1); end

%INPUTS(default values)
%specificPanels('fullFigure')---A string which can be used to specifiy specific panels to be plotted if you don't want to plot the whole figure.
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


%%
mLim = 0.15;
for i = 1:length(behBlks.blks); if strcmp(behBlks.blks(i).exp.subject{1}, 'PC022'); exampleIdx = i; end; end
nBlk = spatialAnalysis.getBlockType(behBlks.blks(exampleIdx), 'norm');

nBlk = prc.filtBlock(nBlk, nBlk.tri.outcome.timeDirChoiceInit(:,1) < 0.5);
moveDiff = nBlk.tri.outcome.timeDirChoiceInit(:,1)-nBlk.tri.outcome.timeDirFirstMove(:,1);
nBlk = prc.filtBlock(nBlk, moveDiff == 0);

grds = prc.getGridsFromBlock(nBlk);

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
cla; hold on
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
meanData = grds.fracRightTurns;
lowBound = grds.fracRightTurnsLowBound;
upBound = grds.fracRightTurnsHighBound;

plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(grds.visValues(1,:)*100, plotData, plt.selectRedBlueColors(grds.audValues(:,1)));
axis square

%%
axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
cla;
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
meanData = grds.timeToChoiceInit;
lowBound = meanData + grds.timeChoiceCrossSE;
upBound = meanData - grds.timeChoiceCrossSE;
plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(grds.visValues(1,:)*100, plotData*1000, plt.selectRedBlueColors(grds.audValues(:,1)));
axis square

%%
numMice = length(behBlks.blks);
allParametersAV = cell2mat(arrayfun(@(x) x.exp.conditionParametersAV{1}, behBlks.blks, 'uni', 0));
[uniParametersAV, ~, rowIdx] = unique(allParametersAV, 'rows');
condFreq = histcounts(rowIdx,1:max(rowIdx)+1)';
cond2Use = uniParametersAV(condFreq>=15,:);

%%
[visGrid, audGrid] = meshgrid(unique(cond2Use(:,2)),unique(cond2Use(:,1)));
condNoNan = double(arrayfun(@(x,y) ismember([x,y], cond2Use, 'rows'),audGrid,visGrid));
condNoNan(condNoNan==0) = nan;
[fracRightTurns, numTrials, timeToChoiceInit, firstMoveTime] = deal(nan*ones([size(visGrid), numMice]));

for i = 1:numMice
    nBlk = spatialAnalysis.getBlockType(behBlks.blks(i), 'norm');
    nBlk = prc.filtBlock(nBlk, nBlk.tri.outcome.timeDirChoiceInit(:,1) < 0.5);

    grds = prc.getGridsFromBlock(nBlk);
    fracRightTurns(:,:,i) = arrayfun(@(x, y) max([nan grds.fracRightTurns(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid).*condNoNan;
    timeToChoiceInit(:,:,i) = arrayfun(@(x,y) max([nan grds.timeToChoiceInit(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid)*1e3.*condNoNan;
    firstMoveTime(:,:,i) = arrayfun(@(x,y) max([nan grds.timeToFirstMove(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid)*1e3.*condNoNan;
    numTrials(:,:,i) = arrayfun(@(x,y) max([nan grds.numTrials(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid).*condNoNan;
end
disp(nansum(numTrials(:)));
%%
axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg);
hold on
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);

meanData = nanmean(fracRightTurns,3);
stdData = nanstd(fracRightTurns,[],3);
numMicePerCond = sum(~isnan(fracRightTurns),3);

lowBound = meanData-stdData./sqrt(numMicePerCond);
upBound = meanData+stdData./sqrt(numMicePerCond);

plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(visGrid(1,:)*100, plotData, plt.selectRedBlueColors(audGrid(:,1)));
axis square
%%
axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg);
cla;
hold on
ylim([170 230])
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);

meanData = nanmean(timeToChoiceInit,3);
stdData = nanstd(timeToChoiceInit,[],3);
numMicePerCond = sum(~isnan(timeToChoiceInit),3);

lowBound = meanData-stdData./sqrt(numMicePerCond);
upBound = meanData+stdData./sqrt(numMicePerCond);

plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(visGrid(1,:)*100, plotData, plt.selectRedBlueColors(audGrid(:,1)));
axis square


%%
axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg);
hold on
ylim([120 180])
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);

meanData = nanmean(firstMoveTime,3);
stdData = nanstd(firstMoveTime,[],3);
numMicePerCond = sum(~isnan(firstMoveTime),3);

lowBound = meanData-stdData./sqrt(numMicePerCond);
upBound = meanData+stdData./sqrt(numMicePerCond);

plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(visGrid(1,:)*100, plotData, plt.selectRedBlueColors(audGrid(:,1)));
axis square
%%%%%%%%
% export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\1_meanBehavior', '-pdf', '-painters');
end