function meanBehavior(behBlks, axHand)
%% This function plots the data panels for figure one of the ms

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

oneMouse = spatialAnalysis('PC022', 'behavior', 0, 1);
nBlk = spatialAnalysis.getBlockType(oneMouse.blks, 'norm');
grds = prc.getGridsFromBlock(nBlk);

axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
hold on
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
meanData = grds.fracRightTurns;
lowBound = grds.fracRightTurnsLowBound;
upBound = grds.fracRightTurnsHighBound;

plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(grds.visValues(1,:)*100, plotData, plt.selectRedBlueColors(grds.audValues(:,1)));
axis square


axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
meanData = grds.timeToThreshMove;
lowBound = meanData + grds.timeToThreshMoveMAD;
upBound = meanData - grds.timeToThreshMoveMAD;
plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(grds.visValues(1,:)*100, plotData*1000, plt.selectRedBlueColors(grds.audValues(:,1)));
axis square

%%
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1); end
numMice = length(behBlks.blks);
allParametersAV = cell2mat(arrayfun(@(x) x.exp.conditionParametersAV{1}, behBlks.blks, 'uni', 0));
[uniParametersAV, ~, rowIdx] = unique(allParametersAV, 'rows');
condFreq = histcounts(rowIdx,1:max(rowIdx)+1)';
cond2Use = uniParametersAV(condFreq>=15,:);

%%
[visGrid, audGrid] = meshgrid(unique(cond2Use(:,2)),unique(cond2Use(:,1)));
condNoNan = double(arrayfun(@(x,y) ismember([x,y], cond2Use, 'rows'),audGrid,visGrid));
condNoNan(condNoNan==0) = nan;
fracRightTurns = nan*ones([size(visGrid), numMice]);
totTrials = 0;
for i = 1:numMice
    nBlk = spatialAnalysis.getBlockType(behBlks.blks(i), 'norm');
    totTrials = totTrials+nBlk.tot.trials;
    grds = prc.getGridsFromBlock(nBlk);
    fracRightTurns(:,:,i) = arrayfun(@(x, y) max([nan grds.fracRightTurns(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid).*condNoNan;
    threshMoveTime(:,:,i) = arrayfun(@(x,y) max([nan grds.timeToThreshMove(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid)*1e3.*condNoNan;
    firstMoveTime(:,:,i) = arrayfun(@(x,y) max([nan grds.timeToFirstMove(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid)*1e3.*condNoNan;
end
disp(totTrials);
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
hold on
ylim([250 330])
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);

meanData = nanmean(threshMoveTime,3);
stdData = nanstd(threshMoveTime,[],3);
numMicePerCond = sum(~isnan(threshMoveTime),3);

lowBound = meanData-stdData./sqrt(numMicePerCond);
upBound = meanData+stdData./sqrt(numMicePerCond);

plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(visGrid(1,:)*100, plotData, plt.selectRedBlueColors(audGrid(:,1)));
axis square


%%%%%%%
axesHandle = plt.tightSubplot(nRows,nCols,6,axesGap,botTopMarg,lftRgtMarg);
hold on
ylim([150 200])
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
export_fig('D:\Dropbox (Personal)\TalksAndApps\Papers\Coen_2020\1_meanBehavior', '-pdf', '-painters');
end