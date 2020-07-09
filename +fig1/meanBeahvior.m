function meanBeahvior(behBlks, axHand)
%% This function plots the data panels for figure one of the ms

%INPUTS(default values)
%specificPanels('fullFigure')---A string which can be used to specifiy specific panels to be plotted if you don't want to plot the whole figure.

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
for i = 1:numMice
    nBlk = spatialAnalysis.getBlockType(behBlks.blks(i), 'norm');
    grds = prc.getGridsFromBlock(nBlk);
    fracRightTurns(:,:,i) = arrayfun(@(x, y) max([nan grds.fracRightTurns(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid).*condNoNan;
    threshMoveTime(:,:,i) = arrayfun(@(x,y) max([nan grds.timeToThreshMove(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid)*1e3.*condNoNan;
%     threshMoveTime(:,:,i) = threshMoveTime(:,:,i)./nanmean(threshMoveTime(2,:,i));
    firstMoveTime(:,:,i) = arrayfun(@(x,y) max([nan grds.timeToFirstMove(grds.visValues==x & grds.audValues==y)]),visGrid,audGrid)*1e3.*condNoNan;
end

%%
close all
axesHandle = plt.tightSubplot(1,3,2,0.1,[0.15 0.15],[0.1 0.1]);
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 1200 400]);
xlabel('Contrast (%)');
ylabel('Rection time (ms)');
hold on
ylim([250 330])
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);

% rT = bsxfun(@rdivide, threshMoveTime, mean(threshMoveTime(2,5,:),1));
meanData = nanmean(threshMoveTime,3);
stdData = nanstd(threshMoveTime,[],3);
numMicePerCond = sum(~isnan(threshMoveTime),3);

lowBound = meanData-stdData./sqrt(numMicePerCond);
upBound = meanData+stdData./sqrt(numMicePerCond);

plotData = cat(3, meanData, lowBound, upBound);
plt.rowsOfGrid(visGrid(1,:)*100, plotData, plt.selectRedBlueColors(audGrid(:,1)));

%%%%%%%
axesHandle = plt.tightSubplot(1,3,3,0.1,[0.15 0.15],[0.1 0.1]);
xlabel('Contrast (%)');
ylabel('Rection time (ms)');
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
%%%%%%%%


axesHandle = plt.tightSubplot(1,3,1,0.1,[0.15 0.15],[0.1 0.1]);
xlabel('Contrast (%)');
ylabel('Frac. Rightward Choices');
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

end