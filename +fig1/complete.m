function meanBeahvior(specificPanels)
%% This function plots the data panels for figure one of the ms

%INPUTS(default values)
%specificPanels('fullFigure')---A string which can be used to specifiy specific panels to be plotted if you don't want to plot the whole figure.

if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behavior', 0, 1); end
figure;
xlabel('Contrast (%)');
ylabel('Frac. Rightward Choices');
hold on

numMice = length(behBlks.blks);
allParametersAV = cell2mat(arrayfun(@(x) x.exp.conditionParametersAV{1}, behBlks.blks, 'uni', 0));
[uniParametersAV, ~, rowIdx] = unique(allParametersAV, 'rows');
condFreq = histcounts(rowIdx,1:max(rowIdx)+1)';
cond2Use = uniParametersAV(condFreq>=15,:);

[visValues, audValues] = meshgrid(unique(cond2Use(:,2)),unique(cond2Use(:,1)));
fracRightTurns = nan*ones([size(visValues), numMice]);
for i = 1:numMice
    nBlk = spatialAnalysis.getBlockType(behBlks.blks(i), 'norm');
    grds = prc.getGridsFromBlock(nBlk);
    fracRightTurns(:,:,i) = arrayfun(@(x, y) max([nan grds.fracRightTurns(grds.visValues==x & grds.audValues==y)]),visValues,audValues);
end

meanData = nanmean(fracRightTurns,3);
stdData = nanstd(fracRightTurns,[],3);
numMicePerCond = sum(~isnan(fracRightTurns),3);

plotData.fracRightTurns = meanData;
plotData.visValues = visValues;
plotData.audValues = audValues;
plotData.fracRightTurnsLowBound = meanData-stdData./sqrt(numMicePerCond);
plotData.fracRightTurnsHighBound = meanData+stdData./sqrt(numMicePerCond);

plotOpt.lineStyle = '-';
plotOpt.errorBars = 1;
plotOpt.errorBars = 1;

plt.fracitonRightChoiceData(plotData, plotOpt);
end