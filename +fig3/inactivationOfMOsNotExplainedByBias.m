function inactivationOfMOsNotExplainedByBias
% if ~exist('uniBlks', 'var'); uniBlks = spatialAnalysis('all', 'uniscan', 1, 1); end
uniMice = {'PC027'; 'PC029'; 'DJ008'; 'DJ006'; 'DJ007'};
for i = 1:length(uniMice)
    uniBlks = spatialAnalysis(uniMice{i}, 'uniscan',1,1);
    mOpt = struct;
    mOpt.contOnly = 1;
    mOpt.nRepeats = 10;
    mOpt.useDif = 0;
    cGLMs = uniBlks.getModelFitsToInactivationData(mOpt);
    
    
    mOpt.contParams = mean(cell2mat(cellfun(@(x) x.prmFits, cGLMs, 'uni', 0)));
    mOpt.useDif = 1;
    mOpt.useGroups = 1;
    mOpt.contOnly = 0;
    mOpt.freeP = [1,1,1,0,1,1]>0;
    
    mOpt.groupIDs = 'mos';
    [~, mGLMsFull{i,1}] = uniBlks.getModelFitsToInactivationData(mOpt);
    mOpt.freeP = [1,0,0,0,0,0]>0;
    [~, mGLMsBias{i,1}] = uniBlks.getModelFitsToInactivationData(mOpt);
end

logLikF = cell2mat(cellfun(@(x) mean(cell2mat(cellfun(@(y) mean(y.logLik),x{1}, 'uni', 0))), mGLMsFull, 'uni', 0));
logLikB = cell2mat(cellfun(@(x) mean(cell2mat(cellfun(@(y) mean(y.logLik),x{1}, 'uni', 0))), mGLMsBias, 'uni', 0));

%%
uniBlks = spatialAnalysis('all', 'uniscan',1,1);
mOpt = struct;
mOpt.contOnly = 1;
mOpt.nRepeats = 10;
mOpt.useDif = 0;
cGLMs = uniBlks.getModelFitsToInactivationData(mOpt);


mOpt.contParams = mean(cell2mat(cellfun(@(x) x.prmFits, cGLMs, 'uni', 0)));
mOpt.useDif = 1;
mOpt.useGroups = 1;
mOpt.contOnly = 0;
mOpt.freeP = [1,1,1,0,1,1]>0;
mOpt.crossVal = 0;
mOpt.crossValFolds = 2;

mOpt.groupIDs = 'mos';
[~, mGLMsFull] = uniBlks.getModelFitsToInactivationData(mOpt);
gridFull = cellfun(@(x) prc.getGridsFromBlock(x.blockData), mGLMsFull{1}, 'uni', 0);
mOpt.freeP = [1,0,0,0,0,0]>0;
[~, mGLMsBias] = uniBlks.getModelFitsToInactivationData(mOpt);
gridBias = cellfun(@(x) prc.getGridsFromBlock(x.blockData), mGLMsFull{1}, 'uni', 0);
%%
figure
axHeight = 250;
axWidth = 250;
nCols = 3;
nRows = 2;
figHeight = nRows*axHeight;
figWidth = nCols*axWidth;
plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);

for i = 1:2
    if i == 1 
        glmPlot = mGLMsFull{1}{1};
        glmPlot.prmFits = mean(cell2mat(cellfun(@(x) x.prmFits, mGLMsFull{1}, 'uni', 0)));
        gridDat = gridFull;
    else
        glmPlot = mGLMsBias{1}{1};
        glmPlot.prmFits = mean(cell2mat(cellfun(@(x) x.prmFits, mGLMsFull{1}, 'uni', 0)));
        gridDat = gridBias;
    end
    
    fracRightTurns = squeeze(nanmean(cell2mat(cellfun(@(x) permute(x.fracRightTurns,[3,1,2]), gridDat, 'uni', 0)),1));
    
    pHatCalculated = glmPlot.calculatepHat(glmPlot.prmFits,'eval');
    [grids.visValues, grids.audValues] = meshgrid(unique(glmPlot.evalPoints(:,1)),unique(glmPlot.evalPoints(:,2)));
    [~, gridIdx] = ismember(glmPlot.evalPoints, [grids.visValues(:), grids.audValues(:)], 'rows');
    plotData = grids.visValues;
    plotData(gridIdx) = pHatCalculated(:,2);
    
    contrastPower = 1;
    ylim([0 1])
    xlim([-1 1])
    midPoint = 0.5;
    
    plotData = log10(plotData./(1-plotData));
    contrastPower = mOpt.contParams(4);
    fracRightTurns = log10(fracRightTurns./(1-fracRightTurns));
    ylim([-3 3])
    xlim([-1 1])
    midPoint = 0;
    
    if i == 2
        plotOpt.lineStyle = '--';
        plotOpt.Marker = 'none';
    else
        plotOpt.lineStyle = '-';
        plotOpt.Marker = 'none';
    end
    
    visValues = (abs(grids.visValues(1,:))).^contrastPower.*sign(grids.visValues(1,:));
    lineColors = plt.selectRedBlueColors(grids.audValues(:,1));
    plt.rowsOfGrid(visValues, plotData, lineColors, plotOpt);
    
    plotOpt.lineStyle = 'none';
    plotOpt.Marker = '.';
    
    visValues = gridDat{1}.visValues(1,:);
    maxContrast = max(abs(visValues(1,:)));
    visValues = abs(visValues./maxContrast).^contrastPower.*sign(visValues);
    plt.rowsOfGrid(visValues, fracRightTurns, lineColors, plotOpt);
    
    box off;
    set(gca, 'xTick', visValues, 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)*100));
    
    xL = xlim; hold on; plot(xL,[midPoint midPoint], '--k', 'linewidth', 1.5);
    yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
end
%%
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\MOsInactivationFullVsBias', '-pdf', '-painters');

end