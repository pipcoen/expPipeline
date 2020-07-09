function viewGLMFits(obj, modelString, cvFolds, plotType)
if ~exist('modelString', 'var'); modelString = 'simpLogSplitVSplitA'; end
if ~exist('cvFolds', 'var'); cvFolds = 0; end
if ~exist('plotType', 'var'); plotType = 'normal'; end

figure;
axesOpt.totalNumOfAxes = length(obj.blks);
axesOpt.btlrMargins = [80 100 80 40];
axesOpt.gapBetweenAxes = [100 60];
axesOpt.numOfRows = min(length(obj.blks), 4);
axesOpt.figureHWRatio = 1.1;
obj.glmFit = cell(length(obj.blks),1);
for i  = 1:length(obj.blks)
    normBlock = spatialAnalysis.getBlockType(obj.blks(i),'norm');
    normBlock = prc.filtBlock(normBlock,~isinf(normBlock.tri.stim.audInitialAzimuth));
    if ~contains(lower(modelString), 'plot'); obj.glmFit{i} = fit.GLMmulti(normBlock, modelString); end
    obj.hand.axes = plt.getAxes(axesOpt, i);
    if ~contains(lower(modelString), 'plot')
        if ~cvFolds; obj.glmFit{i}.fit; end
        if cvFolds; obj.glmFit{i}.fitCV(cvFolds); end
    end
    
    params2use = mean(obj.glmFit{i}.prmFits,1);
    pHatCalculated = obj.glmFit{i}.calculatepHat(params2use,'eval');
    [grids.visValues, grids.audValues] = meshgrid(unique(obj.glmFit{i}.evalPoints(:,1)),unique(obj.glmFit{i}.evalPoints(:,2)));
    [~, gridIdx] = ismember(obj.glmFit{i}.evalPoints, [grids.visValues(:), grids.audValues(:)], 'rows');
    plotData = grids.visValues;
    plotData(gridIdx) = pHatCalculated(:,2);
    plotOpt.lineStyle = '-';
    plotOpt.Marker = 'none';
 
    if strcmp(plotType, 'log')
        contrastPower  = params2use(strcmp(obj.glmFit{i}.prmLabels, 'N'));
        plotData = log(plotData./(1-plotData));
    else
        contrastPower =1;
    end
    visValues = (abs(grids.visValues(1,:))).^contrastPower.*sign(grids.visValues(1,:));
    lineColors = plt.selectRedBlueColors(grids.audValues(:,1));
    plt.rowsOfGrid(visValues, plotData, lineColors, plotOpt);
    
    plotOpt.lineStyle = 'none';
    plotOpt.Marker = '.';
    grids = prc.getGridsFromBlock(normBlock);
    maxContrast = max(abs(grids.visValues(1,:)));
    visValues = abs(grids.visValues(1,:)).^contrastPower.*sign(grids.visValues(1,:))./maxContrast;
    if strcmp(plotType, 'log')
        grids.fracRightTurns = log(grids.fracRightTurns./(1-grids.fracRightTurns));
    end
    plt.rowsOfGrid(visValues, grids.fracRightTurns, lineColors, plotOpt);
    
    xlim([-1 1])
    midPoint = 0.5;
    if strcmp(plotType, 'log')
        maxContrast = ((maxContrast*100).^contrastPower)/100;
        ylim([-6 6])
        midPoint = 0;
    end
    
    box off;
    set(gca, 'xTick', (-1):(1/4):1, 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)*100));
    title(obj.blks(i).exp.subject{1});
    xL = xlim; hold on; plot(xL,[midPoint midPoint], '--k', 'linewidth', 1.5);
    yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
    if length(obj.blks) == 1; set(gcf, 'position', get(gcf, 'position').*[1 0.9 1 1.15]); end
end
figureSize = get(gcf, 'position');
mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];

if strcmp(plotType, 'log')
    plt.suplabel('\fontsize{20} Log(pR/pL)', 'y', mainAxes);
    plt.suplabel('\fontsize{20} Visual Contrast^N', 'x', mainAxes);
    plt.suplabel(['\fontsize{20} ' modelString ': log-axes'], 't', mainAxes);
else
    plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
    plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
    plt.suplabel(['\fontsize{20} ' modelString], 't', mainAxes);
end
obj.hand.figure = [];
end