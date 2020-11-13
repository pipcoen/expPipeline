function viewGLMFits(obj, modelString, cvFolds, plotType, useCurrentAxes)
if ~exist('modelString', 'var') || isempty(modelString); modelString = 'simpLogSplitVSplitA'; end
if ~exist('cvFolds', 'var') || isempty(cvFolds); cvFolds = 0; end
if ~exist('plotType', 'var') || isempty(plotType); plotType = 'normal'; end
if ~exist('useCurrentAxes', 'var'); useCurrentAxes = 0; end
if ~strcmpi(plotType, 'none'); noPlt = 0; else; noPlt = 1; end
if ~strcmpi(plotType, 'only'); onlyPlt = 0; else; onlyPlt = 1; end

if ~useCurrentAxes && ~strcmpi(plotType, 'none'); figure; end
axesOpt.totalNumOfAxes = length(obj.blks);
axesOpt.btlrMargins = [80 100 80 40];
axesOpt.gapBetweenAxes = [100 60];
axesOpt.numOfRows = min(length(obj.blks), 4);
axesOpt.figureHWRatio = 1.1;
if ~onlyPlt; obj.glmFit = cell(length(obj.blks),1); end
for i  = 1:length(obj.blks)

    
    if ~onlyPlt
        normBlock = spatialAnalysis.getBlockType(obj.blks(i),'norm');
        normBlock = prc.filtBlock(normBlock,~isinf(normBlock.tri.stim.audInitialAzimuth));
        disp(normBlock.tot.trials);
        obj.glmFit{i} = fit.GLMmulti(normBlock, modelString);
    else
        normBlock = obj.blks(i);
    end
    if useCurrentAxes; obj.hand.axes = gca; elseif ~noPlt; obj.hand.axes = plt.getAxes(axesOpt, i); end
    if ~onlyPlt
        if ~cvFolds; obj.glmFit{i}.fit; end
        if cvFolds; obj.glmFit{i}.fitCV(cvFolds); end
    end
    if noPlt; return; end
    
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
    
    visDiff = normBlock.tri.stim.visDiff;
    audDiff = normBlock.tri.stim.audDiff;
    responseMade = normBlock.tri.outcome.responseMade;
    [visGrid, audGrid] = meshgrid(unique(visDiff),unique(audDiff));
    maxContrast = max(abs(visGrid(1,:)));
    fracRightTurns = arrayfun(@(x,y) mean(responseMade(ismember([visDiff,audDiff],[x,y],'rows'))==2), visGrid, audGrid);
    
    visValues = abs(visGrid(1,:)).^contrastPower.*sign(visGrid(1,:))./maxContrast;
    if strcmp(plotType, 'log')
        fracRightTurns = log(fracRightTurns./(1-fracRightTurns));
    end
    plt.rowsOfGrid(visValues, fracRightTurns, lineColors, plotOpt);
    
    xlim([-1 1])
    midPoint = 0.5;
    if strcmp(plotType, 'log')
        maxContrast = ((maxContrast*100).^contrastPower)/100;
        ylim([-6 6])
        midPoint = 0;
    end
    
    box off;
    set(gca, 'xTick', (-1):(1/4):1, 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)*100));
    if ~useCurrentAxes; title(obj.blks(i).exp.subject{1}); end
    xL = xlim; hold on; plot(xL,[midPoint midPoint], '--k', 'linewidth', 1.5);
    yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
end
if ~useCurrentAxes
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
end
obj.hand.figure = [];
end