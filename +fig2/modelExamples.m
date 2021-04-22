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
glmBlk = spatialAnalysis('all', 'behavior', 1, 0);
glmBlk.blks = spatialAnalysis.getBlockType(glmBlk.blks, 'norm');
glmBlk.blks = prc.filtBlock(glmBlk.blks, glmBlk.blks.tri.stim.visContrast~=0.06);
paramsAV = glmBlk.blks.exp.conditionParametersAV{find(cellfun(@length, glmBlk.blks.exp.conditionParametersAV)==25,1)};
[vGrid, aGrid] = meshgrid(unique(paramsAV(:,2)),unique(paramsAV(:,1)));

resp = glmBlk.blks.tri.outcome.responseCalc;
vDiff = glmBlk.blks.tri.stim.visDiff;
aDiff = glmBlk.blks.tri.stim.audDiff;
fracRTurnsOrig = arrayfun(@(x,y) mean(resp(vDiff == x & aDiff==y)==2), vGrid, aGrid);
%%
fRightTurnsRep = fracRTurnsOrig*nan;
nRepeats = 10;
for i = 1:nRepeats
    tBlk = glmBlk.copy;
    tBlk.blks = prc.filtBlock(glmBlk.blks, prc.makeFreqUniform(glmBlk.blks.tri.subjectRef));
    tBlk.viewGLMFits('simpLogSplitVSplitA', [], 'none');
    glmRep{i,1} = tBlk.glmFit{1};
    
    resp = tBlk.blks.tri.outcome.responseCalc;
    vDiff = tBlk.blks.tri.stim.visDiff;
    aDiff = tBlk.blks.tri.stim.audDiff;    
    fRightTurnsRep(:,:,i) = arrayfun(@(x,y) mean(resp(vDiff == x & aDiff==y)==2), vGrid, aGrid);
end
%%
glm2Plot = glmRep{1};
glm2Plot.prmFits = mean(cell2mat(cellfun(@(x) x.prmFits,glmRep, 'uni', 0)));
glm2Plot.prmInit = mean(cell2mat(cellfun(@(x) x.prmInit,glmRep, 'uni', 0)));

pHatCalculated = glm2Plot.calculatepHat(glm2Plot.prmFits,'eval');
[grids.visValues, grids.audValues] = meshgrid(unique(glm2Plot.evalPoints(:,1)),unique(glm2Plot.evalPoints(:,2)));
[~, gridIdx] = ismember(glm2Plot.evalPoints, [grids.visValues(:), grids.audValues(:)], 'rows');
plotData = grids.visValues;
plotData(gridIdx) = pHatCalculated(:,2);

for i = [3 6]
    axesHandle = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
    cla;
    fracRightTurns = nanmean(fRightTurnsRep,3);

    contrastPower = 1;
    ylim([0 1])
    xlim([-1 1])
    midPoint = 0.5;
    
    if i == 6
        plotData = log10(plotData./(1-plotData));
        contrastPower = glm2Plot.prmFits(4);
        fracRightTurns = log10(fracRightTurns./(1-fracRightTurns));
        ylim([-2.6 2.6])
        midPoint = 0;
    end
    
    plotOpt.lineStyle = '-';
    plotOpt.Marker = 'none';
    
    visValues = (abs(grids.visValues(1,:))).^contrastPower.*sign(grids.visValues(1,:));
    lineColors = plt.selectRedBlueColors(grids.audValues(:,1));
    plt.rowsOfGrid(visValues, plotData, lineColors, plotOpt);
    
    plotOpt.lineStyle = 'none';
    plotOpt.Marker = '.';
    
    visValues = vGrid(1,:);
    maxContrast = max(abs(visValues(1,:)));
    visValues = abs(visValues./maxContrast).^contrastPower.*sign(visValues);
    plt.rowsOfGrid(visValues, fracRightTurns, lineColors, plotOpt);
    
    box off;
    
    tickVals = -0.80:0.1:0.8;
    tickMarks = sign(tickVals).*(abs(tickVals).^contrastPower)./(maxContrast.^contrastPower);
    tickVals = num2cell(tickVals*100);
    tickVals(2:2:end) = deal({[]});
    set(gca, 'xTick', tickMarks, 'xTickLabel', tickVals);
    
    xL = xlim; hold on; plot(xL,[midPoint midPoint], '--k', 'linewidth', 1.5);
    yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
    axis square;
end
%%
export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\2_modelExamples_equalSubsample', '-pdf', '-painters');
end