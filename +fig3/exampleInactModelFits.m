function exampleInactModelFits(uniBlks)

if ~exist('uniBlks', 'var'); uniBlks = spatialAnalysis('all', 'uniscan', 1, 1); end
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
%%
mOpt.groupIDs = 'mos';
[~, mGLMs] = uniBlks.getModelFitsToInactivationData(mOpt);
mOpt.groupIDs = 'v1';
[~, vGLMs] = uniBlks.getModelFitsToInactivationData(mOpt);
mOpt.groupIDs = 's1';
[~, sGLMs] = uniBlks.getModelFitsToInactivationData(mOpt);
mOpt.groupIDs = 'a1';
[~, aGLMs] = uniBlks.getModelFitsToInactivationData(mOpt);

glmsCont_V1_MOs_A1_S1 = [{cGLMs}; vGLMs; mGLMs; aGLMs; sGLMs];
%%
gridsCont_V1_MOs_A1_S1 = cell(length(glmsCont_V1_MOs_A1_S1),1);
for i = 1:length(gridsCont_V1_MOs_A1_S1)
    gridsCont_V1_MOs_A1_S1{i,1} = cellfun(@(x) prc.getGridsFromBlock(x.blockData,1), glmsCont_V1_MOs_A1_S1{i}, 'uni', 0);
end
%%
axHeight = 250;
axWidth = 250;
nCols = 3;
nRows = 2;
figHeight = nRows*axHeight;
figWidth = nCols*axWidth;

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;

titles = {'Control'; 'V1'; 'MOs'; 'A1'; 'S1'};
for j = 2:2
    figure;
    set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);
    
    for i = 1:5
        plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg); cla;
        
        allGlMs = glmsCont_V1_MOs_A1_S1{i};
        glm2Plot = glmsCont_V1_MOs_A1_S1{i}{1};
        glm2Plot.prmFits = mean(cell2mat(cellfun(@(x) x.prmFits,allGlMs, 'uni', 0)));
        glm2Plot.prmInit = mean(cell2mat(cellfun(@(x) x.prmInit,allGlMs, 'uni', 0)));
        fracRightTurns = squeeze(mean(cell2mat(cellfun(@(x) permute(x.fracRightTurnsComb,[3,1,2]), gridsCont_V1_MOs_A1_S1{i}, 'uni', 0)),1));
        numTrials = round(squeeze(mean(cell2mat(cellfun(@(x) permute(x.numTrials,[3,1,2]), gridsCont_V1_MOs_A1_S1{i}, 'uni', 0)),1)));
        
        if any(fracRightTurns(:) == 1)
            fracRightTurns(fracRightTurns==1) = numTrials(fracRightTurns==1)/(numTrials(fracRightTurns==1)+1);
        end
        
        pHatCalculated = glm2Plot.calculatepHat(glm2Plot.prmFits,'eval');
        [grids.visValues, grids.audValues] = meshgrid(unique(glm2Plot.evalPoints(:,1)),unique(glm2Plot.evalPoints(:,2)));
        [~, gridIdx] = ismember(glm2Plot.evalPoints, [grids.visValues(:), grids.audValues(:)], 'rows');
        plotData = grids.visValues;
        plotData(gridIdx) = pHatCalculated(:,2);
        
        contrastPower = 1;
        ylim([0 1])
        xlim([-1 1])
        midPoint = 0.5;
                
        if j == 2
            plotData = log10(plotData./(1-plotData));
            contrastPower = mOpt.contParams(4);
            fracRightTurns = log10(fracRightTurns./(1-fracRightTurns));
            ylim([-2.5 2.5])
            midPoint = 0;
        end        
        
        plotOpt.lineStyle = '-';
        plotOpt.Marker = 'none';
        
        visValues = (abs(grids.visValues(1,:))).^contrastPower.*sign(grids.visValues(1,:));
        lineColors = plt.selectRedBlueColors(grids.audValues(:,1));
        plt.rowsOfGrid(visValues, plotData, lineColors, plotOpt);
        
        plotOpt.lineStyle = 'none';
        plotOpt.Marker = '.';
        
        visValues = gridsCont_V1_MOs_A1_S1{i}{1}.visValues(1,:);
        maxContrast = max(abs(visValues(1,:)));
        visValues = abs(visValues./maxContrast).^contrastPower.*sign(visValues);
        plt.rowsOfGrid(visValues, fracRightTurns, lineColors, plotOpt);
        
        box off;
        
        tickVals = -0.80:0.1:0.8;
        tickMarks = sign(tickVals).*(abs(tickVals).^contrastPower)./(maxContrast.^contrastPower);
        tickVals = num2cell(tickVals*100);
        tickVals(2:2:end) = deal({[]});
        set(gca, 'xTick', tickMarks, 'xTickLabel', tickVals);
        
        axis square;        
        xL = xlim; hold on; plot(xL,[midPoint midPoint], '--k', 'linewidth', 1.5);
        yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
        title(titles{i});
    end
    
    %"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
    %see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
    if j == 1
%         export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\3_exampleInactModelFits', '-pdf', '-painters');
    elseif j == 2
        export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\3efg_exampleInactModelFits', '-pdf', '-painters');
    end
end