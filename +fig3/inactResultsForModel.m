function inactResultsForModel
load('fig3bModelInactivationResults', 'modelInactivationResults')
figure;
axesOpt.btlrMargins =  [10 30 10 10];
axesOpt.gapBetweenAxes = [10 0];
axesOpt.axesSize = [200 200];
axesOpt.numOfRows = 1;
axesOpt.numOfCols = 6;
axesOpt.totalNumOfAxes = 6;

normEstRepeats = modelInactivationResults.normEstRepeats;
deltaParams = modelInactivationResults.deltaParams;
contParams = modelInactivationResults.contParams;
prmLabels = ({'Bias'; 'vIpsi'; 'vContra'; 'y'; 'aIpsi'; 'aContra'});
for i = find(modelInactivationResults.freeP)
    nShuffles = size(deltaParams{1},1) - normEstRepeats;
    contParams(cellfun(@isempty, contParams)) = deal({nan*ones(max(max(cellfun(@length, contParams))),6)});
    contData = cellfun(@(x) nanmean(x(1:normEstRepeats,i)),deltaParams);
    stdData = cellfun(@(x) nanstd(x(normEstRepeats+1:end,i)),deltaParams);
    sortedData = arrayfun(@(x,y) sort(abs([x{1}(normEstRepeats+1:end,i);y]),'descend'), deltaParams,contData, 'uni', 0);
    
    plt.getAxes(axesOpt, i);
    scanPlot.title = prmLabels{i};
    scanPlot.data = contData./stdData;
    scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(contData), sortedData,'uni', 0));
    scanPlot.colorBarLimits = [-10 10];
    scanPlot.sigLevels = [0.01; 0.001; 0.0001];
    scanPlot.gridXY = modelInactivationResults.gridXY;
    plt.scanningBrainEffects(scanPlot);
end

%"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
%see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_inactResultsForModel', '-pdf', '-painters');
end