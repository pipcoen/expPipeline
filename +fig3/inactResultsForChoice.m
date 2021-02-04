function inactResultsForChoice
% load('fig3aInactResultsForChoice', 'inactResultsForChoice')
load('fig3aInactResultsForChoice', 'inactResultsForChoice')
figure;
subsets = inactResultsForChoice.subsets;
nMice = size(inactResultsForChoice.laserOnData,2);
axesOpt.totalNumOfAxes = length(subsets)*nMice;
axesOpt.btlrMargins =  [10 30 10 10];
axesOpt.gapBetweenAxes = [10 0];
axesOpt.axesSize = [200 200];
axesOpt.numOfRows = 1;
axesOpt.numOfCols = 4;
axesOpt.totalNumOfAxes = 4;


for i = 1:length(subsets)
    for j = 1:nMice
    %%
    nShuffles = size(inactResultsForChoice.shuffLaserOffData{i,j},3);
    contData = inactResultsForChoice.meanContEffects{i,j};
    shuffleData = double(inactResultsForChoice.shuffLaserOnData{i,j} - inactResultsForChoice.shuffLaserOffData{i,j});
    sortedData = cellfun(@squeeze, num2cell(sort(abs(cat(3,shuffleData, contData)),3,'descend'),3), 'uni', 0);
    scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(contData), sortedData,'uni', 0));
    scanPlot.data = contData; disp(contData);
    scanPlot.addTrialNumber = 0;
    sigLevels = (10.^(-2:-1:-10))';
    lastSigLevel = find(sigLevels>min(scanPlot.pVals(:)),1,'last');
    scanPlot.sigLevels = sigLevels(max([1 lastSigLevel-2]):lastSigLevel);
    
    %%Plot the data in the "scanPlot" structure.
    plt.getAxes(axesOpt, ((i-1)*nMice)+j);
    scanPlot.title = subsets{i};
    scanPlot.gridXY = inactResultsForChoice.gridXY{1};
    scanPlot.colorBarLimits = [-0.6 0.6];
    plt.scanningBrainEffects(scanPlot);
    end
end

%"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
%see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
% export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\3_inactResultsForChoice', '-pdf', '-painters');
end