function modelInactivation
%% This function plots the data panels for figure one of the ms
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
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 100 figWidth, figHeight]);

load('D:\Dropbox (Neuropixels)\MouseData\TempFin');
for i = 1:6
    contData = cellfun(@(x) nanmean(x(1:normEstRepeats,i)),deltaParams);
    sortedData = arrayfun(@(x,y) sort(abs([x{1}(normEstRepeats+1:end,i);y]),'descend'), deltaParams,contData, 'uni', 0);
    
    obj.hand.axes = plt.getAxes(axesOpt, i);
    scanPlot.title = prmLabels{i};
    scanPlot.data = contData;
    scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(scanPlot.data), sortedData,'uni', 0));
    scanPlot.data = scanPlot.data;
    disp(scanPlot.data);
    scanPlot.colorBarLimits = [-1 1];
    sigLevels = (10.^(-2:-1:-10))';
    lastSigLevel = find(sigLevels>min(scanPlot.pVals(:)),1,'last');
    scanPlot.sigLevels = sigLevels(max([1 lastSigLevel-2]):lastSigLevel);
    scanPlot.sigLevels = [0.01 0.001 0.0001]';
    plt.scanningBrainEffects(scanPlot);
    axesHandle = plt.tightSubplot(nRows,nCols,i+size(galvoRef,1),axesGap,botTopMarg,lftRgtMarg);


end


%%
export_fig('D:\Dropbox (Personal)\TalksAndApps\Papers\Coen_2020\3_visAndM2Inactivation', '-pdf', '-painters');
end