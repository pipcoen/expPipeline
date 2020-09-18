function modelInactivationM2
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

load('D:\Dropbox (Neuropixels)\MouseData\TempNew');
prmLabels = {'Bias'; 'visScaleIpsi'; 'visScaleConta'; 'N'; 'audScaleIpsi'; 'audScaleContra'};
for i = 1:6
    axesHandle = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
    contData = cellfun(@(x) nanmean(x(1:normEstRepeats,i)),deltaParams);
    sortedData = arrayfun(@(x,y) sort(abs([x{1}(normEstRepeats+1:end,i);y]),'descend'), deltaParams,contData, 'uni', 0);
    
    scanPlot.title = prmLabels{i};
    scanPlot.data = contData;
    scanPlot.pVals = cell2mat(arrayfun(@(x,y) max([find(x==y{1},1) nan])./nShuffles, abs(scanPlot.data), sortedData,'uni', 0));
    scanPlot.data = scanPlot.data;
    disp(scanPlot.data);
    scanPlot.colorBarLimits = [-1 1];
    hist(deltaParams{3,1}(normEstRepeats+1:end,i),100)
    box off;
    hold on
    plot([contData(3,1) contData(3,1)], ylim, 'r', 'linewidth', 2)
    pVal = 10^(ceil(log10(scanPlot.pVals(3,1))));
    if pVal > 0.01
        sigLev = 'NS';
    else
        sigLev = ['p<' num2str(pVal)];
    end
    title([prmLabels{i} ': ' sigLev]);
end
export_fig('D:\Dropbox (Personal)\TalksAndApps\Papers\Coen_2020\FigureParts\3_modelInactivationM2_V2', '-pdf', '-painters');
end