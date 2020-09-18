function modelExamples
%% This function plots the data panels for figure one of the ms
mice2Plot = {'PC011'; 'PC022';'DJ010'};

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
    axesHandle = plt.tightSubplot(nRows,nCols,i,axesGap,botTopMarg,lftRgtMarg);
    behBlks.viewGLMFits('simpLogSplitVSplitA', [], [], 1)
    axis square;

    axesHandle = plt.tightSubplot(nRows,nCols,i+3,axesGap,botTopMarg,lftRgtMarg);
    behBlks.viewGLMFits('simpLogSplitVSplitA', [],'log', 1)
    axis square;
    
end

export_fig('D:\Dropbox (Personal)\TalksAndApps\Papers\Coen_2020\2_modelExamples', '-pdf', '-painters');
end