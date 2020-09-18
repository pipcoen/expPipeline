function visAndM2InactivationBarCharts
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

s.viewInactivationResults('readifgrp', 25000)


export_fig('D:\Dropbox (Personal)\TalksAndApps\Papers\Coen_2020\FigureParts\3_modelInactivationAud_V2', '-pdf', '-painters');
end