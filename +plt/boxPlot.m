function boxPlot(boxPlot, addText)
if ~exist('addText', 'var'); addText = 1; end
if ~isfield(boxPlot, 'plotLabels'); boxPlot.plotLabels = boxPlot.plotData; end
plotData = boxPlot.plotData;
imsc(plotData, boxPlot.axisLimits, boxPlot.colorMap, 'k'); 
daspect([1 1 1]); axis xy;
[xPnts, yPnts] = meshgrid(1:size(plotData,2), 1:size(plotData,1));
if addText
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center'), xPnts, yPnts, plotData)
end
if boxPlot.nSessions == 1; title(sprintf('%s: %d Tri, %s', boxPlot.subject, boxPlot.trialNumber, boxPlot.expDate))
else, title(sprintf('%s: %d Tri, %d Sess', boxPlot.subject, boxPlot.trialNumber, boxPlot.nSessions));
end
set(gca, 'xTick', 1:size(plotData,2), 'xTickLabel', boxPlot.xyValues{1}, 'fontsize', 14)
set(gca, 'yTick', 1:size(plotData,1), 'yTickLabel', boxPlot.xyValues{2}, 'fontsize', 14, 'TickLength', [0, 0])
box off;
end