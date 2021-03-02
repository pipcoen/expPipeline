function boxPlot(boxPlot, addText)
if ~exist('addText', 'var'); addText = 1; end
if ~isfield(boxPlot, 'plotLabels'); boxPlot.plotLabels = boxPlot.plotData; end
if iscell(boxPlot.subject); boxPlot.subject = boxPlot.subject{1}; end
plotData = boxPlot.plotData;

imH = imagesc(plotData);
colormap(boxPlot.colorMap); 
set(imH, 'AlphaData', ~isnan(plotData))
daspect([1 1 1]); axis xy;
set(gca, 'Color', [0, 0, 0])
[xPnts, yPnts] = meshgrid(1:size(plotData,2), 1:size(plotData,1));
if addText
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center'), xPnts, yPnts, plotData)
end
if boxPlot.nExperiments == 1 && length(boxPlot.subject) < 10
    title(sprintf('%s: %d Tri, %s', boxPlot.subject, boxPlot.trialNumber, boxPlot.extraInf))
else, title(sprintf('%s: %d Tri, %d Sess', boxPlot.subject, boxPlot.trialNumber, boxPlot.nExperiments));
end
set(gca, 'xTick', 1:size(plotData,2), 'xTickLabel', boxPlot.xyValues{1}, 'fontsize', 14)
set(gca, 'yTick', 1:size(plotData,1), 'yTickLabel', boxPlot.xyValues{2}, 'fontsize', 14, 'TickLength', [0, 0])
box off;
end