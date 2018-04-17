function scanningBrainEffects(scanPlot)
if ~isfield(scanPlot, 'addTrialNumber'); scanPlot.addTrialNumber = 0; end
plt.allenOutline;
hold on;
plotData = scanPlot.data(:);
MLCoord = scanPlot.gridXY{1}(:);
APCoord = scanPlot.gridXY{2}(:);
nTrials = scanPlot.nTrials(:);
if ~isfield(scanPlot, 'pVals')
    spotSize = APCoord*0+200;
else
    spotSize = scanPlot.pVals(:);
    spotSize = ((APCoord*0+1) + double(spotSize<0.01)*1.5 + double(spotSize<0.001)*1.5 + double(spotSize<0.0001)*1.5)*30;
end
scatter(MLCoord, APCoord, spotSize, plotData,'o', 'filled'); axis equal; drawnow
grid on;
xlim([-5.5 6])
ylim([-5.5 6])
box off; set(gca, 'ycolor', 'w', 'xcolor', 'w', 'xTick', -5:1:5, 'yTick', -5:4, 'gridAlpha', 0.75, 'gridlinestyle', ':', 'GridColor', 'k', 'LineWidth', 1);
plot(0,0, 'pg', 'markersize', 10, 'markerfacecolor', 'g');
colormap(plt.redblue(64));
caxis([-0.7 0.7]);
title(scanPlot.title);
if scanPlot.addTrialNumber
    vIdx = nTrials(:)>0;
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center', 'VerticalAlignment', 'middle'), MLCoord(vIdx), APCoord(vIdx), nTrials(vIdx))
end
end