function scanningBrainEffects(scanPlot)
if ~isfield(scanPlot, 'addTrialNumber'); scanPlot.addTrialNumber = 0; end
plt.allenOutline;
hold on;
plotData = scanPlot.data(:);
MLCoord = scanPlot.gridXY{1}(:);
APCoord = scanPlot.gridXY{2}(:);
if ~isfield(scanPlot, 'pVals')
    spotSize = APCoord*0+200;
else
    spotSize = scanPlot.pVals(:);
    spotSize = ((APCoord*0+1) + double(spotSize<0.001) + double(spotSize<0.0005)*2 + double(spotSize<0.0001)*4);
end
spotSize = spotSize*(200./max(spotSize(:)));
scatter(MLCoord, APCoord, spotSize, plotData,'o', 'filled'); axis equal; drawnow
grid on;
xlim([-5.5 6])
ylim([-5.5 6])
box off; set(gca, 'ycolor', 'w', 'xcolor', 'w', 'xTick', -5:1:5, 'yTick', -5:4, 'gridAlpha', 0.75, 'gridlinestyle', ':', 'GridColor', 'k', 'LineWidth', 1);
% plot(0,0, 'dk', 'markersize', 10, 'markerfacecolor', 'k');
colormap(plt.redblue(64));
if isfield(scanPlot, 'colorBarLimits'); caxis(scanPlot.colorBarLimits); else, caxis([-0.7 0.7]); end

if isfield(scanPlot, 'title'); title(scanPlot.title); end
if scanPlot.addTrialNumber
    nTrials = scanPlot.nTrials(:);
    vIdx = nTrials(:)>0;
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center', 'VerticalAlignment', 'middle'), MLCoord(vIdx), APCoord(vIdx), nTrials(vIdx))
end
end