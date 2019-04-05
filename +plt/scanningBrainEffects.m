function scanningBrainEffects(scanPlot)
if ~isfield(scanPlot, 'addTrialNumber'); scanPlot.addTrialNumber = 0; end
if ~isfield(scanPlot, 'pVals'); scanPlot.pVals = 0; end
plt.allenOutline;
hold on;
plotData = scanPlot.data(:);
MLCoord = scanPlot.gridXY{1}(:);
APCoord = scanPlot.gridXY{2}(:);
pVals = scanPlot.pVals(:);
minSig = min(pVals);
sigLevels = [0.01; 0.001; 0.0001];
if ~scanPlot.pVals || minSig > 0.01
    warning('Pvals do not exist or are all > 0.01');
    spotSize = APCoord*0+200;
else
    spotSize = sum(cell2mat(arrayfun(@(x) pVals<x, sigLevels, 'uni', 0)'),2);
    legendRef = sigLevels(sigLevels>minSig);
    lengendSizes = (1:sum(sigLevels>minSig))'.*(200./max(spotSize(:)));
    spotSize(spotSize==0) = 0.2;
    spotSize = spotSize*(200./max(spotSize(:)));
    legendRef = [legendRef, legendRef*0+4, legendRef*0+6-(0:length(legendRef)-1)' lengendSizes];
end
sigIdx = spotSize~=min(spotSize);
h1 = scatter(MLCoord(sigIdx), APCoord(sigIdx), spotSize(sigIdx), plotData(sigIdx), 'o', 'filled'); axis equal; drawnow
set(h1, 'MArkerEdgeColor', 'k');
scatter(MLCoord(~sigIdx), APCoord(~sigIdx), spotSize(~sigIdx), plotData(~sigIdx), 'o', 'filled'); axis equal; drawnow
if exist('legendRef', 'var')
    scatter(legendRef(:,2), legendRef(:,3), legendRef(:,4), 'k', 'o', 'filled'); axis equal; drawnow
    arrayfun(@(x,y,z) text(x,y, num2str(z), 'VerticalAlignment', 'middle'), legendRef(:,2)+0.5, legendRef(:,3), legendRef(:,1))
end

grid on;
xlim([-5.5 6])
ylim([-5.5 6])
box off; set(gca, 'ycolor', 'w', 'xcolor', 'w', 'xTick', -5:1:5, 'yTick', -5:4, 'gridAlpha', 0.75, 'gridlinestyle', ':', 'GridColor', 'k', 'LineWidth', 1);
colormap(plt.redblue(64));
if isfield(scanPlot, 'colorBarLimits'); caxis(scanPlot.colorBarLimits); else, caxis([-0.7 0.7]); end

if isfield(scanPlot, 'title'); title(scanPlot.title); end
if scanPlot.addTrialNumber
    nTrials = scanPlot.nTrials(:);
    vIdx = nTrials(:)>0;
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center', 'VerticalAlignment', 'middle'), MLCoord(vIdx), APCoord(vIdx), nTrials(vIdx))
end
end