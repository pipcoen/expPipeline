function inactivatedBrainPerturbations(plotData, gridXY, condition, numTrials, addText)
if ~exist('addText', 'var'); addText = 1; end
bregma = [540,0,570];
load allenCorticalBoundaries.mat corticalAreadBoundaries
hold on;
for i =1:length(corticalAreadBoundaries)
    cellfun(@(x) plot((x(:,2)-bregma(3))/100, (bregma(1)-x(:,1))/100,'k'),corticalAreadBoundaries{i});
end

switch condition
    case 'Bias'
        plotData = plotData(:,:,1);
    case 'VR-VL'
        plotData = plotData(:,:,3) - plotData(:,:,2);
    case 'AR-AL'
        plotData =plotData(:,:,4) - plotData(:,:,6);
end
title(condition);
scatter(gridXY{1}(:), gridXY{2}(:), 150, plotData(:),'o','filled'); axis equal;  drawnow
grid on;
xlim([-5.5 5.5])
ylim([-5.5 5])
box off; set(gca, 'ycolor', 'w', 'xcolor', 'w', 'xTick', -5:1:5, 'yTick', -5:4, 'gridAlpha', 0.75, 'gridlinestyle', ':', 'GridColor', 'k', 'LineWidth', 1);
plot(0,0, 'pg', 'markersize', 10, 'markerfacecolor', 'g');

gridXY{1} = gridXY{1}(numTrials~=0);
gridXY{2} = gridXY{2}(numTrials~=0);
numTrials = numTrials(numTrials~=0);
if addText
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center', 'VerticalAlignment', 'middle'), gridXY{1}(:), gridXY{2}(:), numTrials(:))
end
colormap(plt.redblue(64));
caxis(abs(max(plotData(:)))*[-1 1]);
end