function inactivatedBrain(brainPlot, addText)
if ~exist('addText', 'var'); addText = 1; end
if ~isfield(brainPlot, 'brainImage'); brainPlot.brainImage=imread('BrainOutlineBW.png'); end
bregma = [540,0,570];
load allenCorticalBoundaries.mat corticalAreadBoundaries
hold on;
for i =1:length(corticalAreadBoundaries)
    cellfun(@(x) plot((x(:,2)-bregma(3))/100, (bregma(1)-x(:,1))/100,'k'),corticalAreadBoundaries{i});
end
blockDat = {brainPlot.normBlock; brainPlot.laserBlock};
switch brainPlot.condition
    case 'VL'
        filterOpt = cellfun(@(x) x.trialType==2 & x.correctResponse==1, blockDat, 'uni', 0);
        title('Visual Left')
    case 'VR'
        filterOpt = cellfun(@(x) x.trialType==2 & x.correctResponse==2, blockDat, 'uni', 0);
        title('Visual Right')
    case 'AL'
        filterOpt = cellfun(@(x) x.trialType==1 & x.correctResponse==1, blockDat, 'uni', 0);
        title('Auditory Left')
    case 'AR'
        filterOpt = cellfun(@(x) x.trialType==1 & x.correctResponse==2, blockDat, 'uni', 0);
        title('Auditory Right')
    case 'CL'
        filterOpt = cellfun(@(x) x.trialType==4 & x.visInitialAzimuth>0, blockDat, 'uni', 0);
        title('Conflict: Visual Right')
    case 'CR'
        filterOpt = cellfun(@(x) x.trialType==4 & x.visInitialAzimuth<0, blockDat, 'uni', 0);
        title('Conflict: Visual Left')
end
filteredBlocks = cellfun(@(x,y) prc.combineBlocks(x, y), blockDat, filterOpt, 'uni', 0);
[plotData, gridXY] = prc.makeGrid(filteredBlocks{2}, filteredBlocks{2}.responseMade==2, @mean, 'galvouni');
numTrials = prc.makeGrid(filteredBlocks{2}, filteredBlocks{2}.responseMade==2, @length, 'galvouni');

plotData(numTrials==0) = nan;
plotData = plotData - mean(filteredBlocks{1}.responseMade==2);
scatter(gridXY{1}(:), gridXY{2}(:), 150, plotData(:),'o','filled'); axis equal;  drawnow
grid on;
xlim([-5.5 5.5])
ylim([-5.5 4])
box off; set(gca, 'ycolor', 'w', 'xcolor', 'w', 'xTick', -5:1:5, 'yTick', -5:4, 'gridAlpha', 0.75, 'gridlinestyle', ':', 'GridColor', 'k', 'LineWidth', 1);
plot(0,0, 'pg', 'markersize', 10, 'markerfacecolor', 'g');

gridXY{1} = gridXY{1}(numTrials~=0);
gridXY{2} = gridXY{2}(numTrials~=0);
numTrials = numTrials(numTrials~=0);
if addText
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center', 'VerticalAlignment', 'middle'), gridXY{1}(:), gridXY{2}(:), numTrials(:))
end
colormap(plt.redblue(64));
caxis([-0.7 0.7]);
end