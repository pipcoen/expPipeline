function inactivatedBrain(brainPlot, addText)
if ~exist('addText', 'var'); addText = 1; end
if ~isfield(brainPlot, 'brainImage'); brainPlot.brainImage=imread('BrainOutlineBW.png'); end

gridOpt.type = 'galvouni';
imagesc(linspace(-4.5,4.5,1000),linspace(3.75,-5.2,1000),brainPlot.brainImage); axis xy;
hold on;
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
[plotData, gridXY] = prc.makeGrid(filteredBlocks{2}, filteredBlocks{2}.responseMade==2, @mean, gridOpt);
numTrials = prc.makeGrid(filteredBlocks{2}, filteredBlocks{2}.responseMade==2, @length, gridOpt);

plotData(numTrials==0) = nan;
plotData = plotData - mean(filteredBlocks{1}.responseMade==2);
scatter(gridXY{1}(:), gridXY{2}(:), 150, plotData(:),'o','filled'); axis equal;  drawnow
grid on;
xlim([-4.5 4.5])
box off; set(gca, 'ycolor', 'w', 'xcolor', 'w', 'xTick', -4:1:4, 'yTick', -4:3, 'gridAlpha', 0.75, 'gridlinestyle', ':', 'GridColor', 'k', 'LineWidth', 1);
plot(0,0, 'pg', 'markersize', 10, 'markerfacecolor', 'g');

gridXY{1} = gridXY{1}(numTrials~=0);
gridXY{2} = gridXY{2}(numTrials~=0);
numTrials = numTrials(numTrials~=0);
if addText
    arrayfun(@(x,y,z) text(x,y, num2str(round(z*100)/100), 'horizontalalignment', 'center', 'VerticalAlignment', 'middle'), gridXY{1}(:), gridXY{2}(:), numTrials(:))
end
colormap(plt.redblue(64));
caxis([-1 1]);
end