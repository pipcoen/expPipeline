function scatterAltPlotsVisComparison(behBlksOrig, vis2Use, offset)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('behBlksOrig', 'var') || isempty(behBlksOrig); behBlksOrig = spatialAnalysis('all', 'm2ephysgood',0,1,''); end
if ~exist('vis2Use', 'var'); vis2Use = 0.8; end %We use 0.4/40% contrast for comparisons
if ~exist('offset', 'var'); offset = 0; end %We use 0.4/40% contrast for comparisons

blks2Use = behBlksOrig.blks;
if length(blks2Use) == 6; blks2Use(1) = []; end

%pre-assign performance and reaction structures with nans
[perf.vis, perf.visOnly] = deal(nan*ones(length(blks2Use), 1));
[reac.vis, reac.visOnly] = deal(perf.vis);

nMice = length(blks2Use);
for i = 1:nMice
    %"normalize" block--removes timeouts, laser, and nan-response trials
    nBlk = spatialAnalysis.getBlockType(blks2Use(i), 'norm');
    grds = prc.getGridsFromBlock(nBlk);
     
    perf.vis(i) = grds.performance(grds.visValues == vis2Use & grds.audValues == 0);
    perf.visOnly(i) = grds.performance(grds.visValues == vis2Use & isinf(grds.audValues));
    
    reac.vis(i) = mean(grds.reactionTime(abs(grds.visValues) == vis2Use & grds.audValues == 0));
    reac.visOnly(i) = grds.performance(grds.visValues == vis2Use & isinf(grds.audValues));

    tOut.vis(i) = nanmean(grds.fracTimeOutComb(abs(grds.visValues) == vis2Use & grds.audValues == 0));
    tOut.visOnly(i) = nanmean(grds.fracTimeOutComb(abs(grds.visValues) == vis2Use & grds.audValues == 0));
end
if any(sum(isnan([perf.vis, perf.visOnly, reac.vis, reac.visOnly]))) && vis2Use~=0.8
    error('Why are there nans???'); 
end
reac.offset = mean([reac.visOnly reac.vis],2)*offset;
%%
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
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);

axH = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg);
plotAltScatter([perf.vis perf.visOnly], perf.vis*0, axH)
ylim([0.5 1]);
box off;
title(['n = ' num2str(nMice)]);

[~, pVal] = ttest(perf.vis, perf.visOnly);
pVal = round(pVal, 2, 'significant');
text(0, 0.4, ['V vs VOnly: P<' num2str(pVal)]);


axH = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
plotAltScatter([reac.vis reac.visOnly], reac.offset, axH)
ylim([0.1 1.1]);
box off;
set(axH, 'xTickLabel', {'Vis'; 'VisOnly'}, 'XTick', [0.5 1.5])
title(['n = ' num2str(nMice)]);
ylabel('Relative reaction time (s)')
%%
[~, pVal] = ttest(reac.vis, reac.visOnly);
pVal = round(pVal, 2, 'significant');
text(0, -0.6, ['V vs VOnly: P<' num2str(pVal)]);

export_fig(['D:\OneDrive - University College London\Papers\Coen_2021\NeuronRevision\NewFigParts\' ...
    'ComparisonVisOnly'], '-pdf', '-painters');
end


function plotAltScatter(inDat, offset, axH)
nXPnts = size(inDat,2);
inDat = inDat - repmat(offset, 1, nXPnts);
yDat = cell2mat(arrayfun(@(x) [inDat(:,x); mean(inDat(:,x))], 1:nXPnts, 'uni', 0));
xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));

set(axH, 'position', get(axH, 'position').*[1 1 (0.2*nXPnts) 1]);
hold on
for i = 1:nXPnts-1
    cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(:,i:i+1),2), num2cell(yDat(:,i:i+1),2));
end

xDatN = xDat(1:end-1,:);
yDatN = yDat(1:end-1,:);
plot(axH, xDatN, yDatN,'ok', 'MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',5);
plot(axH, xDat(end,:), yDat(end,:),'^c', 'MarkerEdgeColor', 'c','MarkerFaceColor', 'c', 'MarkerSize',6);

xlim([0 nXPnts]);
end