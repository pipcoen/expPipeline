function scatterPlots(behBlksOrig)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('behBlks', 'var'); behBlksOrig = spatialAnalysis('all', 'behavior', 0, 1); end
blks2Use = behBlksOrig.blks;
if length(blks2Use) == 21; blks2Use(5:8) = []; end

%pre-assign performance and reaction structures with nans
[perf.aud, perf.vis, perf.mul] = deal(nan*ones(length(blks2Use), 1));
[reac.aud, reac.vis, reac.coh, reac.con] = deal(perf.aud);
allRTs = cell(length(blks2Use), 1);
mTri = 1; %Can be uses to set a minimum number of trials per experiment

vis2Use = 0.4; %We use 0.4/40% contrast for comparisons
eIdx = strcmp(arrayfun(@(x) x.exp.subject{1}, blks2Use, 'uni', 0), 'PC022');
nMice = length(blks2Use);
for i = 1:nMice
    %"normalize" block--removes timeouts, laser, and nan-response trials
    nBlk = spatialAnalysis.getBlockType(blks2Use(i), 'norm');
    nBlk = prc.filtBlock(nBlk, nBlk.exp.numOfTrials > mTri);    
    grds = prc.getGridsFromBlock(nBlk);
    
    perf.aud(i) = grds.performance(grds.visValues == 0 & grds.audValues > 0);
    perf.vis(i) = grds.performance(grds.visValues == vis2Use & grds.audValues == 0);
    perf.mul(i) = grds.performance(grds.visValues == vis2Use & grds.audValues > 0);
    
    reac.aud(i) = mean(grds.reactionTime(grds.visValues == 0 & grds.audValues ~= 0));
    reac.vis(i) = mean(grds.reactionTime(abs(grds.visValues) == vis2Use & grds.audValues == 0));
    allRTs{i,1} = nBlk.tri.outcome.reactionTime;
    
    cohIdx = grds.visValues.*grds.audValues > 0 &~isnan(grds.reactionTime) & abs(grds.visValues) == vis2Use;
    conIdx = grds.visValues.*grds.audValues < 0 &~isnan(grds.reactionTime) & flipud(cohIdx) & abs(grds.visValues) == vis2Use;
    reac.coh(i) = mean(grds.reactionTime(cohIdx));
    reac.con(i) = mean(grds.reactionTime(conIdx));
end
if any(sum(isnan([perf.aud, perf.vis, perf.mul, reac.aud, reac.vis, reac.coh, reac.con]))); error('Why are there nans???'); end

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
hold on
scatter(axH, perf.aud(~eIdx), perf.vis(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, perf.aud(eIdx), perf.vis(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(perf.aud), mean(perf.vis), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(perf.aud, perf.vis);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
xlim([0.5 1]);
ylim([0.5 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;

%%
axH = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); hold on;
scatter(axH, perf.mul(~eIdx), perf.vis(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, perf.mul(eIdx), perf.vis(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(perf.mul), mean(perf.vis), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(perf.mul, perf.vis);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
xlim([0.5 1]);
ylim([0.5 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;

%%
axH = plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); hold on;
histogram(axH, cell2mat(allRTs),100, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 1);
xlim([0 1.5]);
axis square; 
box off;
%%
axH = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); hold on;
scatter(axH, reac.aud(~eIdx,1), reac.vis(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.aud(eIdx, 1), reac.vis(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(reac.aud), mean(reac.vis), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(reac.aud, reac.vis);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
axis square; 
xlim([0.1 0.35])
ylim([0.1 0.35])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;

%%

axH = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); hold on;
cellfun(@(x,y) plot(x,[y, y], 'k'), num2cell([reac.coh reac.con],2), num2cell(reac.vis));
plot(mean([reac.coh reac.con]), mean(reac.vis).*[1 1], 'k')

scatter(axH, reac.coh(~eIdx,1), reac.vis(~eIdx),25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.con(~eIdx), reac.vis(~eIdx), 25, 'w', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.coh(eIdx), reac.vis(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.con(eIdx), reac.vis(eIdx), 50, 'w', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(reac.coh), mean(reac.vis), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(reac.con), mean(reac.vis), 25, 'w', '^', 'filled', 'MarkerEdgeColor', 'k');


[~, pVal] = ttest(reac.coh, reac.con);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
axis square; 
hold on
xlim([0.1 0.35])
ylim([0.1 0.35])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;

%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\1_scatterPlots', '-pdf', '-painters');
end