function behaviour4Model(ephBehBlks)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('ephBehBlks', 'var'); ephBehBlks = spatialAnalysis('all', 'm2ephysmod', 0, 1, 'raweph'); end

%pre-assign performance and reaction structures with nans
[perf.aud, perf.vis, perf.visO, perf.mul] = deal(nan*ones(length(ephBehBlks), 1));
[reac.aud, reac.vis, reac.visO, reac.coh, reac.con] = deal(perf.aud);
allRTs = cell(length(blks2Use), 1);
mTri = 1; %Can be uses to set a minimum number of trials per experiment

eIdx = strcmp(arrayfun(@(x) x.exp.subject{1}, blks2Use, 'uni', 0), 'PC022');
nMice = length(ephBehBlks.blks);
%%
for i = 1:nMice
    %"normalize" block--removes timeouts, laser, and nan-response trials
    nBlk = spatialAnalysis.getBlockType(ephBehBlks.blks(i), 'norm');
%     nBlk = prc.filtBlock(nBlk, ~isinf(nBlk.tri.stim.audInitialAzimuth));
    grds = prc.getGridsFromBlock(nBlk);
    nVisOnlyTrials(i) = sum(grds.numTrials(grds.visValues == 0.8 & isinf(grds.audValues)));
    
    perf.aud(i) = nanmean(grds.performance(grds.visValues == 0 & grds.audValues > 0));
    perf.vis(i) = grds.performance(grds.visValues == 0.8 & grds.audValues == 0);
    perf.visO(i) = grds.performance(grds.visValues == 0.8 & isinf(grds.audValues));
    
    reac.aud(i) = nanmean(grds.reactionTimeComb(grds.visValues == 0 & grds.audValues == 0))*1000;
    reac.vis(i) = mean(grds.reactionTimeComb(abs(grds.visValues) == 0.8 & grds.audValues==0))*1000;
    reac.visO(i) = mean(grds.reactionTimeComb(abs(grds.visValues) == 0.8 & isinf(grds.audValues)))*1000;
    
%     cohIdx = grds.visValues.*grds.audValues > 0 &~isnan(grds.reactionTime) & abs(grds.visValues) == vis2Use;
%     conIdx = grds.visValues.*grds.audValues < 0 &~isnan(grds.reactionTime) & flipud(cohIdx) & abs(grds.visValues) == vis2Use;
%     reac.coh(i) = mean(grds.reactionTime(cohIdx));
%     reac.con(i) = mean(grds.reactionTime(conIdx));
end
% if any(sum(isnan([perf.aud, perf.vis, perf.mul, reac.aud, reac.vis, reac.coh, reac.con]))); error('Why are there nans???'); end


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
scatter(axH, perf.aud, perf.vis, 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
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
axH = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);
hold on
scatter(axH, reac.aud, reac.vis, 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(reac.aud, reac.vis);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
xlim([150 400]);
ylim([150 400]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;
%%
axH = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg);
hold on
scatter(axH, perf.aud, perf.visO, 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(perf.aud, perf.visO);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
xlim([0.5 1]);
ylim([0.5 1]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;
%%
axH = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg);
hold on
scatter(axH, reac.aud, reac.visO, 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(reac.aud, reac.visO);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
xlim([180 800]);
ylim([180 800]);
axis square; 
hold on
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;
%%
axH = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); hold on;
scatter(axH, reac.aud(~eIdx,1), reac.vis, 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.aud(eIdx, 1), reac.vis, 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
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

scatter(axH, reac.coh(~eIdx,1), reac.vis,25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.con, reac.vis, 25, 'w', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.coh, reac.vis, 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.con, reac.vis, 50, 'w', 's', 'filled', 'MarkerEdgeColor', 'k');
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