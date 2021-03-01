function scatterPlots(behBlksOrig, vis2Use)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('behBlksOrig', 'var') || isempty(behBlksOrig); behBlksOrig = spatialAnalysis('all', 'behavior', 0, 1); end
if ~exist('vis2Use', 'var'); vis2Use = 0.4; end %We use 0.4/40% contrast for comparisons

blks2Use = behBlksOrig.blks;
if length(blks2Use) == 21; blks2Use(5:8) = []; end

%pre-assign performance and reaction structures with nans
[perf.vis, perf.aud, perf.mul] = deal(nan*ones(length(blks2Use), 1));
[reac.vis, reac.aud, reac.coh, reac.con] = deal(perf.vis);
allRTs = cell(length(blks2Use), 1);
mTri = 1; %Can be uses to set a minimum number of trials per experiment

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
    conIdx = grds.visValues.*grds.audValues < 0 &~isnan(grds.reactionTime) & abs(grds.visValues) == vis2Use;
    reac.coh(i) = mean(grds.reactionTime(cohIdx));
    reac.con(i) = mean(grds.reactionTime(conIdx));
end
if any(sum(isnan([perf.vis, perf.aud, perf.mul, reac.vis, reac.aud, reac.coh, reac.con]))) && vis2Use~=0.8
    error('Why are there nans???'); 
end

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
scatter(axH, perf.vis(~eIdx), perf.aud(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, perf.vis(eIdx), perf.aud(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(perf.vis), mean(perf.aud), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(perf.vis, perf.aud);
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
scatter(axH, perf.mul(~eIdx), perf.aud(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, perf.mul(eIdx), perf.aud(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(perf.mul), mean(perf.aud), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(perf.mul, perf.aud);
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
title([num2str(round(mean(cell2mat(allRTs)>0.5)*1000)/10) '% of tri > 0.5s']);
axis square; 
box off;
%%
axH = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); hold on;
scatter(axH, reac.vis(~eIdx,1), reac.aud(~eIdx), 25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.vis(eIdx, 1), reac.aud(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(reac.vis), mean(reac.aud), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
[~, pVal] = ttest(reac.vis, reac.aud);
pVal = round(pVal, 4, 'significant');
title(['n=' num2str(nMice) '   P<' num2str(pVal)]);
axis square; 
xlim([0.1 0.35])
ylim([0.1 0.35])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2);
box off;

%%

axH = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); cla; hold on;
cellfun(@(x,y) plot(x,[y, y], 'k','HandleVisibility','off'), num2cell([reac.coh reac.con],2), num2cell(reac.aud));
plot(mean([reac.coh reac.con]), mean(reac.aud).*[1 1], 'k')

scatter(axH, reac.coh(~eIdx,1), reac.aud(~eIdx),25, 'k', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.con(~eIdx), reac.aud(~eIdx), 25, 'w', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.coh(eIdx), reac.aud(eIdx), 50, 'k', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, reac.con(eIdx), reac.aud(eIdx), 50, 'w', 's', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(reac.coh), mean(reac.aud), 25, 'k', '^', 'filled', 'MarkerEdgeColor', 'k');
scatter(axH, mean(reac.con), mean(reac.aud), 25, 'w', '^', 'filled', 'MarkerEdgeColor', 'k');

[~, pVal] = ttest(reac.coh, reac.con);
pValAll{1,1} = round(pVal, 4, 'significant');
[~, pVal] = ttest(reac.aud, reac.coh);
pValAll{2,1} = round(pVal, 4, 'significant');
[~, pVal] = ttest(reac.aud, reac.con);
pValAll{3,1} = round(pVal, 4, 'significant');
lab2Use = {'Diff'; 'VvsCoh'; 'VvsCon'};

box off;
legendCell = arrayfun(@(x,y) [x{1} ': p <' num2str(round(y{1}, 2, 'significant'))], lab2Use, pValAll,'uni',0);
legend(legendCell, 'location', 'none', 'FontSize', 10, 'Position', get(gca, 'Position').*[1.5, 1.3, 1, 0.3])
legend('boxoff')

axis square; 
hold on
xlim([0.1 0.35])
ylim([0.1 0.35])
plot([min(xlim) max(xlim)], [min(ylim) max(ylim)], '--k', 'linewidth', 2,'HandleVisibility','off');
box off;

%%
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\1_scatterPlots_AudVer', '-pdf', '-painters');
end