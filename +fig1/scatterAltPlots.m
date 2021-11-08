function scatterAltPlots(behBlksOrig, vis2Use, offset)
%Load the block if it doesn't exist. Remove mice that have different parameter values (4 mice of 21)
if ~exist('behBlksOrig', 'var') || isempty(behBlksOrig); behBlksOrig = spatialAnalysis('all', 'behavior', 0, 1, ''); end
if ~exist('vis2Use', 'var'); vis2Use = 0.4; end %We use 0.4/40% contrast for comparisons
if ~exist('offset', 'var'); offset = 1; end %We use 0.4/40% contrast for comparisons

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
reac.offset = mean([reac.aud reac.vis reac.coh reac.con],2)*offset;
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
plotAltScatter([perf.aud perf.vis perf.mul], perf.vis*0, eIdx, axH)
ylim([0.5 1]);
box off;
title(['n = ' num2str(nMice)]);
%%
[~, pVal] = ttest(perf.vis, perf.aud);
pVal = round(pVal, 2, 'significant');
text(0, 0.4, ['A vs V: P<' num2str(pVal)]);

[~, pVal] = ttest(perf.aud, perf.mul);
pVal = round(pVal, 2, 'significant');
text(0, 0.35, ['A vs M: P<' num2str(pVal)]);

[~, pVal] = ttest(perf.vis, perf.mul);
pVal = round(pVal, 2, 'significant');
text(0, 0.3, ['V vs M: P<' num2str(pVal)]);

%%
axH = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); hold on;
plotAltScatter([reac.aud, reac.vis reac.coh reac.con], reac.offset, eIdx, axH)
title(['n=' num2str(nMice)]);
if nargin<2
    ylim([-0.05 0.05]);
else
    ylim([-0.06 0.06]);
end
if offset == 0
    ylim([0.1 0.35]);
end
box off;
%%
[~, pVal] = ttest(reac.aud, reac.vis);
pVal = round(pVal, 2, 'significant');
text(0, min(ylim)-0.08, ['A vs V: P<' num2str(pVal)]);

[~, pVal] = ttest(reac.coh, reac.con);
pVal = round(pVal, 2, 'significant');
text(0, min(ylim)-0.09, ['Coh vs Con: P<' num2str(pVal)]);

[~, pVal] = ttest(reac.aud, reac.coh);
pVal = round(pVal, 2, 'significant');
text(0, min(ylim)-0.10, ['A vs Coh: P<' num2str(pVal)]);

[~, pVal] = ttest(reac.vis, reac.coh);
pVal = round(pVal, 2, 'significant');
text(0, min(ylim)-0.11, ['V vs Coh: P<' num2str(pVal)]);

[~, pVal] = ttest(reac.aud, reac.vis);
pVal = round(pVal, 2, 'significant');
text(0, min(ylim)-0.12, ['A vs Con: P<' num2str(pVal)]);

[~, pVal] = ttest(reac.coh, reac.con);
pVal = round(pVal, 2, 'significant');
text(0, -0.13, ['V vs Con: P<' num2str(pVal)]);

%%
if nargin<2
    export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\1_scatterAltPlots_AudVer', '-pdf', '-painters');
end
end


function plotAltScatter(inDat, offset, eIdx, axH)
nXPnts = size(inDat,2);
inDat = inDat - repmat(offset, 1, nXPnts);
yDat = cell2mat(arrayfun(@(x) [inDat(~eIdx,x); inDat(eIdx,x); mean(inDat(:,x))], 1:nXPnts, 'uni', 0));
xDat = cell2mat(arrayfun(@(x) yDat(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));

set(axH, 'position', get(axH, 'position').*[1 1 (0.2*nXPnts) 1]);
hold on
for i = 1:nXPnts-1
    cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(:,i:i+1),2), num2cell(yDat(:,i:i+1),2));
end

xDatN = xDat(1:end-2,:);
yDatN = yDat(1:end-2,:);
plot(axH, xDatN, yDatN,'ok', 'MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',5);
plot(axH, xDat(end-1,:), yDat(end-1,:),'sc', 'MarkerEdgeColor', 'c','MarkerFaceColor', 'c', 'MarkerSize',6);
plot(axH, xDat(end,:), yDat(end,:),'^c', 'MarkerEdgeColor', 'c','MarkerFaceColor', 'c', 'MarkerSize',6);

xlim([0 nXPnts]);
end