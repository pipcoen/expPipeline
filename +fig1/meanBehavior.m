function meanBehavior(behBlks)
%% This function plots the data panels for figure one of the ms
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behaviour', 0, 1); end
mTri = 1;
behBlks = behBlks.blks;
nMice = length(behBlks);
[rTurns, nTrials, reacT, reacTSE, rTurnsSE] = deal(nan*ones(3, 9, nMice));
visRef = [-1*[0.8, 0.4, 0.2, 0.1] 0 0.1, 0.2, 0.4, 0.8];
eIdx = strcmp(arrayfun(@(x) x.exp.subject{1}, behBlks, 'uni', 0), 'PC022');
%%
for i = 1:nMice
    nBlk = spatialAnalysis.getBlockType(behBlks(i), 'norm');
    nBlk = prc.filtBlock(nBlk, nBlk.exp.numOfTrials > mTri);
    nBlk = prc.filtBlock(nBlk, ismember(nBlk.tri.stim.visDiff, visRef));
    
    keepIdx = ismember(nBlk.exp.conditionParametersAV{1}(:,2), visRef);
    nBlk.exp.conditionParametersAV = cellfun(@(x) x(keepIdx,:), nBlk.exp.conditionParametersAV, 'uni', 0);
    nBlk.exp.conditionLabels = cellfun(@(x) x(keepIdx,:), nBlk.exp.conditionLabels, 'uni', 0);
    grds = prc.getGridsFromBlock(nBlk);
    
    rTurns(:,:,i) = grds.fracRightTurns;
    reacT(:,:,i) = grds.reactionTime;
    nTrials(:,:,i) = grds.numTrials;
    reacTSE(:,:,i) = grds.reactionTimeSE;
    rTurnsSE(:,:,i) = grds.fracRightTurnsSE;
end
fprintf('Example mouse trials: %d\n', (sum(nansum(nTrials(:,:,eIdx)))));
fprintf('Total trials: %d\n', (nansum(nTrials(:))));



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

%%
axesHandle = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
plotData = cat(3, rTurns(:,:,eIdx), rTurns(:,:,eIdx)-rTurnsSE(:,:,eIdx), rTurns(:,:,eIdx)+rTurnsSE(:,:,eIdx));
plt.rowsOfGrid(grds.visValues(1,:)*100, plotData, plt.selectRedBlueColors(grds.audValues(:,1)));
axis square

axesHandle = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg);cla;
ylim([150 300])
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
plotData = cat(3, reacT(:,:,eIdx), reacT(:,:,eIdx)-reacTSE(:,:,eIdx), reacT(:,:,eIdx)+reacTSE(:,:,eIdx));
plt.rowsOfGrid(grds.visValues(1,:)*100, plotData*1000, plt.selectRedBlueColors(grds.audValues(:,1)));
axis square
%%
axesHandle = plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); cla;
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
meanData = mean(rTurns,3);
seData = std(rTurns,[],3)./sqrt(nMice);
plotData = cat(3, meanData, meanData-seData, meanData+seData);
plt.rowsOfGrid(visRef(1,:)*100, plotData, plt.selectRedBlueColors([-60 0 60]));
axis square
%%
axesHandle = plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); cla
ylim([180 260])
xlim([-80 80])
set(gca, 'XTick', [-80 0 80]);
meanData = mean(reacT,3);
seData = std(reacT,[],3)./sqrt(nMice);
plotData = cat(3, meanData, meanData-seData, meanData+seData);
plt.rowsOfGrid(visRef(1,:)*100, plotData*1000, plt.selectRedBlueColors([-60 0 60]));
axis square

%% %%%%%%
% export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\1_meanBehavior', '-pdf', '-painters');
end