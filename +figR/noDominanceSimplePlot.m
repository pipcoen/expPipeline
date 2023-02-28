function noDominanceSimplePlot(behBlks)
%% This function plots the data panels for figure one of the ms
if ~exist('behBlks', 'var'); behBlks = spatialAnalysis('all', 'behaviour', 0, 1, ''); end
behBlks = behBlks.blks;
nMice = length(behBlks);
%%
for i = 1:nMice
    nBlk = spatialAnalysis.getBlockType(behBlks(i), 'norm');
    nBlk = prc.filtBlock(nBlk,~isinf(nBlk.tri.stim.audInitialAzimuth));
    
    grds = prc.getGridsFromBlock(nBlk, 1);
    grds.fracRightTurns = grds.fracRightTurnsComb;
    
    fracAudL = grds.fracRightTurns(grds.audValues < 0 & grds.visValues==0);
    fracAudR = grds.fracRightTurns(grds.audValues > 0 & grds.visValues==0);
    
    fracVisL = grds.fracRightTurns.*(grds.visValues < 0 & grds.audValues==0);
    fracVisL(fracVisL==0) = nan;
    fracVisR = grds.fracRightTurns.*(grds.visValues > 0 & grds.audValues==0);
    fracVisR(fracVisR==0) = nan;

    closestVL = nanmin(abs((1-fracVisL(:)) - fracAudR(:))) == abs((1-fracVisL) - fracAudR);
    closestVR = nanmin(abs((1-fracVisR(:)) - fracAudL(:))) == abs((1-fracVisR) - fracAudL);
    
    res.audR(i,1) = fracAudR;
    res.audL(i,1) = fracAudL;
    res.visR(i,1) = grds.fracRightTurns(closestVR);
    res.visL(i,1) = grds.fracRightTurns(closestVL);
    res.mulAudR(i,1) = grds.fracRightTurns(grds.audValues > 0 &  grds.visValues == grds.visValues(closestVL));
    res.mulAudL(i,1) = grds.fracRightTurns(grds.audValues < 0 &  grds.visValues == grds.visValues(closestVR));
    res.distAudR(i,1) = (1-res.visL(i,1)) - res.audR(i,1);
    res.distAudL(i,1) = (1-res.visR(i,1)) - res.audL(i,1);  
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


axH = plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;
cla

xDatN = num2cell([res.audR*0+1 res.mulAudR*0+2 res.visL*0+3],2);
yDatN = num2cell([res.audR res.mulAudR res.visL],2);

hold on;
cellfun(@(x,y) plot(x(1:2),y(1:2), 'r','HandleVisibility','off'), xDatN, yDatN);
cellfun(@(x,y) plot(x(2:3),y(2:3), 'r','HandleVisibility','off'), xDatN, yDatN);
yDatN = cell2mat(yDatN);
plot(axH, xDatN{1}(1:2), mean(yDatN(:,1:2)),'k', 'lineWidth', 3);
plot(axH, xDatN{1}(2:3), mean(yDatN(:,2:3)),'k', 'lineWidth', 3);

xDatN = num2cell([res.audR*0+5 res.mulAudR*0+6 res.visL*0+7],2);
yDatN = num2cell([res.audL res.mulAudL res.visR],2);
cellfun(@(x,y) plot(x(1:2),y(1:2), 'b','HandleVisibility','off'), xDatN, yDatN);
cellfun(@(x,y) plot(x(2:3),y(2:3), 'b','HandleVisibility','off'), xDatN, yDatN);
yDatN = cell2mat(yDatN);
plot(axH, xDatN{1}(1:2), mean(yDatN(:,1:2)),'k', 'lineWidth', 3);
plot(axH, xDatN{1}(2:3), mean(yDatN(:,2:3)),'k', 'lineWidth', 3);


box off;
ylim([0 1]);
yline(0.5,  '--k');
%% %%%%%%
export_fig('D:\OneDrive\Papers\Coen_2021\CellRevision\FigureParts\EasyDominancePlot', '-pdf', '-painters');
end
