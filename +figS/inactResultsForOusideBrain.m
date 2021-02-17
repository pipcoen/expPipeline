function inactResultsForOusideBrain
load('figS3InactResultsForChoice_IndiMice', 'inactResultsForChoice')
%%
subsets = inactResultsForChoice.subsets;
nMice = size(inactResultsForChoice.laserOnData,2);
[ipsiOut, contraOut] = deal(zeros(nMice, length(subsets)));
for j = 1:nMice
    for i = 1:length(subsets)        
        contData = inactResultsForChoice.meanContEffects{i,j};
        contData = contData(inactResultsForChoice.gridXY{1}{2}==5.5);
        ipsiOut(j,i) = contData(2);
        contraOut(j,i) = contData(end-1);
    end
end
for i = 1:length(subsets)
    [~, p] = ttest(ipsiOut(:,i))
    [~, p] = ttest(contraOut(:,i))
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
plotAltScatter(ipsiOut, ipsiOut(:,1)*0, axH)
ylim([-0.3 0.3]);
box off;
title(['Ipsi-Out n = ' num2str(nMice)]);

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;
set(gcf, 'position', get(gcf, 'position').*[1 1 0 0] + [0 0 figWidth, figHeight]);

axH = plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla;
plotAltScatter(contraOut, ipsiOut(:,1)*0, axH)
ylim([-0.3 0.3]);
box off;
title(['Contra-Out n = ' num2str(nMice)]);
%%

%"contData" is the result from the "normEstRepeats" loops, and "shuffleData" is from the shuffled loops. We then sort these shuffled loops and
%see where the control data appears in the shuffled data. This goes into the "scanPlot" plotting structure, along with the results.
export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\SupX_inactResultsForChoice_OutsideBrain', '-pdf', '-painters');
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