function exampleInactModelFitsSignificantTest
%%
uniMice = {'PC027'; 'PC029'; 'DJ008'; 'DJ006'; 'DJ007'};
logLik = nan*(ones(5,5,length(uniMice)));
for i = 1:length(uniMice)
    s = spatialAnalysis(uniMice{i}, 'uniscan',1,1);
    [logLik(:,:,i), fitRegion, testRegion] = s.compareLLRInactivationAcrossRegions;
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
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 200 figWidth, figHeight]);

%%
plt.tightSubplot(nRows,nCols,1:2,axesGap,botTopMarg,lftRgtMarg); cla;
logLikNR = logLik*0;
for i = 1:5
    for j = 1:5
        if i == j; continue; end
%         [nH(i,j),pVal(i,j)] = ttest(logLik(i,i,:), logLik(j,i,:));
        logLikR(:,:,i) = bsxfun(@minus, logLik(:,:,i), diag(logLik(:,:,i))')*-1;
    end
end

for i = 1:5
    for j = 1:5
        logLikNR(i,j,:) = (logLikR(j,i,:)+logLikR(i,j,:))./2;
        [nH(i,j),pVal(i,j)] = ttest(logLikR(i,i,:), (logLikR(j,i,:)+logLikR(i,j,:))./2);
    end
end

%%
plotDat = mean(logLikR,3);
imagesc(plotDat)
axis square
set(gca, 'XTick', [])
set(gca, 'YTick', 1:5, 'YTickLabels', fitRegion(:,1))
ylabel('Inactivation site 1');
xlabel('Inactivation site 2');
box off
colormap('gray'); colorbar; caxis([-0.4 0])

plt.tightSubplot(nRows,nCols,4:5,axesGap,botTopMarg,lftRgtMarg); cla;
set(gca, 'position', get(gca, 'position').*[1 1.5 1 1])
plotDat = mean(logLikR,3);
plotDat = mean(cat(3,plotDat, plotDat'), 3);
plotDat = triu(plotDat);

imagesc(plotDat)
axis square
set(gca, 'XTick', 1:5, 'XTickLabels', fitRegion(:,1))
set(gca, 'YTick', 1:5, 'YTickLabels', fitRegion(:,1))
ylabel('Inactivation site 1');
xlabel('Inactivation site 2');
box off
colormap('gray'); colorbar; caxis([-0.4 0])

clear offDiag
for i = 1:5
    tDat = logLikR(1:4,1:4,i);
    offDiag(i,1) =  mean(tDat(~eye(4)));
end
[~, pVal] = ttest(offDiag)
%%
plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); cla;
LLRRelNorm = squeeze(permute(logLikR(5,:,:), [3,2,1]))*-1;

yDat = LLRRelNorm(:,1:4);
nXPnts = size(yDat,2);
xDat = cell2mat(arrayfun(@(x) LLRRelNorm(:,1)*0+x-0.5, 1:nXPnts, 'uni', 0));

set(gca, 'position', get(gca, 'position').*[1 1 (0.2*nXPnts) 1]);
hold on
for i = 1:nXPnts-1
    cellfun(@(x,y) plot(x,y, 'k','HandleVisibility','off'), num2cell(xDat(:,i:i+1),2), num2cell(yDat(:,i:i+1),2));
end

plot(gca, xDat, yDat,'ok', 'MarkerEdgeColor', 'none','MarkerFaceColor', 'k', 'MarkerSize',5);
xlim([0 nXPnts]);
set(gca, 'XTick', 0.5:1:3.5, 'XTickLabels', fitRegion(1:4,1))
ylabel('LLR vs no inactivation')
[~,~,stats] = anova2(yDat);
%%
multcompare(stats, 'CType', 'lsd');
%%
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\sigTest4InactivationSites', '-pdf', '-painters');
