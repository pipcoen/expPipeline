function exampleInactModelFitsSignificantTestShuffle
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
load inactCompResults215026;

logLikR = zeros(5,5);
for i = 1:5
    for j = 1:5
        if i == j; continue; end
        idx = ismember(grpOrd, [i,j], 'rows');
        sIdx = normEstRepeats+1:length(logLikSelfTest{idx});
        
        LLRSelf = mean(logLikSelfTest{idx}(1:normEstRepeats));
        LLRTest = mean(logLikTest{idx}(1:normEstRepeats));       
        shuffTest = sort([LLRSelf-LLRTest; logLikSelfTest{idx}(sIdx)-logLikTest{idx}(sIdx)]);
        
        pVal(i,j) = find(shuffTest == LLRSelf-LLRTest)/length(sIdx);        
        logLikR(i,j) = LLRTest - LLRSelf;
    end
end
imagesc(logLikR*-1)
axis square
set(gca, 'XTick', 1:5, 'XTickLabels', posNames)
set(gca, 'YTick', 1:5, 'YTickLabels', posNames)
ylabel('Trained site');
xlabel('Test site');
box off
colormap('gray'); colorbar; caxis([-0.1 0])


plt.tightSubplot(nRows,nCols,4:5,axesGap,botTopMarg,lftRgtMarg); cla;
set(gca, 'position', get(gca, 'position').*[1 1.5 1 1])
plotDat = mean(logLikR*-1,3);
plotDat = mean(cat(3,plotDat, plotDat'), 3);
plotDat = triu(plotDat);

imagesc(plotDat)
axis square
set(gca, 'XTick', 1:5, 'XTickLabels', posNames)
set(gca, 'YTick', 1:5, 'YTickLabels', posNames)
ylabel('Inactivation site 1');
xlabel('Inactivation site 2');
box off
colormap('gray'); colorbar; caxis([-0.1 0])

%%
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\sigTest4InactivationSitesShuffle', '-pdf', '-painters');
