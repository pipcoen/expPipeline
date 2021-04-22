function probesAndCellsInBrain
%% This function plots the data panels for figure one of the ms
s = spatialAnalysis('all', 'm2ephysgood', 1, 0, 'eph', 'multiSpaceWorld');
s.blks = prc.filtBlock(s.blks,~(s.blks.pen.numOfClusters==44));
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
kil.loadAtlas;
atlas.tv = tv;
atlas.st = st;
atlas.av = av;
s.blks  = kil.getClusterLoactionsInAllenSpace(s.blks, [], atlas);

probeLength = 3840;
[probeLines, pIdx, probeRef] = unique(s.blks.pen.calcLine, 'rows');
probeCount = groupcounts(probeRef);
probeTips = s.blks.pen.calcTip(pIdx,:);
scalingFactor = s.blks.pen.scalingFactor(pIdx);
colC = jet(6);
%%
plt.tightSubplot(nRows,nCols,1,axesGap,botTopMarg,lftRgtMarg); cla;
axH = gca;
hold on
kil.plotAllenOutlines('sag', {'MOs'}, axH)
for i = 1:length(pIdx)
    startPoint = probeTips(i,:)-(probeLength/10).*probeLines(i,:)*scalingFactor(i);
    probeStEn = [startPoint' probeTips(i,:)'];
    plot(axH, probeStEn(3,:), probeStEn(2,:), 'color', colC(probeCount(i),:), 'linewidth', 1.5)
end
colormap(jet(6))
colorbar
%%
plt.tightSubplot(nRows,nCols,2,axesGap,botTopMarg,lftRgtMarg); cla;
axH = gca;
hold on
kil.plotAllenOutlines('cor', {'MOs'}, axH)
for i = 1:length(pIdx)
    startPoint = probeTips(i,:)-(probeLength/10).*probeLines(i,:)*scalingFactor(i);
    probeStEn = [startPoint' probeTips(i,:)'];
    plot(axH, probeStEn(1,:), probeStEn(2,:), 'color', colC(probeCount(i),:), 'linewidth', 1.5)
end

plt.tightSubplot(nRows,nCols,3,axesGap,botTopMarg,lftRgtMarg); cla;
axH = gca;
hold on
kil.plotAllenOutlines('top', {'MOs'}, axH)
for i = 1:length(pIdx)
    startPoint = probeTips(i,:)-(probeLength/10).*probeLines(i,:)*scalingFactor(i);
    probeStEn = [startPoint' probeTips(i,:)'];
    plot(axH, probeStEn(1,:), probeStEn(3,:), 'color', colC(probeCount(i),:), 'linewidth', 1.5)
end

plt.tightSubplot(nRows,nCols,4,axesGap,botTopMarg,lftRgtMarg); cla;
colorLabels = repmat([0.5 0.5 0.5],s.blks.tot.clusters,1);
areaIdx = contains(s.blks.clu.parent, {'MOs'; 'FRP'});
colorLabels(areaIdx,:) = repmat([0 0 0], sum(areaIdx), 1);
plt.clustersInBrain(s.blks.clu, colorLabels,1);
%%
lfpEx = 49;
plt.tightSubplot(nRows,nCols,5,axesGap,botTopMarg,lftRgtMarg); cla;
powerSpectra = zscore(flipud(log10(s.blks.pen.lfpPowerSpectra{lfpEx}.powerSpectra(1:100,:)')));
freqPoints = s.blks.pen.lfpPowerSpectra{lfpEx}.freqPoints(1:100);

refChannels = [5; 44; 81; 120; 157; 196; 233; 272; 309; 348];
chan2Plot = setdiff(1:384, refChannels);

imagesc(freqPoints,1:length(chan2Plot), powerSpectra(chan2Plot,:));
colormap default
colorbar;
%%
export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\SupX_ProbesAndCellsInBrain', '-pdf', '-painters');
end