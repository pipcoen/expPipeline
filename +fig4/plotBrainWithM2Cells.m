function plotBrainWithM2Cells
%% This function plots the data panels for figure one of the ms
s = spatialAnalysis('all', 'm2ephysgood', 1, 1);
%%
kil.loadAtlas;
atlas.tv = tv;
atlas.st = st;
atlas.av = av;
s.blks  = kil.getClusterLoactionsInAllenSpace(s.blks, [], atlas);
colorLabels = repmat([0.5 0.5 0.5],s.blks.tot.clusters,1);
areaIdx = contains(s.blks.clu.parent, {'MOs'; 'FRP'});
colorLabels(areaIdx,:) = repmat([0 0 0], sum(areaIdx), 1);
plt.clustersInBrain(s.blks.clu, colorLabels);

export_fig('D:\OneDrive\Papers\Coen_2020\FigureParts\4_BrainWithCells', '-pdf', '-painters');
end