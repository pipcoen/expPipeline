function reducedBlock = filtClusters(blk, clusterCriterion)
%% Function to combine and filter blk files.
nSpikes = length(blk.ephSpikeAmps);
nClusters = length(blk.ephClusterAmps);
fieldNames = fields(blk);
reducedBlock = blk;

spikeCriterion = ismember(blk.ephSpikeCluster, find(clusterCriterion));
for fieldName = fieldNames'
    currField = fieldName{1};
    tDat = blk.(currField);
    if size(tDat,1) == nSpikes; reducedBlock.(currField) = tDat(spikeCriterion,:); end
    if size(tDat,1) == nClusters; reducedBlock.(currField) = tDat(clusterCriterion,:); end
end
[~,reducedBlock.ephSpikeCluster] = ismember(reducedBlock.ephSpikeCluster, find(clusterCriterion));
end
