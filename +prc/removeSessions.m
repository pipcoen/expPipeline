function reducedBlock = removeSessions(blk, sessions2remove)
%% Function to combine and filter blk files.
nTrials = length(blk.sessionIdx);
nSessions = length(unique(blk.sessionIdx));
reducedBlock = blk;
fieldNames = fields(blk);
for fieldName = fieldNames'
    currField = fieldName{1};
    tDat = blk.(currField);
    if size(tDat,1) == nSessions; reducedBlock.(currField) = tDat(~ismember(1:nSessions, sessions2remove),:); end
    if size(tDat,1) == nTrials; reducedBlock.(currField) = tDat(~ismember(blk.sessionIdx, sessions2remove),:); end
end
reducedBlock.nSessions = length(unique(reducedBlock.sessionIdx));
reducedBlock.sessionIdx = cumsum(diff([0;reducedBlock.sessionIdx])>0);

if isfield(blk, 'ephSpikeSite')
    reducedBlock = kil.filtClusters(reducedBlock, ~ismember(blk.ephClusterSession, sessions2remove));
    reducedBlock.ephClusterSession = cumsum(diff([0;reducedBlock.ephClusterSession])>0);
    reducedBlock.ephSpikeSession = cumsum(diff([0;reducedBlock.ephSpikeSession])>0);
end
end
