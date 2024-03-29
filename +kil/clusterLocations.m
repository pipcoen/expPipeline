function clusterLocations = getClusterLoactionsInAllenSpace(clusterDepths, clusterPenetrationIdx)
probLength = 3840;
ephysRecord = load(prc.pathFinder('ephysrecord'), 'ephysRecord'); ephysRecord = ephysRecord.ephysRecord;
penetrations = unique(clusterPenetrationIdx);
penetrationDetails = ephysRecord(cell2mat(arrayfun(@(x) find([ephysRecord.penetrationIdx] == x),penetrations, 'uni', 0)));
tipLocations = cell2mat(arrayfun(@(x) penetrationDetails(x==penetrations).calcTip,clusterPenetrationIdx, 'uni', 0));
scaleFactors = cell2mat(arrayfun(@(x) penetrationDetails(x==penetrations).scalingFactor,clusterPenetrationIdx, 'uni', 0));
probeVectors = cell2mat(arrayfun(@(x) penetrationDetails(x==penetrations).calcLine(:)',clusterPenetrationIdx, 'uni', 0));
clusterOffsets = repmat(scaleFactors.*(probLength-clusterDepths),1,3);
clusterLocations.position = round(bsxfun(@minus, tipLocations,(clusterOffsets.*probeVectors)/10));
kil.loadAtlas;
clusterLocations.areaIdx = arrayfun(@(x,y,z) av(x,y,z), clusterLocations.position(:,3),clusterLocations.position(:,2),clusterLocations.position(:,1));
clusterLocations.areaName = st.acronym(clusterLocations.areaIdx);
parentIdx = st.parent_structure_id(clusterLocations.areaIdx);
parentIdx(parentIdx==0) = 997;
[~,parentIdx] = ismember(parentIdx,st.id);
clusterLocations.parent = st.acronym(parentIdx);
end