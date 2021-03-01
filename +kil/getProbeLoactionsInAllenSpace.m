function getProbeLoactionsInAllenSpace(blk)
kil.loadAtlas;
probeLength = 3840;
i = 1;

%%
endPoint = blk.pen.calcTip(i,:);
startPoint = endPoint-(probeLength/10).*blk.pen.calcLine(i,:)*blk.pen.scalingFactor(i);
pointsAlongProbe = round(cell2mat(arrayfun(@(x,y) linspace(x,y,round(probeLength/5))', currStartPoint, currEndPoint, 'uni', 0)));

cDat = num2cell(guiData.refColorMap(length(guiData.handles.recordingProbe)),2);
pIdx = round(linspace(1,length(pointsAlongProbe),length(guiData.handles.recordingProbe)+1));
pltPnts = num2cell([pointsAlongProbe(pIdx(1:end-1),:),pointsAlongProbe(pIdx(2:end),:)],2);
cellfun(@(x,y,z) set(x,'XData',y([3 6]),'YData',y([1 4]),'ZData',y([2 5]),'color',z), num2cell(guiData.handles.recordingProbe)',pltPnts,cDat)

pointsAlongProbe = round(cell2mat(arrayfun(@(x,y) linspace(x,y,round(guiData.probeLength/5))', currStartPoint, currEndPoint, 'uni', 0)));
pointsAlongProbe(any(bsxfun(@gt, pointsAlongProbe, size(guiData.av)),2) | any(pointsAlongProbe<=0,2),:) = 1;
probeAreas = arrayfun(@(x,y,z) guiData.av(x,y,z), pointsAlongProbe(:,3),pointsAlongProbe(:,2),pointsAlongProbe(:,1));
probeAreaBoundaries = find(diff([1e5; double(probeAreas); 1e5])~=0);
probeAreaCenters = (probeAreaBoundaries(2:end) + probeAreaBoundaries(1:end-1))/2;
probeAreaCenters(probeAreaCenters>length(probeAreas)) = length(probeAreas);
probeAreaLabels = guiData.st.safe_name(probeAreas(round(probeAreaCenters)));
[~, ~, uniLocations] = unique(probeAreaLabels, 'stable');
newAreaBoundaries = find(diff([-1; uniLocations;-1])~=0);
probeAreaLabels = probeAreaLabels(newAreaBoundaries(1:end-1));
for i = 1:length(newAreaBoundaries)-1
    idx2Avergage = newAreaBoundaries(i):newAreaBoundaries(i+1)-1;
    probeAreaCenters(idx2Avergage) = mean(probeAreaCenters(idx2Avergage));
end
probeAreaCenters = probeAreaCenters(newAreaBoundaries(1:end-1));

ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
ephysRecordForBrain = ephysRecord(contains({ephysRecord.subject}', ephysRecord(guiData.currPenetration).subject));
averageScalingForBrain = round(nanmean([ephysRecordForBrain.scalingFactor])*100)/100;

% Update the probe areas
yyaxis(guiData.handles.regionsAxes,'right');
set(guiData.handles.probeRegionsPlot,'YData',(1:length(probeAreas))*5,'CData',probeAreas);
set(guiData.handles.regionsAxes,'YTick',probeAreaCenters*5,'YTickLabels',probeAreaLabels);
set(guiData.probeDetails, 'String', sprintf('SF: %0.2f (%s)', guiData.currScalingFactor, num2str(averageScalingForBrain)));
% Upload guiData
guidata(probeDepthGui, guiData);
end





