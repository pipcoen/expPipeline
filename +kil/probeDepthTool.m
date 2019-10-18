function probeDepthTool(subjects)
%%
if ~iscell(subjects); subjects = {subjects}; end
allenPath = prc.pathFinder('allenAtlas');
cmap = load([fileparts(which('allenCCFbregma')) '\allen_ccf_colormap_2017.mat']);
expList = load(prc.pathFinder('expList')); expList = expList.expList;
%%
ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
selectedRecords = ephysRecord(contains({ephysRecord.subject}', subjects));
expList = expList(contains({expList.subject}', subjects) & [expList.excluded]' ~=1 & contains({expList.expType}', 'eph'));

guiData = struct;
guiData.expList = expList;
guiData.penetrations = find(contains({ephysRecord.subject}', subjects));
guiData.tv = readNPY([allenPath 'template_volume_10um.npy']); % grey-scale "background signal intensity"
guiData.av = readNPY([allenPath 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
guiData.st = loadStructureTree([allenPath 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
guiData.cmap = cmap.cmap; % Atlas colormap
guiData.bregma = [540,0,570];
guiData.probeLength = 3840; % Length of probe (um)
guiData.plottedStructuresIdx = []; % Plotted structures
guiData.probeAngle = [0;90]; % Probe angles in ML/DV

probeDepthGui = figure('color','w');
plotR = 3;
plotC = 6;
guiData.nPnts = 20;
guiData.refColorMap = @lines;

% Set up the atlas axes
brainSliceAxes = subplot(plotR,plotC,[1:plotR, (1:plotR)+plotC, (1:plotR)+plotC*2]); hold on;
slicePlot = surface('EdgeColor','none'); % Slice on 3D atlas
[~, guiData.handles.brainOutline] = plotBrainGrid([],brainSliceAxes);
hold(brainSliceAxes,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[apMax,dvMax,mlMax] = size(guiData.tv);
zoomLev = 50;
xlim([zoomLev,apMax-zoomLev]);
ylim([zoomLev,mlMax-zoomLev]);
zlim([zoomLev,dvMax-zoomLev])
probeGuideLine = plot3(rand(1), rand(1), rand(1), '--k', 'linewidth', 1.5);
histologyPoints = plot3(rand(1), rand(1), rand(1), '.b','MarkerSize',35);
recordingProbe = arrayfun(@(x) plot3(rand(1), rand(1), rand(1), '-m', 'linewidth', 8), 1:guiData.nPnts);

% Make 3D rotation the default state (toggle on/off with 'r')
h = rotate3d(brainSliceAxes);
h.Enable = 'on';
% Update the slice whenever a rotation is completed
h.ActionPostCallback = @updateSlice;

% Set functions for key presses
hManager = uigetmodemanager(probeDepthGui);
[hManager.WindowListenerHandles.Enabled] = deal(false);
set(probeDepthGui,'KeyPressFcn',@keyPress);
set(probeDepthGui,'KeyReleaseFcn',@keyRelease);

% Set up the probe area axes
regionsAxes = subplot(plotR,plotC,4:plotC:plotR*plotC);
regionsAxes.ActivePositionProperty = 'position';
set(regionsAxes,'FontSize',11);
yyaxis(regionsAxes,'left');
probeRegionsPlot = image(0);
set(regionsAxes,'XTick','','YLim',[0,3840],'YColor','k','YDir','reverse');
ylabel(regionsAxes,'Depth (\mum)');
colormap(regionsAxes,guiData.cmap);
caxis([1,size(guiData.cmap,1)])
yyaxis(regionsAxes,'right');
set(regionsAxes,'XTick','','YLim',[0,3840],'YColor','k','YDir','reverse');

% Set up units axis
unitAxes = subplot(plotR,plotC,5:plotC:plotR*plotC);
hold(unitAxes,'on');
set(unitAxes,'YDir','reverse', 'yAxisLocation', 'right', 'xcolor', 'w','XTick','','YLim',[0,3840],'YColor','w', 'color', 'none','YTick','');
guideCellScatter = plot(unitAxes,0,0,'.r','MarkerSize',10);
cellScatter = plot(unitAxes,0,0,'.k','MarkerSize',15);


cDat = num2cell(guiData.refColorMap(guiData.nPnts),2);
linIdx = round(linspace(0, 3840, guiData.nPnts+1));
arrayfun(@(x,y) plot([0 0],[linIdx(x) linIdx(x+1)],'color',y{1}, 'linewidth', 6), (1:guiData.nPnts)',cDat)
xlim([-0.1,1]);
ylabel('Depth (\mum)')


% CorrAxis
corrAxis = subplot(plotR,plotC,6:plotC:plotR*plotC);
hold(corrAxis,'on');
set(corrAxis,'XTick','','YLim',[0,3840],'xlim',[1,200], 'YTick', '','YDir','reverse');
lfpSpectrumPlot = imagesc(0);
guiData.currPlotIdx = 1;

linkaxes([unitAxes,regionsAxes],'y');

% Position the axes
set(brainSliceAxes,'Position',[0.1,0.1,0.5,0.9]);
set(regionsAxes,'Position',[0.7,0.1,0.03,0.8]);
set(unitAxes,'Position',[0.7,0.1,0.03,0.8]);
set(corrAxis,'Position',[0.8,0.1,0.12,0.8]);

% Probe information
[~, order2Process] = sortrows([isnan([selectedRecords.scalingFactor]') [selectedRecords.probePainted]'], 'descend');
guiData.penetrations = guiData.penetrations(order2Process);
guiData.currPenetration = guiData.penetrations(1);
guiData.currProbeVector = [];
guiData.currEntryPoint = [];
guiData.probeStepSize = 20;
guiData.currScalingFactor = 1;

%set up all handles
guiData.handles.regionsAxes = regionsAxes;
guiData.handles.unitAxes = unitAxes;
guiData.handles.cellScatter = cellScatter;
guiData.handles.guideCellScatter = guideCellScatter;
guiData.handles.brainSliceAxes = brainSliceAxes;
guiData.handles.sliceVolume = 'av'; % The volume shown in the slice
guiData.handles.probeGuideLine = probeGuideLine;
guiData.handles.recordingProbe = recordingProbe;
guiData.handles.histologyPoints = histologyPoints;
guiData.handles.slicePlot = slicePlot; % Slice on 3D atlas
guiData.handles.probeRegionsPlot = probeRegionsPlot;
guiData.handles.corrAxis = corrAxis;
guiData.handles.lfpSpectrumPlot = lfpSpectrumPlot;
guiData.handles.lfpSpectrumLine = lfpSpectrumPlot;
guiData.title = annotation('textbox', [0.25, 0.95, 0.5, 0], 'string', 'My Text', 'EdgeColor', 'none',...
    'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
guiData.probeDetails = annotation('textbox', [0.25, 0.925, 0.5, 0], 'string', 'My Text', 'EdgeColor', 'none',...
    'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
guidata(probeDepthGui, guiData);
%%
switchPenetration(probeDepthGui);
end

function keyPress(probeDepthGui,eventdata)
guiData = guidata(probeDepthGui);
switch eventdata.Key
    case 'uparrow'
        guiData.currTipLocation = guiData.currTipLocation-(guiData.probeStepSize/10).*guiData.currProbeVector';
    case 'downarrow'
        guiData.currTipLocation = guiData.currTipLocation+(guiData.probeStepSize/10).*guiData.currProbeVector';
    case 't'
        if strcmp(guiData.handles.sliceVolume, 'av'); guiData.handles.sliceVolume = 'tv';
        else, guiData.handles.sliceVolume = 'av';
        end
    case 's'
        ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
        ephysRecord(guiData.currPenetration).calcLine = guiData.currProbeVector';
        ephysRecord(guiData.currPenetration).calcTip = guiData.currTipLocation;
        ephysRecord(guiData.currPenetration).scalingFactor = guiData.currScalingFactor;
        save(prc.pathFinder('ephysrecord'), 'ephysRecord');
        set(guiData.title, 'String', strrep(get(guiData.title, 'String'), 'UNSAVED', 'SAVED'));
    case 'c'
        ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
        recordIdx = [ephysRecord.penetrationIdx] == guiData.currPenetration;
        ephysRecord(recordIdx).calcLine = nan;
        ephysRecord(recordIdx).calcTip = nan;
        ephysRecord(recordIdx).scalingFactor = nan;
        save(prc.pathFinder('ephysrecord'), 'ephysRecord');
        set(guiData.title, 'String', strrep(get(guiData.title, 'String'), 'SAVED', 'UNSAVED'));
    case 'equal'
        guiData.penetrations = circshift(guiData.penetrations, -1);
        guiData.currPenetration = guiData.penetrations(1);
    case 'hyphen'
        guiData.penetrations = circshift(guiData.penetrations, 1);
        guiData.currPenetration = guiData.penetrations(1);
    case 'z'
        guiData.currScalingFactor = guiData.currScalingFactor+0.01;
    case 'x'
        guiData.currScalingFactor = guiData.currScalingFactor-0.01;
    case 'p'
        guiData.currPlotIdx = ((guiData.currPlotIdx-1)*-1)+2;
        plotData = guiData.spectrumPlotXYZ(guiData.currPlotIdx,:);
        if guiData.currPlotIdx == 2; set(guiData.handles.corrAxis,'XLim', [0 3840]);
        else, set(guiData.handles.corrAxis,'XLim', [0 200]);
        end
        set(guiData.handles.lfpSpectrumPlot, 'XData', plotData{1}, 'YData', plotData{2}, 'CData', plotData{3});
end
guidata(probeDepthGui, guiData);
end

function keyRelease(probeDepthGui,eventdata)
switch eventdata.Key
    case {'rightarrow','leftarrow','uparrow','downarrow','z','x','s'}
        updateProbePosition(probeDepthGui);
    case {'t'}
        updateSlice(probeDepthGui);
    case {'equal', 'hyphen', 'c'}
        switchPenetration(probeDepthGui);
end
end

function updateSlice(probeDepthGui,varargin)
guiData = guidata(probeDepthGui);
axes(guiData.handles.brainSliceAxes);
currCamPos = campos;
currEntryPoint = guiData.currEntryPoint([3 1 2]);
currProbeVector = guiData.currProbeVector([3 1 2])';

% Get probe-camera vector
probeCameraVector = currEntryPoint - currCamPos;
% Get the vector to plot the plane in (along with probe vector)
planePlotVector = cross(probeCameraVector,currProbeVector);
% Get the normal vector of the plane
normVector = cross(planePlotVector,currProbeVector);
% Get the plane offset through the probe
planeOffset = -(normVector*currEntryPoint');

% Define a plane of points to index (the plane grid is defined based on the which cardinal plan is most
% orthogonal to the plotted plane. this is janky but it works)
slicePixSpace = 3;
[~,camPlane] = max(abs(normVector./norm(normVector)));

switch camPlane
    case 1
        [plnY,plnZ] = meshgrid(1:slicePixSpace:size(guiData.tv,3),1:slicePixSpace:size(guiData.tv,2));
        plnX = (normVector(2)*plnY+normVector(3)*plnZ + planeOffset)/-normVector(1);
    case 2
        [plnX,plnZ] = meshgrid(1:slicePixSpace:size(guiData.tv,1),1:slicePixSpace:size(guiData.tv,2));
        plnY = (normVector(1)*plnX+normVector(3)*plnZ + planeOffset)/-normVector(2);
    case 3
        [plnX,plnY] = meshgrid(1:slicePixSpace:size(guiData.tv,1),1:slicePixSpace:size(guiData.tv,3));
        plnZ = (normVector(1)*plnX+normVector(2)*plnY + planeOffset)/-normVector(3);
end

% Get the coordiates on the plane
xIdx = round(plnX);
yIdx = round(plnY);
zIdx = round(plnZ);
vSize = size(guiData.tv);

% Find plane coordinates in bounds with the volume
useIdx = min(cat(3,xIdx,yIdx,zIdx),[],3)>0 & xIdx<vSize(1) & yIdx<vSize(3) & zIdx<vSize(2);
currSliceIdx = sub2ind(size(guiData.tv),xIdx(useIdx),zIdx(useIdx),yIdx(useIdx));

% Find plane coordinates that contain brain
currSliceIsbrain = false(size(useIdx));
currSliceIsbrain(useIdx) = guiData.av(currSliceIdx) > 1;

% Index coordinates in bounds + with brain
finalPixIdx = sub2ind(size(guiData.tv),xIdx(currSliceIsbrain),zIdx(currSliceIsbrain),yIdx(currSliceIsbrain));

% Grab pixels from (selected) volume
currSlice = nan(size(useIdx));
switch guiData.handles.sliceVolume
    case 'tv'
        currSlice(currSliceIsbrain) = guiData.tv(finalPixIdx);
        colormap(guiData.handles.brainSliceAxes,'gray');
        caxis([0,255]);
    case 'av'
        currSlice(currSliceIsbrain) = guiData.av(finalPixIdx);
        colormap(guiData.handles.brainSliceAxes,guiData.cmap);
        caxis([1,size(guiData.cmap,1)]);
end

% Update the slice display
set(guiData.handles.slicePlot,'XData',plnX,'YData',plnY,'ZData',plnZ,'CData',currSlice);

% Upload guiData
guidata(probeDepthGui, guiData);
end


function switchPenetration(probeDepthGui)
guiData = guidata(probeDepthGui);
expList = guiData.expList;
ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
penetrationDetails = ephysRecord(guiData.currPenetration);

ephysRecordForBrain = ephysRecord(contains({ephysRecord.subject}', penetrationDetails.subject));
exp4Pen = expList(strcmp({expList.subject}', penetrationDetails.subject) & strcmp({expList.expDate}', penetrationDetails.expDate));
if strcmp(penetrationDetails.expNum, 'all'); exp4Pen = exp4Pen(1);
else, exp4Pen = exp4Pen(strcmp({exp4Pen.expNum}', penetrationDetails.expNum));
end

blk = prc.combineBlocks(prc.getDataFromDates(penetrationDetails.subject, penetrationDetails.expDate, exp4Pen.expDef,'eph'));
blk = prc.filtBlock(blk, blk.pen.ephysRecordIdx==guiData.currPenetration, 'penetration');
clusterDepths = blk.clu.depths;
spikeTimes = blk.spk.times;
spikeCluster = blk.spk.clusterNumber;
%%
powerSpectra = flipud(log(blk.pen.lfpPowerSpectra{1}.powerSpectra'));
channelMean = mean(powerSpectra,2);
channels2Remove = abs(channelMean-mean(channelMean))>2;
powerSpectra(channels2Remove,:) = median(powerSpectra(:));
keptMedians = mat2gray(median(powerSpectra(~channels2Remove,:),2));
freqPoints = blk.pen.lfpPowerSpectra{1}.freqPoints;
idx2highlight = sub2ind(size(powerSpectra), find(~channels2Remove),round((max(freqPoints)-1)*keptMedians)+1);
powerSpectra(idx2highlight) = nan;


guiData.spectrumPlotXYZ{1,1} = freqPoints;
guiData.spectrumPlotXYZ{1,2} = 1:10:guiData.probeLength;
guiData.spectrumPlotXYZ{1,3} = powerSpectra;

numCorrGroups = 25;
groupEdgeDepth = linspace(0,guiData.probeLength,numCorrGroups+1);
groupDepths = discretize(clusterDepths,groupEdgeDepth);
groupDepthsCenters = groupEdgeDepth(1:end-1)+(diff(groupEdgeDepth)/2);
uniqueDepths = 1:length(groupEdgeDepth)-1;

spikeBinning = 0.01; % seconds
corrEdges = nanmin(spikeTimes):spikeBinning:nanmax(spikeTimes);
binnedSpikesDepth = zeros(length(uniqueDepths),length(corrEdges)-1);
for currDepth = 1:length(uniqueDepths)
    binnedSpikesDepth(currDepth,:) = histcounts(spikeTimes(ismember(spikeCluster,find(groupDepths == uniqueDepths(currDepth)))), corrEdges);
end
muaCorr = corrcoef(binnedSpikesDepth');
muaCorr(eye(numCorrGroups)>0) = nan;

guiData.spectrumPlotXYZ{2,1} = groupDepthsCenters;
guiData.spectrumPlotXYZ{2,2} = groupDepthsCenters;
guiData.spectrumPlotXYZ{2,3} = muaCorr;

plotData = guiData.spectrumPlotXYZ(guiData.currPlotIdx,:);
set(guiData.handles.lfpSpectrumPlot, 'XData', plotData{1}, 'YData', plotData{2}, 'CData', plotData{3});
%%
spikeJitter = abs(randn(length(clusterDepths), 1)/2.5)+0.1;
set(guiData.handles.cellScatter, 'XData', spikeJitter, 'YData', clusterDepths);
set(guiData.handles.unitAxes, 'xlim', [0 1]);

%%
if ~penetrationDetails.probePainted
    refProbeIdx = find([ephysRecord.probePainted]' & [ephysRecord.histIdx]'==penetrationDetails.histIdx & contains({ephysRecord.subject}', penetrationDetails.subject));
    refProbeDetails = ephysRecord(refProbeIdx);
    exp4Pen = expList(strcmp({expList.subject}', refProbeDetails.subject) & strcmp({expList.expDate}', refProbeDetails.expDate));
    if strcmp(refProbeDetails.expNum, 'all'); exp4Pen = exp4Pen(1);
    else, exp4Pen = exp4Pen(strcmp({exp4Pen.expNum}', refProbeDetails.expNum));
    end
    blk = prc.combineBlocks(prc.getDataFromDates(refProbeDetails.subject, refProbeDetails.expDate, exp4Pen.expDef,'eph'));
    blk = prc.filtBlock(blk, blk.pen.ephysRecordIdx==refProbeIdx, 'penetration');
    clusterDepths = blk.clu.depths;
    spikeJitter = abs(randn(length(clusterDepths), 1)/2.5)+0.1;
    
    guiData.guideTipLocation = refProbeDetails.calcTip;
    guiData.guideCellDepths = clusterDepths;
    guiData.currScalingFactor = refProbeDetails.scalingFactor;
    set(guiData.handles.guideCellScatter, 'XData', spikeJitter, 'YData', clusterDepths);
    set(guiData.title, 'String', sprintf('%s: %s Penetration %d--Probe was NOT painted', penetrationDetails.subject,penetrationDetails.expDate,guiData.currPenetration));
else
    guiData.currScalingFactor = penetrationDetails.scalingFactor;
    if isnan(guiData.currScalingFactor) && isnan(nanmean([ephysRecordForBrain.scalingFactor])); guiData.currScalingFactor = 1; 
    elseif isnan(guiData.currScalingFactor); guiData.currScalingFactor = (nanmean([ephysRecordForBrain.scalingFactor]));
    end
    guiData.guideTipLocation = nan;
    set(guiData.handles.guideCellScatter, 'XData', [], 'YData', []);
    set(guiData.title, 'String', sprintf('%s: %s Penetration %d--Probe was painted', penetrationDetails.subject,penetrationDetails.expDate,guiData.currPenetration));
end

histologyData = load([prc.pathFinder('probepathdata', penetrationDetails.subject) num2str(penetrationDetails.histIdx)]);
histologyPoints = histologyData.probe_ccf;
histologyPoints = [histologyPoints(:,3) histologyPoints(:,2) histologyPoints(:,1)];
probeCenter = mean(histologyPoints,1);
xyz = bsxfun(@minus,histologyPoints,probeCenter);
[~,~,V] = svd(xyz,0);
guiData.currProbeVector = V(:,1);
if guiData.currProbeVector(2) < 0; guiData.currProbeVector = guiData.currProbeVector*-1; end
guidePnts = round(bsxfun(@plus,bsxfun(@times,(-1000:1000)',guiData.currProbeVector'),probeCenter));
guidePnts = unique(guidePnts, 'rows', 'stable');
guidePnts(any(bsxfun(@gt, guidePnts, fliplr(size(guiData.av))),2) | any(guidePnts<=0,2),:) = [];
pntsBrainSpace = sub2ind(size(guiData.av),guidePnts(:,3), guidePnts(:,2), guidePnts(:,1));
pointsInBrain = guidePnts(guiData.av(pntsBrainSpace)>1,:);
guiData.currEntryPoint = pointsInBrain(pointsInBrain(:,2)==min(pointsInBrain(:,2)),:);
guiData.currEntryPoint = guiData.currEntryPoint(1,:);

guiData.currTipLocation = penetrationDetails.calcTip;
if isnan(guiData.currTipLocation)
    guiData.currTipLocation = guiData.currEntryPoint+(penetrationDetails.estDepth/10).*guiData.currProbeVector';
    set(guiData.title, 'String', [get(guiData.title, 'String') ': Probe Location UNSAVED...']);
else, set(guiData.title, 'String', [get(guiData.title, 'String') ': Probe Location SAVED']);
end

set(guiData.handles.probeGuideLine, 'XData',guidePnts(:,3), 'YData',guidePnts(:,1), 'ZData',guidePnts(:,2))
set(guiData.handles.histologyPoints, 'XData',histologyPoints(:,3), 'YData',histologyPoints(:,1), 'ZData',histologyPoints(:,2))

guidata(probeDepthGui, guiData);
updateProbePosition(probeDepthGui);
updateSlice(probeDepthGui);
end


function updateProbePosition(probeDepthGui)
guiData = guidata(probeDepthGui);
%%
if ~isnan(guiData.guideTipLocation)
    tipDistance = sqrt(sum((guiData.currTipLocation-guiData.guideTipLocation).^2))*10;
    tipDistance = tipDistance*sign(guiData.currTipLocation(2)-guiData.guideTipLocation(2));
    set(guiData.handles.guideCellScatter, 'YData',guiData.guideCellDepths-tipDistance);
end

currEndPoint = guiData.currTipLocation;
currStartPoint = currEndPoint-(guiData.probeLength/10).*guiData.currProbeVector'*guiData.currScalingFactor;
pointsAlongProbe = round(cell2mat(arrayfun(@(x,y) linspace(x,y,round(guiData.probeLength/5))', currStartPoint, currEndPoint, 'uni', 0)));

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





