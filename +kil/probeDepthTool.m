function probeDepthTool(subjects)
%%
if ~iscell(subjects); subjects = {subjects}; end
allenPath = prc.pathFinder('allenAtlas');
cmap = load([fileparts(which('allenCCFbregma')) '\allen_ccf_colormap_2017.mat']);
expList = load(prc.pathFinder('expList')); expList = expList.expList;
%%
expList = expList(contains({expList.subject}', subjects) & [expList.excluded]' ~=1 & contains({expList.expType}', 'eph'));
expList = prc.updatePaths(expList);
ephysListPath = prc.pathFinder('ephysrecord');

guiData = struct;
guiData.expList = expList;
guiData.ephysRec = table2struct(readtable(ephysListPath));
guiData.penetrations = unique(cell2mat(cellfun(@(x) reshape([x.penetrationIdx],[],1), {expList.expDets}', 'uni', 0)));
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

% Set up the atlas axes
brainSliceAxes = subplot(plotR,plotC,[1:plotR, (1:plotR)+plotC, (1:plotR)+plotC*2]); hold on;
slicePlot = surface('EdgeColor','none'); % Slice on 3D atlas
[~, guiData.handles.brainOutline] = plotBrainGrid([],brainSliceAxes);
hold(brainSliceAxes,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[apMax,dvMax,mlMax] = size(guiData.tv);
xlim([-10,apMax+10]);
ylim([-10,mlMax+10]);
zlim([-10,dvMax+10])
probeGuideLine = plot3(rand(1), rand(1), rand(1), '--b', 'linewidth', 2.5);
histologyPoints = plot3(rand(1), rand(1), rand(1), '.b','MarkerSize',20);
recordingProbe = plot3(rand(1), rand(1), rand(1), '-m', 'linewidth', 4);

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
title(regionsAxes,'Probe areas');

% Set up units axis
unitAxes = subplot(plotR,plotC,5:plotC:plotR*plotC);
hold(unitAxes,'on');
set(unitAxes,'YDir','reverse', 'yAxisLocation', 'right', 'xcolor', 'w','XTick','','YLim',[0,3840],'YColor','w', 'color', 'none','YTick','');
cellDensityLine = plot(unitAxes,1,1,'-r', 'linewidth', 3);
xlim([-0.1,1]);
ylabel('Depth (\mum)')

% CorrAxis
corrAxis = subplot(plotR,plotC,6:plotC:plotR*plotC);
hold(corrAxis,'on');
set(corrAxis,'XTick','','YLim',[0,3840],'xlim',[0,3840], 'YTick', '');
lfpSpectrumPlot = image(0);
% set(CorrAxis,'YDir','reverse', 'yAxisLocation', 'right', 'xcolor', 'w','XTick','','YLim',[0,3840],'YColor','w', 'color', 'none','YTick','');

linkaxes([unitAxes,regionsAxes],'y');

% Position the axes
set(brainSliceAxes,'Position',[0.1,0.1,0.5,0.9]);
set(regionsAxes,'Position',[0.7,0.1,0.03,0.8]);
set(unitAxes,'Position',[0.7,0.1,0.03,0.8]);
set(corrAxis,'Position',[0.8,0.1,0.12,0.8]);

% Probe information
guiData.currPenetration = guiData.penetrations(1);
guiData.currProbeDepth = [];
guiData.currProbeVector = [];
guiData.currEntryPoint = [];

%set up all handles
guiData.handles.regionsAxes = regionsAxes;
guiData.handles.unitAxes = unitAxes;
guiData.handles.cellDensityLine = cellDensityLine;
guiData.handles.brainSliceAxes = brainSliceAxes;
guiData.handles.sliceVolume = 'av'; % The volume shown in the slice
guiData.handles.probeGuideLine = probeGuideLine;
guiData.handles.recordingProbe = recordingProbe;
guiData.handles.histologyPoints = histologyPoints;
guiData.handles.slicePlot = slicePlot; % Slice on 3D atlas
guiData.handles.probeRegionsPlot = probeRegionsPlot;
guiData.handles.corrAxis = corrAxis;
guiData.handles.lfpSpectrumPlot = lfpSpectrumPlot;
guidata(probeDepthGui, guiData);
%%
switchPenetration(probeDepthGui);
end

function keyPress(probeDepthGui,eventdata)
guiData = guidata(probeDepthGui);
switch eventdata.Key
    case 'uparrow'
        guiData.currProbeDepth = guiData.currProbeDepth-20;
    case 'downarrow'
        guiData.currProbeDepth = guiData.currProbeDepth+20;
    case 't'
        if strcmp(guiData.handles.sliceVolume, 'av'); guiData.handles.sliceVolume = 'tv';
        else, guiData.handles.sliceVolume = 'av';
        end
    case 'n'
        guiData.currPenetration = guiData.penetrations(find(guiData.penetrations==guiData.currPenetration)+1);
        disp( guiData.currPenetration);
end
guidata(probeDepthGui, guiData);
end

function keyRelease(probeDepthGui,eventdata)
switch eventdata.Key
    case {'rightarrow','leftarrow','uparrow','downarrow'}
        updateProbePosition(probeDepthGui);
    case {'t'}
        updateSlice(probeDepthGui);
    case {'n'}
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

keyExperiments = find(cellfun(@(x) any([x.penetrationIdx]==guiData.currPenetration), {guiData.expList.expDets}'));
ephInfo = load(guiData.expList(keyExperiments(1)).processedData, 'eph');
ephInfo = prc.filtStruct(ephInfo.eph, ephInfo.eph.clusterSite==guiData.currPenetration);
ephInfo = prc.filtStruct(ephInfo, ephInfo.spikeSite==guiData.currPenetration);
clusterDepths = ephInfo.clusterDepths;
%%

%%
% numCorrGroups = 20;
% groupEdgeDepth = linspace(0,guiData.probeLength,numCorrGroups+1);
% groupDepths = discretize(clusterDepths,groupEdgeDepth);
% groupDepthsCenters = groupEdgeDepth(1:end-1)+(diff(groupEdgeDepth)/2);
% uniqueDepths = 1:length(groupEdgeDepth)-1;
% 
% spikeBinning = 0.01; % seconds
% corrEdges = nanmin(ephInfo.spikeTimes):spikeBinning:nanmax(ephInfo.spikeTimes);
% binnedSpikesDepth = zeros(length(uniqueDepths),length(corrEdges)-1);
% for currDepth = 1:length(uniqueDepths)
%     binnedSpikesDepth(currDepth,:) = histcounts(ephInfo.spikeTimes(ismember(ephInfo.spikeCluster,find(groupDepths == uniqueDepths(currDepth)))), corrEdges);
% end
% muaCorr = corrcoef(binnedSpikesDepth');
% muaCorr(eye(numCorrGroups)>0) = nan;
% imagesc(guiData.handles.corrAxis, groupDepthsCenters, groupDepthsCenters, muaCorr);
% set(guiData.handles.corrAxis, 'YDir', 'reverse')

%%
set(guiData.handles.lfpSpectrumPlot, 'XData', 'YData', 'Z, 
%%
[cellNumbersAlongProbe, binEdges] = histcounts(clusterDepths, 0:100:guiData.probeLength);
binCenters = mean([binEdges(1:end-1); binEdges(2:end)]);
queryPoints = 0:10:guiData.probeLength;
cellNumbersAlongProbe = interp1(binCenters,cellNumbersAlongProbe,queryPoints,'linear', 'extrap');
set(guiData.handles.cellDensityLine, 'XData', cellNumbersAlongProbe, 'YData', queryPoints);
set(guiData.handles.unitAxes, 'xlim', [0 max(cellNumbersAlongProbe)+2]);

penetrationDetails = guiData.ephysRec([guiData.ephysRec.penetrationIdx]==guiData.currPenetration);
histologyData = load(prc.pathFinder('probepathdata', penetrationDetails.subject));
histologyPoints = histologyData.pointList.pointList{penetrationDetails.histIdx,1};
guiData.currProbeDepth = penetrationDetails.calcDepth;
if isnan(guiData.currProbeDepth); guiData.currProbeDepth = penetrationDetails.estDepth; end

probeCenter = mean(histologyPoints,1);
xyz = bsxfun(@minus,histologyPoints,probeCenter);
[~,~,V] = svd(xyz,0);
guiData.currProbeVector = V(:,1);
if guiData.currProbeVector(2) < 0; guiData.currProbeVector = guiData.currProbeVector*-1; end

guidePnts = round(bsxfun(@plus,bsxfun(@times,(-1000:1000)',guiData.currProbeVector'),probeCenter));
guidePnts = unique(guidePnts, 'rows', 'stable');
guidePnts(any(bsxfun(@gt, guidePnts, size(guiData.av)),2) | any(guidePnts<=0,2),:) = [];
pntsBrainSpace = sub2ind(size(guiData.av),guidePnts(:,3), guidePnts(:,2), guidePnts(:,1));
pointsInBrain = guidePnts(guiData.av(pntsBrainSpace)>1,:);
guiData.currEntryPoint = pointsInBrain(pointsInBrain(:,2)==min(pointsInBrain(:,2)),:);

set(guiData.handles.probeGuideLine, 'XData',guidePnts(:,3), 'YData',guidePnts(:,1), 'ZData',guidePnts(:,2))
set(guiData.handles.histologyPoints, 'XData',histologyPoints(:,3), 'YData',histologyPoints(:,1), 'ZData',histologyPoints(:,2))

guidata(probeDepthGui, guiData);
updateProbePosition(probeDepthGui);
updateSlice(probeDepthGui);
end


function updateProbePosition(probeDepthGui)
guiData = guidata(probeDepthGui);
%%
currEndPoint = guiData.currEntryPoint+(guiData.currProbeDepth/10).*guiData.currProbeVector';
currStartPoint = currEndPoint-(guiData.probeLength/10).*guiData.currProbeVector';
recordingProbe = [currEndPoint; currStartPoint];

set(guiData.handles.recordingProbe, 'XData',recordingProbe(:,3), 'YData',recordingProbe(:,1), 'ZData',recordingProbe(:,2))
pointsAlongProbe = round(cell2mat(arrayfun(@(x,y) linspace(x,y,round(guiData.probeLength/5))', currStartPoint, currEndPoint, 'uni', 0)));
pointsAlongProbe(any(bsxfun(@gt, pointsAlongProbe, size(guiData.av)),2) | any(pointsAlongProbe<=0,2),:) = 1;
probeAreas = arrayfun(@(x,y,z) guiData.av(x,y,z), pointsAlongProbe(:,3),pointsAlongProbe(:,2),pointsAlongProbe(:,1));
probeAreaBoundaries = find(diff([1e5; double(probeAreas); 1e5])~=0);
probeAreaCenters = (probeAreaBoundaries(2:end) + probeAreaBoundaries(1:end-1))/2;
probeAreaLabels = guiData.st.safe_name(probeAreas(round(probeAreaCenters)));
[~, ~, uniLocations] = unique(probeAreaLabels, 'stable');
newAreaBoundaries = find(diff([-1; uniLocations;-1])~=0);
probeAreaLabels = probeAreaLabels(newAreaBoundaries(1:end-1));
for i = 1:length(newAreaBoundaries)-1
    idx2Avergage = newAreaBoundaries(i):newAreaBoundaries(i+1)-1;
    probeAreaCenters(idx2Avergage) = mean(probeAreaCenters(idx2Avergage));
end
probeAreaCenters = probeAreaCenters(newAreaBoundaries(1:end-1));

% Update the probe areas
yyaxis(guiData.handles.regionsAxes,'right');
set(guiData.handles.probeRegionsPlot,'YData',(1:length(probeAreas))*5,'CData',probeAreas);
set(guiData.handles.regionsAxes,'YTick',probeAreaCenters*5,'YTickLabels',probeAreaLabels);

% Upload guiData
guidata(probeDepthGui, guiData);
end

















