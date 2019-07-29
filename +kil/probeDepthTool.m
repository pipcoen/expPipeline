function probeDepthTool(subjects)
%%
if ~iscell(subjects); subjects = {subjects}; end
allenPath = prc.pathFinder('allenAtlas');
cmap = load([fileparts(which('allenCCFbregma')) '\allen_ccf_colormap_2017.mat']);
expList = load(prc.pathFinder('expList')); expList = expList.expList;
%%
expList = expList(contains({expList.subject}', subjects) & [expList.excluded]' ~=1 & contains({expList.expType}', 'eph'));
expList = prc.updatePaths(expList);

guiData = struct;
guiData.sessionData = cellfun(@load, {expList.processedData}', 'uni', 0);
guiData.subjects = subjects;
guiData.probeDataPaths = cellfun(@(x) prc.pathFinder('probepathdata', x), subjects, 'uni', 0);
guiData.tv = readNPY([allenPath 'template_volume_10um.npy']); % grey-scale "background signal intensity"
guiData.av = readNPY([allenPath 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
guiData.st = loadStructureTree([allenPath 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
guiData.cmap = cmap.cmap; % Atlas colormap
guiData.bregma = [540,0,570];
guiData.probeLength = 3.84; % Length of probe
guiData.plottedStructuresIdx = []; % Plotted structures
guiData.probeAngle = [0;90]; % Probe angles in ML/DV
%%
probeDepthGui = figure('color','w');

% Set up the atlas axes
guiData.handles.brainSliceAxis = subplot(3,5,[1:3, 6:8, 11:12]); hold on;
[~, guiData.handles.brainOutline] = plotBrainGrid([],guiData.handles.brainSliceAxis);
hold(guiData.handles.brainSliceAxis,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[apMax,dvMax,mlMax] = size(guiData.tv);
xlim([-10,apMax+10]);
ylim([-10,mlMax+10]);
zlim([-10,dvMax+10])
guiData.handles.sliceVolume = 'tv'; % The volume shown in the slice

% Make 3D rotation the default state (toggle on/off with 'r')
h = rotate3d(guiData.handles.brainSliceAxis);
h.Enable = 'on';
% Update the slice whenever a rotation is completed
h.ActionPostCallback = @updateSlice;
% axis off;

% Set up the probe area axes
guiData.handles.regionsAxis = subplot(3,5,4:5:15);
guiData.handles.regionsAxis.ActivePositionProperty = 'position';
set(guiData.handles.regionsAxis,'FontSize',11);
yyaxis(guiData.handles.regionsAxis,'left');
guiData.probeRegionsPlot = image(0);
set(guiData.handles.regionsAxis,'XTick','','YLim',[0,3840],'YColor','k','YDir','reverse');
ylabel(guiData.handles.regionsAxis,'Depth (\mum)');
colormap(guiData.handles.regionsAxis,guiData.cmap);
caxis([1,size(guiData.cmap,1)])
yyaxis(guiData.handles.regionsAxis,'right');
set(guiData.handles.regionsAxis,'XTick','','YLim',[0,3840],'YColor','k','YDir','reverse');
title(guiData.handles.regionsAxis,'Probe areas');

% Set up units axis
guiData.handles.unitAxes = subplot(3,5,5:5:15);
set(guiData.handles.unitAxes,'YDir','reverse', 'yAxisLocation', 'right', 'xcolor', 'w','XTick','','YLim',[0,3840],'YColor','k');
xlim([-0.1,1]);
ylabel('Depth (\mum)')


% Position the axes
set(guiData.handles.brainSliceAxis,'Position',[0.1,0.1,0.5,0.9]);
set(guiData.handles.regionsAxis,'Position',[0.7,0.1,0.03,0.8]);

% Probe information
guiData.currSuject = 1;
guiData.currProbe = 1;
guiData.currDeopthOffset = 0;

switchProbe(probeDepthGui);
% updateSlice(probeDepthGui);
guidata(probeDepthGui, guiData);
end

function updateSlice(probeDepthGui,varargin)
% Get guidata
guiData = guidata(probeDepthGui);
% Only update the slice if it's visible
% Get current position of camera
currCamPos = campos;

    % Get probe vector
    probe_ref_top = [guiData.handles.probe_ref_line.XData(1), ...
        guiData.handles.probe_ref_line.YData(1),guiData.handles.probe_ref_line.ZData(1)];
    probe_ref_bottom = [guiData.handles.probe_ref_line.XData(2), ...
        guiData.handles.probe_ref_line.YData(2),guiData.handles.probe_ref_line.ZData(2)];
    probe_vector = probe_ref_top - probe_ref_bottom;
    
    % Get probe-camera vector
    probe_camera_vector = probe_ref_top - currCamPos;
    
    % Get the vector to plot the plane in (along with probe vector)
    plot_vector = cross(probe_camera_vector,probe_vector);
    
    % Get the normal vector of the plane
    normal_vector = cross(plot_vector,probe_vector);
    
    % Get the plane offset through the probe
    plane_offset = -(normal_vector*probe_ref_top');
    
    % Define a plane of points to index
    % (the plane grid is defined based on the which cardinal plan is most
    % orthogonal to the plotted plane. this is janky but it works)
    slice_px_space = 3;
    %[~,cam_plane] = max(abs((campos - camtarget)./norm(campos - camtarget)));
    
    [~,cam_plane] = max(abs(normal_vector./norm(normal_vector)));
    
    switch cam_plane
        
        case 1
            [plane_y,plane_z] = meshgrid(1:slice_px_space:size(guiData.tv,3),1:slice_px_space:size(guiData.tv,2));
            plane_x = ...
                (normal_vector(2)*plane_y+normal_vector(3)*plane_z + plane_offset)/ ...
                -normal_vector(1);
            
        case 2
            [plane_x,plane_z] = meshgrid(1:slice_px_space:size(guiData.tv,1),1:slice_px_space:size(guiData.tv,2));
            plane_y = ...
                (normal_vector(1)*plane_x+normal_vector(3)*plane_z + plane_offset)/ ...
                -normal_vector(2);
            
        case 3
            [plane_x,plane_y] = meshgrid(1:slice_px_space:size(guiData.tv,1),1:slice_px_space:size(guiData.tv,3));
            plane_z = ...
                (normal_vector(1)*plane_x+normal_vector(2)*plane_y + plane_offset)/ ...
                -normal_vector(3);
            
    end
    
    % Get the coordiates on the plane
    x_idx = round(plane_x);
    y_idx = round(plane_y);
    z_idx = round(plane_z);
    
    % Find plane coordinates in bounds with the volume
    use_xd = x_idx > 0 & x_idx < size(guiData.tv,1);
    use_yd = y_idx > 0 & y_idx < size(guiData.tv,3);
    use_zd = z_idx > 0 & z_idx < size(guiData.tv,2);
    use_idx = use_xd & use_yd & use_zd;
    
    curr_slice_idx = sub2ind(size(guiData.tv),x_idx(use_idx),z_idx(use_idx),y_idx(use_idx));
    
    % Find plane coordinates that contain brain
    curr_slice_isbrain = false(size(use_idx));
    curr_slice_isbrain(use_idx) = guiData.av(curr_slice_idx) > 1;
    
    % Index coordinates in bounds + with brain
    grab_pix_idx = sub2ind(size(guiData.tv),x_idx(curr_slice_isbrain),z_idx(curr_slice_isbrain),y_idx(curr_slice_isbrain));
    
    % Grab pixels from (selected) volume
    curr_slice = nan(size(use_idx));
    switch guiData.handles.slice_volume
        case 'tv'
            curr_slice(curr_slice_isbrain) = guiData.tv(grab_pix_idx);
            colormap(guiData.handles.axes_atlas,'gray');
            caxis([0,255]);
        case 'av'
            curr_slice(curr_slice_isbrain) = guiData.av(grab_pix_idx);
            colormap(guiData.handles.axes_atlas,guiData.cmap);
            caxis([1,size(guiData.cmap,1)]);
    end
    
    % Update the slice display
    set(guiData.handles.slice_plot,'XData',plane_x,'YData',plane_y,'ZData',plane_z,'CData',curr_slice);
    
    % Upload guiData
    guidata(probeDepthGui, guiData);
end

function switchProbe(probeDepthGui)
% Load histology points from P Shamash's program and draw best-fit line
% (AP_cortexlab_filename option: 'probe_histology')
guiData = guidata(probeDepthGui);
%%
guiData.currProbeData = load(guiData.probeDataPaths{guiData.currBrain});
guiData.currProbeData = guiData.currProbeData.pointList.pointList;
histologyPoints = guiData.currProbeData{guiData.currProbe,1};

probeCenter = mean(histologyPoints,1);
xyz = bsxfun(@minus,histologyPoints,probeCenter);
[~,~,V] = svd(xyz,0);
histologyProbeDirection = V(:,1);
pnts = round(bsxfun(@plus,bsxfun(@times,(-1000:1000)',histologyProbeDirection'),probeCenter));
% pnts = pnts(:, [3 2 1]);
pnts = unique(pnts, 'rows', 'stable');
pnts(any(bsxfun(@gt, pnts, size(guiData.av)),2) | any(pnts<=0,2),:) = [];


pntsBrainSpace = sub2ind(size(guiData.av),pnts(:,3), pnts(:,2), pnts(:,1));
pointsInBrain = pnts(guiData.av(pntsBrainSpace)>1,:);
% pointsInBrain = pnts;
%% 
plot3(histologyPoints(:,3),histologyPoints(:,1),histologyPoints(:,2),'.b','MarkerSize',20);
plot3(pointsInBrain(:,3), pointsInBrain(:,1), pointsInBrain(:,2), '*')
line(pntsBrainSpace(:,3),pntsBrainSpace(:,1),pntsBrainSpace(:,2),'color','k--','linewidth',2)
%%
guidata(probeDepthGui, guiData);
end


















