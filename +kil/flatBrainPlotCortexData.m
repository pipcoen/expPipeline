function flatBrainPlotCortexData
%% Flattening the mouse cortex

% TO DO BEFORE RELEASING:
% remove spurious variables and refactor
% turn it into a few functions...

% TO DO ONE DAY (GUIDED BY AN ATLAS)
% add Olfactory areas (but not the bulb) and Entorhinal cortex
% id_entorhinal = StructureTree.id(strcmp(StructureTree.name,'Entorhinal area')); % it's 909
% ii_entorhinal = contains(StructureTree.structure_id_path, ['/' num2str(id_entorhinal) '/']);
% StructureTree(ii_entorhinal,:)
% id_olfactory = StructureTree.id(strcmp(StructureTree.name,'Olfactory areas')); % it's 690
% ii_olfactory = contains(StructureTree.structure_id_path, ['/' num2str(id_olfactory) '/']);
% StructureTree(ii_olfactory,:)
% I would probably drop the olfactory bulb
% notice that I don't need a volume with 5 layers. All I need is layers 1
% and 2 and equivalents in other types of cortex.

%% load volumes and structure tree
allenPath = prc.pathFinder('allenAtlas');

% A table that explains the labels
StructureTree = loadStructureTree([allenPath '/structure_tree_safe_2017.csv']); 
nStructures = size(StructureTree,1);

% Template Volume: gray-scale "background signal intensity" [AP,DV,ML]
% TemplateVolume = readNPY('../Allen/template_volume_10um.npy'); 

% Annotation Volume: the number at each pixel labels the area [AP,DV,ML]
AnnotatedVolume = readNPY([allenPath '/annotation_volume_10um_by_index.npy']); 
nn = size(AnnotatedVolume);

% Allen CCF-bregma transform (from eyeballing Paxinos->CCF)
% this also is a function: AllenCCFBregma
bregma = [540,0,570]; % [AP,DV,ML]

%% Inspect the annotation volume

iX = 870; % a place of high curvature -- used to be bregma(1);
iZ = bregma(3)-200;

ShowSections( AnnotatedVolume, iX, iZ );
colormap(allen_ccf_colormap('2017'));


%% RUN THIS ONLY ONCE - make a volume with layers (from 1 to 5)

% this should work? but it takes a long time
% LayerVolume = layers(AnnotatedVolume);
 
% id_isocortex = StructureTree.id(strcmp(StructureTree.name,'Isocortex')); % the id of isocortex is 315
% ii_isocortex = contains(StructureTree.structure_id_path, ['/' num2str(id_isocortex) '/']);
% 
% LayerNames = {'layer 1','layer 2','layer 4','layer 5','layer 6'};
% 
% LayerVolume = zeros(size(AnnotatedVolume),'uint8');
% for iLayer = 1:5
%     ii_layer = contains(StructureTree.name,LayerNames{iLayer},'IgnoreCase',true)&ii_isocortex;
%     LayerRows = find(ii_layer); % rows of st that have this layer
%     for iLayerRow = 1:length(LayerRows)
%         iRow = LayerRows(iLayerRow);
%         disp(StructureTree.name(iRow));
%         LayerVolume(AnnotatedVolume==iRow)=iLayer;
%     end
% end
% 
% save('Data/LayerVolume','LayerVolume');

%% Load LayerVolume

load('Data/LayerVolume');

% hack: put a 6 in the non-cortex parts of the brain
LayerVolume = LayerVolume + uint8(AnnotatedVolume>1);
LayerVolume(LayerVolume==1) = uint8(7);
LayerVolume = LayerVolume -1;

ShowSections( LayerVolume, iX, iZ );

%% Assign a cortical area to every point (zero if not in cortex)

% I edited STructureTree to have cortical areas only
% save CortexTree CortexTree
load('Data/CortexTree');

LayerRows = unique(AnnotatedVolume(LayerVolume<6&LayerVolume>0));
% these are rows of the table StructureTree
% they correspond to index + 1
% StructureTree(LayerRows,:)

areas = zeros(nStructures,1);
for iRow = LayerRows' % iRow = 8 is one of them
    if any( CortexTree.index == iRow-1 )
        % ind is in CortexTree already
        areas(iRow) = iRow;
    else
        NewRow = iRow;
        while(~any(CortexTree.index==NewRow-1))
            NewRow = find(StructureTree.id==StructureTree.parent_structure_id(NewRow));
        end
        areas(iRow) = NewRow;
    end
end

iiCorticalAreas = unique(areas(areas>0));
nCorticalAreas = length(iiCorticalAreas);

CorticalAreas = StructureTree(iiCorticalAreas,:);

%% Downsample it and smooth it

DecFac = 5; % 10 can run on my laptop (32 GB), 5 requires a workstation (128 GB)

% decimate it 
DecLayerVolume = double(LayerVolume(1:DecFac:nn(1),1:DecFac:nn(2),1:DecFac:nn(3)));
dnn = nn/DecFac;

%% find all the edges of layer 2 with layer 1 or with the outside 

CortexSurf = false(dnn);
for ind = find(DecLayerVolume == 2)'
    [i,j,k]=ind2sub(dnn,ind);
    Comparisons = DecLayerVolume(i-1:i+1, j-1:j+1, k-1:k+1)<2;
    CortexSurf(ind)=any(Comparisons(:)); 
end
clear ind i j k

ax = ShowSections( CortexSurf, iX, iZ, DecFac );
colormap(allen_ccf_colormap('2017'));

% drop one hemisphere
CortexSurf(:,:,(dnn(3)/2+1):end) = false;

%% 3D graphics using isosurface

fv = isosurface(double(CortexSurf),0.5);

figure; clf
p = patch( fv);
p.FaceColor = 'red';
p.EdgeColor = 'none';
clear p
lighting gouraud
axis equal

view(3); 
axis tight
camlight headlight

%% Find coordinates of points in CortexSurf

nSurfPoints = nnz(CortexSurf); 
SurfXYZ = zeros(3,nSurfPoints);
[SurfXYZ(1,:),SurfXYZ(2,:),SurfXYZ(3,:)] = ind2sub( dnn, find(CortexSurf) );

%% Assign an area to every point on the surface 

% the rows in StructureTree for each entry in SurfXYZ
SurfRows = areas(AnnotatedVolume(sub2ind(nn,...
    SurfXYZ(1,:)*DecFac,...
    SurfXYZ(2,:)*DecFac,...
    SurfXYZ(3,:)*DecFac)));

SurfAreas = zeros(nSurfPoints,1); 
for iA = 1:nCorticalAreas
    SurfAreas( SurfRows==iiCorticalAreas(iA) ) = iA;
end
SurfAreas(SurfAreas==0) = nCorticalAreas+1; % so the colormap works well
clear SurfRows

%% Load Pip's colormap

% addpath('Data');
MatteoColormap; % loads cMap
PipColorMap = [ cell2mat(cMap(:,2)); [255 255 255]]/255;
clear cMap

%% Plot the surface in 3D with areas colored appropriately

figure; clf
scatter3(SurfXYZ(1,:),SurfXYZ(2,:),SurfXYZ(3,:),20, SurfAreas, 'filled'); hold on
axis equal
axis off
colormap(PipColorMap);

%% Quickest paths to surface

figure; clf
s = scatter3(SurfXYZ(1,:),SurfXYZ(2,:),SurfXYZ(3,:),20, SurfAreas, 'filled'); hold on
axis equal
axis off
colormap(PipColorMap);
set(s, 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.2)

% ShowSections( DecLayerVolume, iX, iZ, DecFac );

% paths to the surface
StartPoints = [[37,37,iZ/10];[41,33,iZ/10];[29,29,iZ/10];[37,45,iZ/10]]*10/DecFac; % or you could choose all the points in L6

for iStartPoint = 1:4  
    start = StartPoints(iStartPoint,:);
    % Find the closest SurfXYZ point
    [~,iSurfPoint] = min( vecnorm(SurfXYZ' - repmat(start,[nSurfPoints, 1]),2,2) );
    plot3([start(1) SurfXYZ(1,iSurfPoint)],...
        [start(2) SurfXYZ(2,iSurfPoint)],...
        [start(3) SurfXYZ(3,iSurfPoint)],'k-','LineWidth',3);
end
            
title('Paths to the surface');

%% Isomap to flatten it -- RUN THIS ONCE

% save X X % I guess X was SurfXYZ
% 
% % Run this on ZULTRA
% load X
% addpath('Isomap');
% [Xstrain,XGraph] = Run3dIsomap(X);
% save Xstrain50 Xstrain
% save XGraph50 XGraph

%% OR BETTER LOAD THIS

load('Data/Xstrain50'); % loads XStrain
load('Data/XGraph50'); % loads XGraph

SurfXY = Xstrain; clear Xstrain

%% calculate the 2D medoid of each area and add it as a field in CorticalAreas

Medoid2D = zeros(nCorticalAreas,2);
for iA = 1:nCorticalAreas
    pp = (SurfAreas ==  iA);
    [~,Medoid2D(iA,:)] = kmedoids( SurfXY(:,pp)', 1 );
end
CorticalAreas = addvars(CorticalAreas,Medoid2D);

%% Plot the flattened surface in 2D with areas colored appropriately

figure; clf
p = plot(XGraph); hold on;
p.XData = SurfXY(1,:);
p.YData = SurfXY(2,:);
p.EdgeColor = 'k';

scatter(SurfXY(1,:),SurfXY(2,:),20,SurfAreas, 'filled'); hold on
axis('equal'); axis tight; axis off
colormap(PipColorMap);

for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), CorticalAreas.acronym(iA),'hori','center','vert','middle', 'color', 'k' ); 
end
title('Flattened projection of the 3D points -- as a graph');
clear p

%% Now do a Delaunay triangulation on the 2-D data

TriangIn2D = delaunayTriangulation( SurfXY(1,:)',  SurfXY(2,:)' );

figure; clf;
triplot(TriangIn2D);
axis equal; axis tight; axis off
title('Delaunay triangulation in 2D');

%% Look at the triangulation in 3-D

% plot the triangulation in 3D coords

TriangIn3D = triangulation(TriangIn2D.ConnectivityList,SurfXYZ');

figure;
trisurf(TriangIn3D, 'FaceColor','cyan','FaceAlpha',0.8);
axis equal
camlight 
lighting gouraud
title('Delaunay triangulation in 3D');

%% The boundary of the region

BoundaryPoly = polyshape(SurfXY(:, boundary( SurfXY', 1))');
% in  world coordinates

%% Now define a 2D cartesian grid and ensure it is inside the cortex

CartScale = 2; % used to be 4;

ni = CartScale * ceil(range(SurfXY(2,:)));
nj = CartScale * ceil(range(SurfXY(1,:)));

WorldLims = zeros(2,2);
[WorldLims(1,1),WorldLims(1,2)] =  bounds(SurfXY(1,:));
[WorldLims(2,1),WorldLims(2,2)] =  bounds(SurfXY(2,:));

Ref2D = imref2d( [ni nj], WorldLims(1,:), WorldLims(2,:) );

% [foo1, foo2] = Ref2D.intrinsicToWorld( 1,1 )
% [foo1, foo2] = Ref2D.intrinsicToWorld( nj,ni )
% [foo1, foo2] = Ref2D.intrinsicToWorld( 1, ni )
% [foo1, foo2] = Ref2D.intrinsicToWorld( nj, 1 )

ii = 1:Ref2D.ImageSize(1); % floor(min(SurfXY(1,:))):1/CartScale:ceil(max(SurfXY(1,:)));
jj = 1:Ref2D.ImageSize(2); % floor(min(SurfXY(2,:))):1/CartScale:ceil(max(SurfXY(2,:)));
% ni = length(ii); 
% nj = length(jj);

[iii,jjj] = meshgrid(1:Ref2D.ImageSize(1),1:Ref2D.ImageSize(2));
CombsJI = [jjj(:), iii(:)]; % all combinations of [i j]

Cart2D = 0*CombsJI;
[Cart2D(:,1),Cart2D(:,2)] = Ref2D.intrinsicToWorld( CombsJI(:,1), CombsJI(:,2) );

% keep only the cartesian points that are in the cortex
Cart2D = Cart2D( isinterior(BoundaryPoly,Cart2D),:);
nCart2D = size(Cart2D,1);

figure; clf
plot(BoundaryPoly); hold on 
% plot( Cart2D(:,1), Cart2D(:,2), 'r.' );
title('Cartesian grid in 2D');
axis equal

%% Assign barycentric coordinates to the cartesian grid

% assign each point to a triangle
Cart2DTriangs = pointLocation(TriangIn2D,Cart2D);

% barycentric coordinates of each point
Bary2D = cartesianToBarycentric(TriangIn2D,Cart2DTriangs,Cart2D);

% figure; clf;
% triplot(TriangIn2D); hold on 
% axis equal; axis tight; axis off
% plot( Cart2D(:,1), Cart2D(:,2), 'k.' )

%% now plot the grid in 3D using their barycentric coords

Cart3D = barycentricToCartesian(TriangIn3D,Cart2DTriangs,Bary2D);
Cart3D = round(Cart3D*DecFac); % scale back up to full resolution

CartRows = areas(AnnotatedVolume(sub2ind(nn,Cart3D(:,1),Cart3D(:,2),Cart3D(:,3))));
CartAreas = CartRows;
for iA = 1:nCorticalAreas
    CartAreas( CartRows==iiCorticalAreas(iA) ) = iA;
end
CartAreas(CartAreas==0) = nCorticalAreas+1; % for graphical purposes

figure; clf
scatter3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),5,CartAreas);
colormap(PipColorMap);
axis equal
title('Cartesian 2D grid plotted in 3D');

%% Compare the areas of the cortical areas in 2D vs 3D

SizeIn2D = zeros(nCorticalAreas,1);
SizeIn3D = zeros(nCorticalAreas,1);
for iA = 1:nCorticalAreas
    SizeIn3D(iA) = nnz(SurfAreas ==  iA);
    SizeIn2D(iA) = nnz(CartAreas ==  iA);
end

figure; plot( SizeIn2D, SizeIn3D, 'ko' );
xlabel('Size in 2D');
ylabel('Size in 3D');

%% Plot a nice flattened map of the cortex, with labels

figure; clf
scatter(Cart2D(:,1),Cart2D(:,2),20, CartAreas, 'filled');
colormap(PipColorMap);
axis equal
title('Cartesian 2D grid plotted in 2D');
for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), ...
        CorticalAreas.acronym(iA),'hori','center','vert','middle', 'color', 'k' ); 
end

%% Make a 2D image of the cortex

[subs1,subs2] = Ref2D.worldToSubscript( Cart2D(:,1), Cart2D(:,2) );

AreaImg = zeros(Ref2D.ImageSize); 
AreaImg( sub2ind(Ref2D.ImageSize,subs1,subs2 )) = CartAreas;
AreaImg(AreaImg==0)= nCorticalAreas+1; % to work well with the colormap

figure; clf; 
image(Ref2D.XWorldLimits,Ref2D.YWorldLimits, AreaImg); 
axis equal; axis tight; axis xy
colormap(PipColorMap); 
title('2D map plotted as an image');



%% Clean up the noise and make a nice 2D map of the areas

CleanImg = AreaImg;
for iA =  [21, 26, 1:nCorticalAreas, 24]
    bw = logical(CleanImg == iA);
    CleanImg(bw) = 0;  % unassigned
    bw = imclose( bw, strel('disk',5));
    CleanImg(bw) = iA;
end

CleanImg(CleanImg==0)= nCorticalAreas+1;

figure; clf; 
image( Ref2D.XWorldLimits,Ref2D.YWorldLimits, CleanImg); 
axis equal; axis tight; axis xy
colormap(PipColorMap); 
title('Smoothed 2D map plotted as an image');

% now these are in the wrong place...

for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), [CorticalAreas.acronym(iA)],'hori','center','vert','middle', 'color', 'k' ); 
end

% safety check
% hold on; plot(BoundaryPoly); 
clear AnnotatedVolume LayerVolume;
save('D:\Dropbox (Neuropixels)\MouseData\matteoBrainData');


