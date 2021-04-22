
%% Flattening the mouse cortex

% To do:
% % make curves for area borders and project them to 3D and 2D
% come up with color scheme for 2D cortex
% (segmentation of regions could be done with superpixel3 maybe?)

%% load volumes and structure tree

addpath(genpath('../allenCCF'));
addpath(genpath('../npy-matlab'));

% A table that explains the labels
StructureTree = loadStructureTree('../Allen/structure_tree_safe_2017.csv'); 
nStructures = size(StructureTree,1);

% Template Volume: gray-scale "background signal intensity" [AP,DV,ML]
% TemplateVolume = readNPY('../Allen/template_volume_10um.npy'); 

% Annotation Volume: the number at each pixel labels the area [AP,DV,ML]
AnnotatedVolume = readNPY('../Allen/annotation_volume_10um_by_index.npy'); 
nn = size(AnnotatedVolume);

% Allen CCF-bregma transform (from eyeballing Paxinos->CCF)
% this also is a function: AllenCCFBregma
bregma = [540,0,570]; % [AP,DV,ML]

%% Inspect the annotation volume

iX = 870; % a place of high curvature -- used to be bregma(1);
iZ = bregma(3)-200;

ShowSections( AnnotatedVolume, iX, iZ );
colormap(allen_ccf_colormap('2017'));

%% Assign a layer to every structure (NaN if not in cortex)

% THiS STUFF IS NOT NEEDED, IS IT?

% id_isocortex = StructureTree.id(strcmp(StructureTree.name,'Isocortex')); % the id of isocortex is 315
% ii_isocortex = contains(StructureTree.structure_id_path, ['/' num2str(id_isocortex) '/']);
% 
% layers = nan(nStructures,1);
% LayerNames = {'layer 1','layer 2','layer 4','layer 5','layer 6'};
% for iLayer = 1:5
%     ii_layer = contains(StructureTree.name,LayerNames{iLayer},'IgnoreCase',true)&ii_isocortex;
%     layers(ii_layer) = iLayer;
% end

%% RUN THIS ONLY ONCE - make a volume with layers (from 1 to 5)

% % this should work? but it takes a long time
% % LayerVolume = layers(AnnotatedVolume);
 
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

%% Assign a cortical area to every structure (zero if not in cortex)

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

% I should give this a name!
StructureTree(iiCorticalAreas,:)

%% Downsample it and smooth it

DecFac = 5; % 10 can run on my laptop (32 GB), 5 requires a workstation (128 GB)

% decimate it 
DecLayerVolume = double(LayerVolume(1:DecFac:nn(1),1:DecFac:nn(2),1:DecFac:nn(3)));
dnn = nn/DecFac;


%% find the gradients

% blur it
H = fspecial3('ellipsoid',(30/DecFac)*[1 1 1]); 
SmoothDecLayerVolume = imfilter(DecLayerVolume,H,'replicate');

ShowSections( SmoothDecLayerVolume, iX, iZ, DecFac );

[ DnGradY, DnGradX, DnGradZ] = gradient(SmoothDecLayerVolume,DecFac); 

% I no longer need the smooth volume
clear H SmoothDecLayerVolume

% % set the gradients to zero outside cortex
DnGradX = DnGradX.*(DecLayerVolume<6).*(DecLayerVolume>0);
DnGradY = DnGradY.*(DecLayerVolume<6).*(DecLayerVolume>0);
DnGradZ = DnGradZ.*(DecLayerVolume<6).*(DecLayerVolume>0);

% figure; clf; ax = [];
% ax(1) = subplot(1,3,1); imagesc( squeeze( DnGradX(iX/DecFac,:,:) ) ); axis image
% ax(2) = subplot(1,3,2); imagesc( squeeze( DnGradY(iX/DecFac,:,:) ) ); axis image
% ax(3) = subplot(1,3,3); imagesc( squeeze( DnGradZ(iX/DecFac,:,:) ) ); axis image
% colormap bone

ax = ShowSections( DecLayerVolume, iX, iZ, DecFac );
 
xx = 1:4:dnn(1);
yy = 1:4:dnn(2);
zz = 1:4:dnn(3);

axes(ax(1)); quiver( zz, yy, squeeze(DnGradZ(iX/DecFac,yy,zz)), squeeze(DnGradY(iX/DecFac,yy,zz)) )
axes(ax(2)); quiver( xx, yy, squeeze(DnGradX(xx,yy,iZ/DecFac))', squeeze(DnGradY(xx,yy,iZ/DecFac))' )
clear xx yy zz
% colormap hot

% paths to the surface
StartPoints = [[37,37];[41,33];[29,29];[37,45]]*10/DecFac; % or you could choose all the points in L6
for iStartPoint = 1:4  
    start = StartPoints(iStartPoint,:);
    StLnIn  = streamline( DnGradX(:,:,iZ/DecFac)', DnGradY(:,:,iZ/DecFac)',start(1),start(2),[2 100]);
    StLnOut = streamline(-DnGradX(:,:,iZ/DecFac)',-DnGradY(:,:,iZ/DecFac)',start(1),start(2),[2 100]);
end
clear iStartPoint start StartPoints

%% find all the edges of layer 2 with layer 1 or with the outside 

CortexSurf = false(dnn);
for ind = find(DecLayerVolume == 2)'
    [i,j,k]=ind2sub(dnn,ind);
    Comparisons = DecLayerVolume(i-1:i+1, j-1:j+1, k-1:k+1)<2;
    CortexSurf(ind)=any(Comparisons(:)); 
end

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
lighting gouraud
axis equal

view(3); 
axis tight
camlight headlight

%% Load Pip's colormap

addpath('Data');
MatteoColormap; % loads cMap
PipColorMap = [ cell2mat(cMap(:,2)); [255 255 255]]/255;
clear cMap

%% another graphic in 3D

nSurfPoints = nnz(CortexSurf)
SurfCoords = zeros(3,nSurfPoints);
[SurfCoords(1,:),SurfCoords(2,:),SurfCoords(3,:)] = ind2sub( dnn, find(CortexSurf) );

% the rows in StructureTree for each entry in X
SurfRows = areas(AnnotatedVolume(sub2ind(nn,...
    SurfCoords(1,:)*DecFac,...
    SurfCoords(2,:)*DecFac,...
    SurfCoords(3,:)*DecFac)));

SurfAreas = SurfRows;
for iCorticalArea = 1:nCorticalAreas
    SurfAreas( SurfRows==iiCorticalAreas(iCorticalArea) ) = iCorticalArea;
end
SurfAreas(SurfAreas==0) = nCorticalAreas+1;

figure; clf
scatter3(SurfCoords(1,:),SurfCoords(2,:),SurfCoords(3,:),20, SurfAreas, 'filled');
axis equal
axis off

colormap(PipColorMap);

%% Isomap to flatten it -- RUN THIS ONCE

% save X X
% 
% % Run this on ZULTRA
% load X
% addpath('Isomap');
% [Xstrain,XGraph] = Run3dIsomap(X);
% save Xstrain50 Xstrain
% save XGraph50 XGraph

%% OR LOAD THIS
load('Data/Xstrain50');
load('Data/XGraph50');

%% Graphics

figure; clf
p = plot(XGraph); hold on;
p.XData = Xstrain(1,:);
p.YData = Xstrain(2,:);
p.EdgeColor = 'k';

scatter(Xstrain(1,:),Xstrain(2,:),20,SurfAreas, 'filled'); hold on
axis('equal'); axis tight; axis off

colormap(PipColorMap);

for iCorticalArea = 1:nCorticalAreas
    Acronym = StructureTree.acronym( iiCorticalAreas(iCorticalArea) );
    pp = (SurfRows == iiCorticalAreas(iCorticalArea));
    [~,c] = kmedoids( Xstrain(:,pp)', 1 );
    text( c(1),c(2), Acronym,'hori','center','vert','middle', 'color', 'k' ); 
end

%% Now do a Delaunay triangulation on the 2-D data

TriangIn2D = delaunayTriangulation( Xstrain(1,:)',  Xstrain(2,:)' );

figure; clf;
triplot(TriangIn2D);
axis equal; axis tight; axis off

%% Look at the triangulation in 3-D

% plot the triangulation in 3D coords

TriangIn3D = triangulation(TriangIn2D.ConnectivityList,SurfCoords');

figure;
trisurf(TriangIn3D, 'FaceColor','cyan','FaceAlpha',0.8);
axis equal
camlight 
lighting gouraud
title('Delaunay triangulation in 3D');

%% Now define a cartesian grid and ensure it is inside the 2D cortex

CartScale = 4;
ii = floor(min(Xstrain(1,:))):1/CartScale:ceil(max(Xstrain(1,:)));
jj = floor(min(Xstrain(2,:))):1/CartScale:ceil(max(Xstrain(2,:)));

[jjj,iii] = meshgrid(ii,jj);
Cart2D = [jjj(:), iii(:)];

BoundaryPoints = boundary( Xstrain', 1);
ps = polyshape(Xstrain(:, BoundaryPoints)');
in = isinterior(ps,Cart2D);

Cart2D = Cart2D(in,:);

% add the boundary points, it will be useful to have them
% Cart2D = [ Cart2D; Xstrain(:, BoundaryPoints)' ];

figure; clf
plot(ps); hold on 
plot( Cart2D(:,1), Cart2D(:,2), 'r.' );
title('Cartesian grid in 2D');

%% Assign barycentric coordinates to the cartesian grid

% assign each point to a triangle
IDs = pointLocation(TriangIn2D,Cart2D);

% figure; clf;
% triplot(TriangIn2D); hold on 
% axis equal; axis tight; axis off
% plot( Cart2D(:,1), Cart2D(:,2), 'k.' )

Bary2D = cartesianToBarycentric(TriangIn2D,IDs,Cart2D);

%% now plot the grid in 3D using their barycentric coords

Cart3D = barycentricToCartesian(TriangIn3D,IDs,Bary2D);
Cart3D = round(Cart3D*DecFac); % scale back up to full resolution

figure;
plot3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),'ko');

CartRows = areas(AnnotatedVolume(sub2ind(nn,Cart3D(:,1),Cart3D(:,2),Cart3D(:,3))));
CartAreas = CartRows;
for iCorticalArea = 1:nCorticalAreas
    CartAreas( CartRows==iiCorticalAreas(iCorticalArea) ) = iCorticalArea;
end

figure; clf
scatter3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),20, CartAreas, 'filled');
colormap(PipColorMap);
axis equal
title('Cartesian 2D grid plotted in 3D');

%% Make a nice flattened map of the cortex, with labels

figure; clf
scatter(Cart2D(:,1),Cart2D(:,2),20, CartAreas, 'filled');
colormap(PipColorMap);
axis equal
title('Cartesian 2D grid plotted in 2D');

%% Make an image of it

iiCart2D = Cart2D;
iiCart2D(:,1) = CartScale*(Cart2D(:,1) - min(Cart2D(:,1)) )+1;
iiCart2D(:,2) = CartScale*(Cart2D(:,2) - min(Cart2D(:,2)) )+1;

AreaImg = 0*iii;
s = sub2ind(size(iii),iiCart2D(:,2),iiCart2D(:,1));
AreaImg(s) = CartAreas;

AreaInfo = regionprops(AreaImg);
for iArea = 1:nCorticalAreas
    ThisAcronym = StructureTree.acronym(iiCorticalAreas(iArea));
    AreaInfo(iArea).Acronym = ThisAcronym{1};
end

% clean up the noise
CleanImg = AreaImg;
for iArea =  [21, 26, 1:nCorticalAreas, 24]
    bw = logical(CleanImg == iArea);
    CleanImg(bw) = 0;  % unassigned
    bw = imclose( bw, strel('disk',5));
    CleanImg(bw) = iArea;
end


    
% need to do this for graphics (so background can be white)
AreaImg(AreaImg==0)= nCorticalAreas+1;
CleanImg(CleanImg==0)= nCorticalAreas+1;




figure; clf; 
image(CleanImg); axis equal; axis image; axis off
colormap(PipColorMap); 
set(gca,'ydir','normal');
for iA = 1:nCorticalAreas
    text( AreaInfo(iA).Centroid(1), AreaInfo(iA).Centroid(2), [num2str(iA) ' ' AreaInfo(iA).Acronym], 'hori', 'center', 'vert', 'middle' );
end

%% Get vector graphics for the boundaries

% this "almost" works, except that some areas like ACAv and appear twice,
% and ILA has an issue. And VISpor too maybe?

% this works fine but it does not produce vector graphics -- perhaps it's
% good enough?

% mask = boundarymask(CleanImg);
% figure; imagesc(mask)
% set(gca,'ydir','normal');
% axis equal
% colormap(1-bone)

figure; clf; 
image(CleanImg); axis equal; axis image; axis off
colormap(PipColorMap); 
set(gca,'ydir','normal');
for iA = 1:nCorticalAreas
    text( AreaInfo(iA).Centroid(1), AreaInfo(iA).Centroid(2), [num2str(iA) ' ' AreaInfo(iA).Acronym], 'hori', 'center', 'vert', 'middle' );
end

for iArea = 1:nCorticalAreas
    bw = logical(CleanImg == iArea);
    [B,L] = bwboundaries(bw,'noholes'); hold on;
    CleanImg(L==1) = iArea; % this modifies the image BUT IT CAN LEAVE DETRITUS OUTSIDE
    k = 1;
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 3);
    drawnow
    % input('')
end
set(gca,'ydir','normal');


% Note to self: for outlines of areas in 3D, use the info shown here:
% https://uk.mathworks.com/help/images/ref/reducepoly.html
%%


b = boundarymask( AreaImg, 4 );


%%

% ACAd and RSPv are split in two...
% medoid works better than centroid...

I = uint8(AreaImg == 1);
bw = imbinarize(I);
[B,L] = bwboundaries(bw,'noholes');

imshow(I)
hold on;
k = 1;
boundary = B{k};
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
hold off;
%%
% Use |reducepoly| to reduce the number of points defining the coin boundary. 
p = [boundary(:,2) boundary(:,1)];
tolerance = 0.02; % choose suitable tolerance
p_reduced = reducepoly(p,tolerance);
%%
% Compare the original polygon overlaid over the reduced polygon and see
% how well the shape defined by fewer vertices matches the original polygon.
hf = figure;
ha = axes('parent',hf,'box','on','Ydir','reverse');
axis equal
% Original data.
line(p(:,1),p(:,2),'parent',ha,...
      'color',[1 0.5 0],'linestyle','-','linewidth',1.5,...
      'marker','o','markersize',4.5)
% Reduced data.
line(p_reduced(:,1),p_reduced(:,2),'parent',ha,...
       'color',[0 0 1],'linestyle','-','linewidth',2,...
       'marker','o','markersize',5);
legend('Original points','Reduced points');
title('Douglas-Peucker algorithm'); 




%%

% in principle now every point in 3D can be assigned to a point in 2D?

% point in cortex -> point on surface (using gradients)
% point on surface -> barymetric coords (using 3D triang)
% barymetric coords -> point on plane (using 2D triang)


%% and reversing the direction: 

% points in 2D should appear in 3D

% also, use that thing about superpixels to find smooth contours? and then plot
% them in 3D



