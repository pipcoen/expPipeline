
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
nRows = size(StructureTree,1);

% Template Volume: gray-scale "background signal intensity" [AP,DV,ML]
TemplateVolume = readNPY('../Allen/template_volume_10um.npy'); 

% Annotation Volume: the number at each pixel labels the area [AP,DV,ML]
AnnotatedVolume = readNPY('../Allen/annotation_volume_10um_by_index.npy'); 

% Allen CCF-bregma transform (from eyeballing Paxinos->CCF)
% this also is a function: AllenCCFBregma
bregma = [540,0,570]; % [AP,DV,ML]

%% I edited STructureTree to have cortical areas only

% save CortexTree CortexTree

load CortexTree

%% Inspect the annotation volume

iX = 870; % a place of high curvature -- used to be bregma(1);
iZ = bregma(3)-200;

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( AnnotatedVolume( iX,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( AnnotatedVolume( :,:,iZ))' ); axis image
colormap(allen_ccf_colormap('2017'));
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');

%% Assign a layer to every structure (NaN if not in cortex)

id_isocortex = StructureTree.id(strcmp(StructureTree.name,'Isocortex')); % the id of isocortex is 315
ii_isocortex = contains(StructureTree.structure_id_path, ['/' num2str(id_isocortex) '/']);

layers = nan(nRows,1);
LayerNames = {'layer 1','layer 2','layer 4','layer 5','layer 6'};
for iLayer = 1:5
    ii_layer = contains(StructureTree.name,LayerNames{iLayer},'IgnoreCase',true)&ii_isocortex;
    layers(ii_layer) = iLayer;
end


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
% save('LayerVolume','LayerVolume');

%% Load LayerVolume

load LayerVolume;
nn = size(LayerVolume);

% hack: put a 6 in the non-cortex parts of the brain
LayerVolume = LayerVolume + uint8(AnnotatedVolume>1);
LayerVolume(LayerVolume==1) = uint8(7);
LayerVolume = LayerVolume -1;

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( LayerVolume( iX,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( LayerVolume( :,:,iZ))' ); axis image
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');

%% Assign a cortical area to every structure (zero if not in cortex)

LayerRows = unique(AnnotatedVolume(LayerVolume<6&LayerVolume>0));
% these are rows of the table StructureTree
% they correspond to index + 1
% StructureTree(LayerRows,:)

areas = zeros(nRows,1);
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

UniqueAreas = unique(areas(areas>0));
nUniqueAreas = length(UniqueAreas);

StructureTree(UniqueAreas,:)

%% Downsample it and smooth it

DecFac = 5; % 10 is safe, but 5 finds the whole cortex.But 5 crashes isomap

% decimate it 
DecLayerVolume = double(LayerVolume(1:DecFac:nn(1),1:DecFac:nn(2),1:DecFac:nn(3)));
dnn = nn/DecFac;

% blur it
H = fspecial3('ellipsoid',(30/DecFac)*[1 1 1]); 
SmoothDecLayerVolume = imfilter(DecLayerVolume,H,'replicate');

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( SmoothDecLayerVolume( iX/DecFac,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( SmoothDecLayerVolume( :,:,iZ/DecFac))' ); axis image
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');
axes(ax(1)); hold on; plot(iZ/DecFac*[1 1], [1 dnn(2)])
axes(ax(2)); hold on; plot(iX/DecFac*[1 1], [1 dnn(2)])

%% find the gradients

[ DnGradY, DnGradX, DnGradZ] = gradient(SmoothDecLayerVolume,DecFac); 

% % set the gradients to zero outside cortex
DnGradX = DnGradX.*(DecLayerVolume<6).*(DecLayerVolume>0);
DnGradY = DnGradY.*(DecLayerVolume<6).*(DecLayerVolume>0);
DnGradZ = DnGradZ.*(DecLayerVolume<6).*(DecLayerVolume>0);

figure; clf; ax = [];
ax(1) = subplot(1,3,1); imagesc( squeeze( DnGradX(iX/DecFac,:,:) ) ); axis image
ax(2) = subplot(1,3,2); imagesc( squeeze( DnGradY(iX/DecFac,:,:) ) ); axis image
ax(3) = subplot(1,3,3); imagesc( squeeze( DnGradZ(iX/DecFac,:,:) ) ); axis image
colormap bone

figure; clf; ax = [];
xx = 1:4:dnn(1);
yy = 1:4:dnn(2);
zz = 1:4:dnn(3);
ax(1) = subplot(1,2,1); imagesc( squeeze( DecLayerVolume(iX/DecFac,:,:) ) ); hold on
quiver( zz, yy, squeeze(DnGradZ(iX/DecFac,yy,zz)), squeeze(DnGradY(iX/DecFac,yy,zz)) )
axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( DecLayerVolume(:,:,iZ/DecFac) )' ); hold on
quiver( xx, yy, squeeze(DnGradX(xx,yy,iZ/DecFac))', squeeze(DnGradY(xx,yy,iZ/DecFac))' )
axis image
colormap hot

%% a path to surface

StartPoints = [[37,37];[41,33];[29,29];[37,45]]*10/DecFac; % or you could choose all the points in L6

for iPoint = 1:4  
    start = StartPoints(iPoint,:);
    StLnIn  = streamline( DnGradX(:,:,iZ/DecFac)', DnGradY(:,:,iZ/DecFac)',start(1),start(2),[2 100]);
    StLnOut = streamline(-DnGradX(:,:,iZ/DecFac)',-DnGradY(:,:,iZ/DecFac)',start(1),start(2),[2 100]);
end

%% find all the edges of layer 2 with layer 1 or with the outside 

CortexSurf = false(dnn);
for ind = find(DecLayerVolume == 2)'
    [i,j,k]=ind2sub(dnn,ind);
    Comparisons = DecLayerVolume(i-1:i+1, j-1:j+1, k-1:k+1)<2;
    CortexSurf(ind)=any(Comparisons(:)); 
end

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( CortexSurf( iX/DecFac,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( CortexSurf( :,:,iZ/DecFac))' ); axis image
colormap(allen_ccf_colormap('2017'));
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');

% drop one hemisphere
CortexSurf(:,:,(dnn(3)/2+1):end) = false;

%% 3D graphics using isosurface

fv = isosurface(double(CortexSurf),0.5);

figure; clf
p = patch( fv);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight left
lighting gouraud

%% another graphic in 3D

n = nnz(CortexSurf);

[xi,yi,zi] = ind2sub( dnn, find(CortexSurf) );
X = [xi,yi,zi]';

% the rows in StructureTree for each entry in X
XRows = areas(AnnotatedVolume(sub2ind(nn,xi*DecFac,yi*DecFac,zi*DecFac)));

XAreas = XRows;
for iUniqueArea = 1:nUniqueAreas
    XAreas( XRows==UniqueAreas(iUniqueArea) ) = iUniqueArea;
end

figure; clf
scatter3(X(1,:),X(2,:),X(3,:),20, XAreas, 'filled');
axis equal
axis off

AreaColorMap = [0.9,0.9,0.9;colorcube(nUniqueAreas)];
colormap(AreaColorMap);

%% Isomap to flatten it -- RUN THIS ONCE

save X X

% Run this on ZULTRA
load X
addpath('Isomap');
[Xstrain,XGraph] = Run3dIsomap(X);
save Xstrain Xstrain
save XGraph XGraph

% [Xstrain,XGraph] = Run3dIsomap(X,XRows); 
% % colormap(allen_ccf_colormap('2017'));
% colormap(AreaColorMap);

%% OR LOAD THIS
load Xstrain50
load XGraph50


%% Graphics

figure; clf
p = plot(XGraph); hold on;
p.XData = Xstrain(1,:);
p.YData = Xstrain(2,:);
p.EdgeColor = 'k';

scatter(Xstrain(1,:),Xstrain(2,:),20,XAreas, 'filled'); hold on
axis('equal'); axis tight; axis off

colormap(AreaColorMap);

for iUniqueArea = 1:nUniqueAreas
    Acronym = StructureTree.acronym( UniqueAreas(iUniqueArea) );
    pp = (XRows == UniqueAreas(iUniqueArea));
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

TriangIn3D = triangulation(TriangIn2D.ConnectivityList,X');

figure;
trisurf(TriangIn3D, 'FaceColor','cyan','FaceAlpha',0.8);
axis equal
camlight 
lighting gouraud

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
plot( Cart2D(:,1), Cart2D(:,2), 'r.' )

%% Assign barycentric coordinates to the cartesian grid

% assign each point to a triangle
IDs = pointLocation(TriangIn2D,Cart2D);

figure; clf;
triplot(TriangIn2D); hold on 
axis equal; axis tight; axis off
plot( Cart2D(:,1), Cart2D(:,2), 'ko' )

Bary2D = cartesianToBarycentric(TriangIn2D,IDs,Cart2D);

%% now plot the grid in 3D using their barycentric coords

Cart3D = barycentricToCartesian(TriangIn3D,IDs,Bary2D);
Cart3D = round(Cart3D*DecFac); % scale back up to full resolution

figure;
plot3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),'ko');

CartRows = areas(AnnotatedVolume(sub2ind(nn,Cart3D(:,1),Cart3D(:,2),Cart3D(:,3))));
CartAreas = CartRows;
for iUniqueArea = 1:nUniqueAreas
    CartAreas( CartRows==UniqueAreas(iUniqueArea) ) = iUniqueArea;
end

figure; clf
scatter3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),20, CartAreas, 'filled');
colormap(AreaColorMap);
axis equal

%% Make a nice flattened map of the cortex, with labels

figure; clf
scatter(Cart2D(:,1),Cart2D(:,2),20, CartAreas, 'filled');
colormap(AreaColorMap);
axis equal

iiCart2D = Cart2D;
iiCart2D(:,1) = CartScale*(Cart2D(:,1) - min(Cart2D(:,1)) )+1;
iiCart2D(:,2) = CartScale*(Cart2D(:,2) - min(Cart2D(:,2)) )+1;

MyImg = 0*iii;
s = sub2ind(size(iii),iiCart2D(:,2),iiCart2D(:,1));
MyImg(s) = CartAreas;

% the colors are wrong!! <--------------------------
figure; image(MyImg);
colormap(AreaColorMap); axis equal; axis image; axis off
set(gca,'ydir','normal');

b = boundarymask( MyImg, 4 );
 
figure; imagesc(1-b), colormap bone
axis image
set(gca,'ydir','normal');

% I should do some simple dilations and erosions beforehand... <--------


%% Note to self: for outlines of areas in 3D, use the info shown here:

% https://uk.mathworks.com/help/images/ref/reducepoly.html

%% 

for iArea = 1:nUniqueAreas
    disp(StructureTree.acronym(UniqueAreas(iArea)));
    I = uint16(MyImg == iArea);
    bw = imbinarize(I);
    bw = imclose(bw,strel('disk',3));
    [B,L] = bwboundaries(bw,'noholes'); hold on;
    k = 1;
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end
set(gca,'ydir','normal');

I = uint8(MyImg == 1);
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



