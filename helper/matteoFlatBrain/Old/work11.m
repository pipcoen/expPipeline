
%% 2021-02-20 seeing if I can do gradients of cortical layers

% To do:
% make a surface for the surface of the cortex
% every point in surface should have a way in (use streamline) 
% make more slice views (use montage? use slice?)
% identify potential difficult spots
% project edge of layer 1 to 2D (using a mesh maybe?)
% make curves for area borders and project them to 3D and 2D
% come up with color scheme for 2D cortex
% (segmentation of regions could be done with superpixel3 maybe?)

% Not to do:
% - don't go beyond isocortex (no layers)
% - use 3D volumentric image processing tools from Matlab
%% load volumes and structure tree

addpath(genpath('../allenCCF'));
addpath(genpath('../npy-matlab'));

% A table that explains the labels
StructureTree = loadStructureTree('../Allen/structure_tree_safe_2017.csv'); 

% Template Volume: gray-scale "background signal intensity" [AP,DV,ML]
TemplateVolume = readNPY('../Allen/template_volume_10um.npy'); 

% Annotation Volume: the number at each pixel labels the area [AP,DV,ML]
AnnotatedVolume = readNPY('../Allen/annotation_volume_10um_by_index.npy'); 

% Allen CCF-bregma transform (from eyeballing Paxinos->CCF)
% this also is a function: AllenCCFBregma
bregma = [540,0,570]; % [AP,DV,ML]

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

%% RUN THIS ONLY ONCE - make a volume with layers (from 1 to 5)

% StructureTree(strcmp(StructureTree.name,'Isocortex'),:) % this is at depth = 5
% 
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

%% I edited STructureTree to have cortical areas only

% save CortexTree CortexTree

load CortexTree

%% RUN THIS ONLY ONCE - make a volume with layers (from 1 to 5)

StructureTree(strcmp(StructureTree.name,'Isocortex'),:) % this is at depth = 5

id_isocortex = StructureTree.id(strcmp(StructureTree.name,'Isocortex')); % the id of isocortex is 315
ii_isocortex = contains(StructureTree.structure_id_path, ['/' num2str(id_isocortex) '/']);

% now I would like to assign a layer to every one of those (NaN if not
% relevant), and a cortical area to every one of those (NaN if not
% relevant)

LayerNames = {'layer 1','layer 2','layer 4','layer 5','layer 6'};

layers = nan(size(StructureTree,1),1);
for iLayer = 1:5
    ii_layer = contains(StructureTree.name,LayerNames{iLayer},'IgnoreCase',true)&ii_isocortex;
    layers(ii_layer) = iLayer;
end

% areas = zeros(size(StructureTree,1),1);
% for ind = find(ii_isocortex)' % ind = 8 is one
%     newind = ind;
%     while(StructureTree.depth(newind) > 8)
%         newind = find(StructureTree.id==StructureTree.parent_structure_id(newind));
%     end
%     areas(ind) = newind;
% end
% nnz(unique(areas))
% StructureTree.name( unique(areas(areas>0)) )

% this should work? but it takes a long time
% LayerVolume = layers(AnnotatedVolume);

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

%% For each point in cortex, assign a cortical area

nRows = size(StructureTree,1);
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

StructureTree(UniqueAreas,:)

%% Downsample it and smooth it

DecFac = 10; % 10 is safe, but 5 finds the whole cortex.But 5 crashes isomap

% decimate it 
DecLayerVolume = double(LayerVolume(1:DecFac:nn(1),1:DecFac:nn(2),1:DecFac:nn(3)));

% blur it
H = fspecial3('ellipsoid',(30/DecFac)*[1 1 1]); 
SmoothDecLayerVolume = imfilter(DecLayerVolume,H,'replicate');

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( SmoothDecLayerVolume( iX/DecFac,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( SmoothDecLayerVolume( :,:,iZ/DecFac))' ); axis image
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');
axes(ax(1)); hold on; plot(iZ/DecFac*[1 1], [1 nn(2)/DecFac])
axes(ax(2)); hold on; plot(iX/DecFac*[1 1], [1 nn(2)/DecFac])

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
xx = 1:4:nn(1)/DecFac;
yy = 1:4:nn(2)/DecFac;
zz = 1:4:nn(3)/DecFac;
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

CortexBoundary = false(nn/DecFac);
for ind = find(DecLayerVolume == 2)'
    [i,j,k]=ind2sub(nn/DecFac,ind);
    Comparisons = DecLayerVolume(i-1:i+1, j-1:j+1, k-1:k+1)<2;
    CortexBoundary(ind)=any(Comparisons(:)); 
end

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( CortexBoundary( iX/DecFac,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( CortexBoundary( :,:,iZ/DecFac))' ); axis image
colormap(allen_ccf_colormap('2017'));
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');


% drop one hemisphere
CortexBoundary(:,:,(nn(3)/DecFac/2+1):end) = false;

fv = isosurface(double(CortexBoundary),0.99);

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

n = nnz(CortexBoundary);

inds = find(CortexBoundary); 
[xi,yi,zi] = ind2sub( size(CortexBoundary), inds );

X = [xi,yi,zi]';

v = areas(AnnotatedVolume(sub2ind(nn,xi*DecFac,yi*DecFac,zi*DecFac)));

u = v;
for iUniqueArea = 1:nUniqueAreas
    u( v==UniqueAreas(iUniqueArea) ) = iUniqueArea;
end

figure; clf
scatter3(X(1,:),X(2,:),X(3,:),20, u, 'filled');
colormap(colorcube(length(unique(areas))));
axis equal

nUniqueAreas = length(UniqueAreas);
AreaColorMap = [0.8,0.8,0.8;colorcube(max(UniqueAreas))];
colormap(AreaColorMap);

%% Isomap to flatten it -- RUN THIS ONCE

addpath('Isomap');
[Xstrain,XGraph] = Run3dIsomap(X,v); 

save Xstrain Xstrain
save XGraph XGraph

%% OR LOAD THIS
load Xstrain
load XGraph

% colormap(allen_ccf_colormap('2017'));
colormap(AreaColorMap);

%% Graphics

figure; clf
p = plot(XGraph); hold on;
p.XData = Xstrain(1,:);
p.YData = Xstrain(2,:);
p.EdgeColor = 'k';

scatter(Xstrain(1,:),Xstrain(2,:),20,u, 'filled'); hold on
axis('equal'); axis tight; axis off

colormap(AreaColorMap);

for iUniqueArea = 1:nUniqueAreas
    Acronym = StructureTree.acronym( UniqueAreas(iUniqueArea) );
    pp = (v == UniqueAreas(iUniqueArea));
    [~,c] = kmedoids( Xstrain(:,pp)', 1 );
    text( c(1),c(2), Acronym,'hori','center','vert','middle', 'color', 'k' ); 
end
    

%% Now do a Delaunay triangulation on the 2-D data

TriangIn2D = delaunayTriangulation( Xstrain(1,:)',  Xstrain(2,:)' );

figure; clf;
triplot(TriangIn2D);
axis equal; axis tight; axis off

%% Now define a grid and ensure it is inside the 2D cortex

ii = floor(min(Xstrain(1,:))):0.5:ceil(max(Xstrain(1,:)));
jj = floor(min(Xstrain(2,:))):0.5:ceil(max(Xstrain(2,:)));

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

%% find its barycentric coordinates

% assign each point to a triangle
IDs = pointLocation(TriangIn2D,Cart2D);

figure; clf;
triplot(TriangIn2D); hold on 
axis equal; axis tight; axis off
plot( Cart2D(:,1), Cart2D(:,2), 'ko' )

Bary2D = cartesianToBarycentric(TriangIn2D,IDs,Cart2D);

%% Then look at the triangulation on the 3-D data

% plot the triangulation in 3D coords

TriangIn3D = triangulation(TriangIn2D.ConnectivityList,X');

figure;
trisurf(TriangIn3D, 'FaceColor','cyan','FaceAlpha',0.8);
axis equal
camlight 
lighting gouraud

%% now plot the grid in 3D using their barycentric coords

Cart3D = barycentricToCartesian(TriangIn3D,IDs,Bary2D);
Cart3D = round(Cart3D*DecFac); % scale back up to full resolution

figure;
plot3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),'ko');

v = AnnotatedVolume(sub2ind(nn,Cart3D(:,1),Cart3D(:,2),Cart3D(:,3)));

figure; clf
scatter3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),20, v, 'filled');
colormap(colorcube(length(unique(areas))));
axis equal

%% Make a nice flattened map of the cortex, with labels

figure; clf
scatter(Cart2D(:,1),Cart2D(:,2),20, v, 'filled');
colormap(colorcube(length(unique(v))));
axis equal

    
    

% in principle now every point in 3D can be assigned to a point in 2D?

% point in cortex -> point on surface (using gradients)
% point on surface -> barymetric coords (using 3D triang)
% barymetric coords -> point on plane (using 2D triang)


%% and reversing the direction: 

% points in 2D should appear in 3D

% also, use that thing about superpixels to find smooth contours? and then plot
% them in 3D



