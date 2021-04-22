
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

%% SKIP THIS SECTION these are the functions that would be called 

% I haven't looked at this one at all: allen_ccf_npx(tv,av,st);

% file_save_location = 'C:\Histology\Mouse1'; % where will the probe locations be saved
% probe_name = 'test'; % name probe to avoid overwriting
% 
% f = allenAtlasBrowser(TemplateVolume, AnnotatedVolume, StructureTree, file_save_location, probe_name);
% 
% plotBrainGrid();

%% example of navigation through the structure tree

% this section can be skipped too

idx = 85; % 417 is piriform cortex
StructureTree.structure_id_path(idx)

idx1 = find(StructureTree.id == StructureTree.parent_structure_id(idx ));
idx2 = find(StructureTree.id == StructureTree.parent_structure_id(idx1));
idx3 = find(StructureTree.id == StructureTree.parent_structure_id(idx2));
idx4 = find(StructureTree.id == StructureTree.parent_structure_id(idx3));
idx5 = find(StructureTree.id == StructureTree.parent_structure_id(idx4));
idx6 = find(StructureTree.id == StructureTree.parent_structure_id(idx5));

StructureTree.name(idx )
StructureTree.name(idx1)
StructureTree.name(idx2)
StructureTree.name(idx3)
StructureTree.name(idx4) % isocortex
StructureTree.name(idx5)
StructureTree.name(idx6)

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

%% make a volume with layers (from 1 to 5)

% RUN THIS ONCE ONLY
StructureTree(strcmp(StructureTree.name,'Isocortex'),:) % this is at depth = 5

id_isocortex = StructureTree.id(strcmp(StructureTree.name,'Isocortex')); % the id of isocortex is 315
ii_isocortex = contains(StructureTree.structure_id_path, ['/' num2str(id_isocortex) '/']);

LayerNames = {'layer 1','layer 2','layer 4','layer 5','layer 6'};

LayerVolume = zeros(size(AnnotatedVolume),'uint8');
for iLayer = 1:5
    ii_layer = contains(StructureTree.name,LayerNames{iLayer},'IgnoreCase',true)&ii_isocortex;
    LayerRows = find(ii_layer); % rows of st that have this layer
    for iLayerRow = 1:length(LayerRows)
        iRow = LayerRows(iLayerRow);
        disp(StructureTree.name(iRow));
        LayerVolume(AnnotatedVolume==iRow)=iLayer;
    end
end

save('LayerVolume','LayerVolume');

%% load LayerVolume

load LayerVolume;
[nx,ny,nz] = size(LayerVolume);

% hack: put a 6 in the non-cortex parts of the brain
LayerVolume = LayerVolume + uint8(AnnotatedVolume>1);
LayerVolume(LayerVolume==1) = uint8(7);
LayerVolume = LayerVolume -1;

%% Inspect LayerVolume

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( LayerVolume( iX,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( LayerVolume( :,:,iZ))' ); axis image

%% downsample it and smooth it

% decimate it brutally
DecLayerVolume = double(LayerVolume(1:10:nx,1:10:ny,1:10:nz));

% blur it
H = fspecial3('ellipsoid',[3 3 3]); % was 4
SmoothDecLayerVolume = imfilter(DecLayerVolume,H,'replicate');

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( SmoothDecLayerVolume( iX/10,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( SmoothDecLayerVolume( :,:,iZ/10))' ); axis image
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');
axes(ax(1)); hold on; plot(iZ/10*[1 1], [1 ny/10])
axes(ax(2)); hold on; plot(iX/10*[1 1], [1 ny/10])
%% find the gradients

[ DnGradY, DnGradX, DnGradZ] = gradient(SmoothDecLayerVolume,10);


% % set the gradients to zero outside cortex
DnGradX = DnGradX.*(DecLayerVolume<6).*(DecLayerVolume>0);
DnGradY = DnGradY.*(DecLayerVolume<6).*(DecLayerVolume>0);
DnGradZ = DnGradZ.*(DecLayerVolume<6).*(DecLayerVolume>0);

figure; clf; ax = [];
ax(1) = subplot(1,3,1); imagesc( squeeze( DnGradX(iX/10,:,:) ) ); axis image
ax(2) = subplot(1,3,2); imagesc( squeeze( DnGradY(iX/10,:,:) ) ); axis image
ax(3) = subplot(1,3,3); imagesc( squeeze( DnGradZ(iX/10,:,:) ) ); axis image
colormap bone

figure; clf; ax = [];
xx = 1:4:nx/10;
yy = 1:4:ny/10;
zz = 1:4:nz/10;
ax(1) = subplot(1,2,1); imagesc( squeeze( DecLayerVolume(iX/10,:,:) ) ); hold on
quiver( zz, yy, squeeze(DnGradZ(iX/10,yy,zz)), squeeze(DnGradY(iX/10,yy,zz)) )
axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( DecLayerVolume(:,:,iZ/10) )' ); hold on
quiver( xx, yy, squeeze(DnGradX(xx,yy,iZ/10))', squeeze(DnGradY(xx,yy,iZ/10))' )
axis image
colormap hot


%% do I understand gradient??

% img = false(16,20);
% img(8,10) = true;
% img = imdilate(img,strel('disk',5));
% 
% [dy,dx] = gradient(double(img));
% 
% figure; clf; ax = [];
% ax(1) = subplot(1,3,1); imagesc(img,[-1 1]); axis image
% ax(2) = subplot(1,3,2); imagesc(dx ,[-1 1]); axis image
% ax(3) = subplot(1,3,3); imagesc(dy ,[-1 1]); axis image
% colormap bone
% axes(ax(1)); hold on; quiver( dx, dy, 'r-' )
% 
% xx = 1:size(img,1);
% yy = 1:size(img,2);
% axes(ax(1)); hold on; quiver( yy, xx, dy, dx, 'r-' )


%% a path to surface

StartPoints = {[37,37],[41,33],[29,29],[37,45]}; % or you could choose all the points in L6

for iPoint = 1:4  
    start = StartPoints{iPoint};
    StLnIn  = streamline( DnGradX(:,:,iZ/10)', DnGradY(:,:,iZ/10)',start(1),start(2),[2 100]);
    StLnOut = streamline(-DnGradX(:,:,iZ/10)',-DnGradY(:,:,iZ/10)',start(1),start(2),[2 100]);
end

%% find all the edges of layer 1

[ DY, DX, DZ] = gradient(DecLayerVolume);
D = sqrt(DX.^2+DY.^2+DZ.^2);
figure; clf;
imagesc( squeeze( D(iX/10,:,:) ) ); axis image

CortexBoundary = (DecLayerVolume == 1)&(D>0.4);

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( CortexBoundary( iX/10,:,:) ) ); axis image
ax(2) = subplot(1,2,2); imagesc( squeeze( CortexBoundary( :,:,iZ/10))' ); axis image
colormap(allen_ccf_colormap('2017'));
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');

figure; clf;
imagesc( squeeze( CortexBoundary(iX/10,:,:) ) ); axis image

fv = isosurface(double(CortexBoundary),0.99);

% 
% SmoothCortexBoundary = imfilter(double(CortexBoundary),H,'replicate');
% fv = isosurface(SmoothCortexBoundary,0.1);

% % this does not work... I was going to look for the crossing 
% % of layer 1 w layer 2 but it's not quite working
% 
% fv = isosurface(SmoothDecLayerVolume.*(DecLayerVolume<=2).*(DecLayerVolume>1),2);

figure; clf
p = patch( fv);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud

% try the convex hull
% 
% inds = find(CortexBoundary); 
% [xi,yi,zi] = ind2sub( size(CortexBoundary), inds );
% 
% k = convhull( xi, yi, zi );
% 
% figure; clf
% trisurf(k,xi,yi,zi,'FaceColor','cyan')

% promising? maybe on a single hemisphere it would be ok
% but maybe it's more a matter for dilations and erosions

DilCortexBoundary = imdilate(CortexBoundary,strel('sphere',2));
fv = isosurface(DilCortexBoundary,0.5);

% there seems to be a particular issue at iX = 870; (17.5 87 61)

%% try Delaunay triangulation

% I will want to remove the triangles that "close" the cortex

inds = find(CortexBoundary); 
[xi,yi,zi] = ind2sub( size(CortexBoundary), inds );

dt = delaunayTriangulation( xi, yi, zi );

% this takes forever...
% figure;
% tetramesh(dt,'FaceAlpha',0.3);

[F,P] = freeBoundary(dt);
trisurf(F,P(:,1),P(:,2),P(:,3), ...
       'FaceColor','cyan','FaceAlpha',0.8);
axis equal

%% so the ConvexHull was the thing to do... but on a single hemisphere!


inds = find(CortexBoundary(:,:,1:nz/20)); 
[xi,yi,zi] = ind2sub( size(CortexBoundary), inds );

% dt = delaunayTriangulation( xi, yi, zi );
% C = convexHull(dt);

k = convhull( xi, yi, zi );
% figure; clf
% trisurf(k,xi,yi,zi,'FaceColor','cyan')
%  axis equal
 
t = triangulation(k,[xi,yi,zi]);

figure; clf
trisurf(t,'FaceColor','cyan')
 axis equal
camlight 
lighting gouraud


% the points that I don't want to connect are those where the connection
% goes through brain (but not isocortex)

% take the midpoint of every edge and see if it is in brain

% or the midpoint of every face, using incenter

cts = round( incenter(t)*10 ); % times 10 to upsample
ncts = size(cts,1);

rr = zeros(ncts,1);
for ic = 1:ncts
    rr(ic) = LayerVolume( cts(ic,1), cts(ic,2), cts(ic,3) );
end

% remove the triangles where the face goes through noncortex
t = triangulation(k(rr<6,:),[xi,yi,zi]);

figure; clf
trisurf(t,'FaceColor','red','EdgeColor','none')
 axis equal
camlight 
lighting gouraud

% this point is spurious

37 43 57
find( xi==37 & yi==43 & zi==57 )

%% should go back to looking at the cortex boundary (before triangulation)