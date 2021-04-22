
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

%% find all the edges of layer 1 with layer 2

CortexBoundary = false(nn/DecFac);
for ind = find(DecLayerVolume == 1)'
    [i,j,k]=ind2sub(nn/DecFac,ind);
    Comparisons = DecLayerVolume(i-1:i+1, j-1:j+1, k-1:k+1)==2;
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


%% Run isomap to flatten it

addpath('Isomap');

n = nnz(CortexBoundary);

inds = find(CortexBoundary); 
[xi,yi,zi] = ind2sub( size(CortexBoundary), inds );

X = [xi,yi,zi]';

%%

Xstrain = Run3dIsomap(X); % with no graphics

%%

v = AnnotatedVolume(sub2ind(nn,xi*DecFac,yi*DecFac,zi*DecFac));

Xstrain = Run3dIsomap(X,v);

% colormap(allen_ccf_colormap('2017'));
colormap(colorcube(length(unique(v))));

%%

AnnotatedVolume(sub2ind(nn,48*DecFac,67*DecFac,8*DecFac))
%%

% run the code in IsomapExample

figure; clf
p = plot(MyGraph);
p.EdgeColor = 'k';
p.XData = Xstrain1(1,:);
p.YData = Xstrain1(2,:);

hold on;
scatter(Xstrain1(1,:),Xstrain1(2,:),ms,v, 'filled');
colormap jet(256);
axis('equal'); axis('off');

colormap(colorcube(length(unique(v))));

%% Now do a Delaunay triangulation on the 2-D data

TriangIn2D = delaunayTriangulation( Xstrain1(1,:)',  Xstrain1(2,:)' );

figure; triplot(TriangIn2D);

%% Now define a grid and find its barycentric coordinates

ii = linspace( floor(min(Xstrain1(1,:))),ceil(max(Xstrain1(1,:))), 50 );
jj = linspace( floor(min(Xstrain1(2,:))),ceil(max(Xstrain1(2,:))), 50 );

[jjj,iii] = meshgrid(ii,jj);
Cart2D = [jjj(:), iii(:)];

hold on 
plot( Cart2D(:,1), Cart2D(:,2), 'ko' )

% should remove the triangles that are spurious
dd = distances(MyGraph,'Method','unweighted');
nConns = size(TriangIn2D.ConnectivityList,1);

maxdist = zeros(nConns,1);
for iConn = 1:nConns
    s = TriangIn2D.ConnectivityList(iConn,:);
    maxdist(iConn) = max(max(dd(s,s)));
end


figure; hist(maxdist(:))

ShortOnes = maxdist<5;

TriangIn2D = triangulation(TriangIn2D.ConnectivityList(ShortOnes,:),TriangIn2D.Points);

figure; clf
triplot(TriangIn2D); hold on 

% where the incenters are outside the cortex...

% now I should find the barycentric coordinates of all of them (the ones
% inside the cortex, not the ones outside the cortex...)


% assign each point to a triangle
IDs = pointLocation(TriangIn2D,Cart2D);

% remove the ones that are not in the triangulation
Outsiders = isnan(IDs);
Cart2D(Outsiders,:) = [];
IDs(Outsiders) = [];

figure; clf
triplot(TriangIn2D); hold on 
plot( Cart2D(:,1), Cart2D(:,2), 'ko' )

Bary2D = cartesianToBarycentric(TriangIn2D,IDs,Cart2D);

% (but only for inside points)

% also, should find the outside points, which define a 1-D ring

%% Then look at the triangulation on the 3-D data

% plot the triangulation in 3D coords

TriangIn3D = triangulation(TriangIn2D.ConnectivityList,X');


figure;
trisurf(TriangIn3D, 'FaceColor','cyan','FaceAlpha',0.8);
axis equal
camlight 
lighting gouraud

Cart3D = barycentricToCartesian(TriangIn3D,IDs,Bary2D);


%% now plot the grid in 3D using their barycentric coords

figure;
plot3(Cart3D(:,1),Cart3D(:,2),Cart3D(:,3),'ko');

%% in principle now every point in 3D can be assigned to a point in 2D?

% point in cortex -> point on surface (using gradients)
% point on surface -> barymetric coords (using 3D triang)
% barymetric coords -> point on plane (using 2D triang)


%% and reversing the direction: 

% points in 2D should appear in 3D

% also, use that thing about superpixels to find smooth contours? and then plot
% them in 3D



%% Older code











[F,P] = freeBoundary(TriangIn2D);






%% NOW I NEED TO GET RID OF SOME TRIANGLES


% THINGS ARE MESSY FROM HERE ON

% so the ConvexHull was the thing to do... but on a single hemisphere!


inds = find(CortexBoundary(:,:,1:nn(3)/20)); 
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

cts = round( incenter(t)*DecFac ); % times DecFac to upsample
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

% 37 43 57
find( xi==37 & yi==43 & zi==57 )

%% should go back to looking at the cortex boundary (before triangulation)