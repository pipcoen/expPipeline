
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

CorticalAreas = StructureTree(iiCorticalAreas,:);

%% Downsample it and smooth it

DecFac = 5; % 10 can run on my laptop (32 GB), 5 requires a workstation (128 GB)

% decimate it 
DecLayerVolume = double(LayerVolume(1:DecFac:nn(1),1:DecFac:nn(2),1:DecFac:nn(3)));
dnn = nn/DecFac;


%% find the gradients

% blur it
H = fspecial3('ellipsoid',(40/DecFac)*[1 1 1]); % was 30 and it was probably fine
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

%% Plot the surface in 3D with areas colored appropriately

nSurfPoints = nnz(CortexSurf); 
SurfXYZ = zeros(3,nSurfPoints);
[SurfXYZ(1,:),SurfXYZ(2,:),SurfXYZ(3,:)] = ind2sub( dnn, find(CortexSurf) );

% the rows in StructureTree for each entry in X
SurfRows = areas(AnnotatedVolume(sub2ind(nn,...
    SurfXYZ(1,:)*DecFac,...
    SurfXYZ(2,:)*DecFac,...
    SurfXYZ(3,:)*DecFac)));

SurfAreas = SurfRows;
for iA = 1:nCorticalAreas
    SurfAreas( SurfRows==iiCorticalAreas(iA) ) = iA;
end
SurfAreas(SurfAreas==0) = nCorticalAreas+1; % so the colormap works well

% Graphics
figure; clf
scatter3(SurfXYZ(1,:),SurfXYZ(2,:),SurfXYZ(3,:),20, SurfAreas, 'filled');
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

load('Data/Xstrain50'); % loads XStrain
load('Data/XGraph50'); % loads XGraph

%% calculate the 2D medoid of each area and add it as a field in CorticalAreas

Medoid2D = zeros(nCorticalAreas,2);
for iA = 1:nCorticalAreas
    pp = (SurfAreas ==  iA);
    [~,Medoid2D(iA,:)] = kmedoids( Xstrain(:,pp)', 1 );
end
CorticalAreas = addvars(CorticalAreas,Medoid2D);

%% Plot the flattened surface in 2D with areas colored appropriately

figure; clf
p = plot(XGraph); hold on;
p.XData = Xstrain(1,:);
p.YData = Xstrain(2,:);
p.EdgeColor = 'k';

scatter(Xstrain(1,:),Xstrain(2,:),20,SurfAreas, 'filled'); hold on
axis('equal'); axis tight; axis off
colormap(PipColorMap);

for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), CorticalAreas.acronym(iA),'hori','center','vert','middle', 'color', 'k' ); 
end
title('Flattened projection of the 3D points');

%% Now do a Delaunay triangulation on the 2-D data

TriangIn2D = delaunayTriangulation( Xstrain(1,:)',  Xstrain(2,:)' );

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

%% Now define a 2D cartesian grid and ensure it is inside the cortex

CartScale = 4;
ii = floor(min(Xstrain(1,:))):1/CartScale:ceil(max(Xstrain(1,:)));
jj = floor(min(Xstrain(2,:))):1/CartScale:ceil(max(Xstrain(2,:)));

[jjj,iii] = meshgrid(ii,jj);
Cart2D = [jjj(:), iii(:)];

BoundaryPoints = boundary( Xstrain', 1);
BoundaryPoly = polyshape(Xstrain(:, BoundaryPoints)');

Cart2D = Cart2D( isinterior(BoundaryPoly,Cart2D),:);

% add the boundary points, it will be useful to have them
% Cart2D = [ Cart2D; Xstrain(:, BoundaryPoints)' ];

figure; clf
plot(BoundaryPoly); hold on 
plot( Cart2D(:,1), Cart2D(:,2), 'r.' );
title('Cartesian grid in 2D');

% should use the 2D cartesian grid to calculate the area in 2D and compare
% it to the area in 3D for each cortical area.

%% Assign barycentric coordinates to the cartesian grid

% assign each point to a triangle
Cart2DTriangs = pointLocation(TriangIn2D,Cart2D);
% find the barycentri coordinates of each point
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

%% Plot a nice flattened map of the cortex, with labels

figure; clf
scatter(Cart2D(:,1),Cart2D(:,2),20, CartAreas, 'filled');
colormap(PipColorMap);
axis equal
title('Cartesian 2D grid plotted in 2D');
for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), CorticalAreas.acronym(iA),'hori','center','vert','middle', 'color', 'k' ); 
end

%% Turn it into an image

iiCart2D = Cart2D;
iiCart2D(:,1) = CartScale*(Cart2D(:,1) - min(Cart2D(:,1)) )+1;
iiCart2D(:,2) = CartScale*(Cart2D(:,2) - min(Cart2D(:,2)) )+1;

AreaImg = 0*iii;
AreaImg( sub2ind(size(iii),iiCart2D(:,2),iiCart2D(:,1)) ) = CartAreas;
clear iiCart2D

% clean up the noise
CleanImg = AreaImg;
for iA =  [21, 26, 1:nCorticalAreas, 24]
    bw = logical(CleanImg == iA);
    CleanImg(bw) = 0;  % unassigned
    bw = imclose( bw, strel('disk',5));
    CleanImg(bw) = iA;
end

% need to do this for graphics (so background can be white)
AreaImg(AreaImg==0)= nCorticalAreas+1;
CleanImg(CleanImg==0)= nCorticalAreas+1;

figure; clf; 
image( ii, jj, CleanImg); axis equal; axis image; axis off
colormap(PipColorMap); 
set(gca,'ydir','normal');
for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), [iA CorticalAreas.acronym(iA)],'hori','center','vert','middle', 'color', 'k' ); 
end
title('Smoothed 2D map plotted as an image');

%% Find the 2D boundaries of the areas

Boundary2D = cell(nCorticalAreas,1);
for iA = 1:nCorticalAreas
    bw = logical(CleanImg == iA);
    [B,~] = bwboundaries(bw,'noholes'); hold on;
    boundary = B{1}; % there can be multiple boundaries: choose the longest
    for iBoundary = 2:numel(B)
            if numel(B{iBoundary})> numel(boundary)
                boundary = B{iBoundary};
            end
    end
    Boundary2D{iA} = boundary/CartScale - 1 + min(Cart2D(:,[2 1])); % rescale to the image coordinates
end
CorticalAreas = addvars(CorticalAreas,Boundary2D);

figure; clf; 
image(ii, jj, CleanImg); axis equal; axis image; axis off
colormap(PipColorMap); 
set(gca,'ydir','normal');
for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), CorticalAreas.acronym(iA),'hori','center','vert','middle', 'color', 'k' ); hold on
    plot( CorticalAreas.Boundary2D{iA}(:,2), CorticalAreas.Boundary2D{iA}(:,1), 'b', 'LineWidth', 3);
end
set(gca,'ydir','normal');

%% Go from 3D points to 2D points in Pip's recordings

% in principle now every point in 3D can be assigned to a point in 2D?

% point in cortex -> point on surface (using gradients)
% point on surface -> barymetric coords (using 3D triang)
% barymetric coords -> point on plane (using 2D triang)

% let's try it. 

% a mat file with the location of Pip's cells.  
load('PipCells/cellData4Matteo.mat');

% "cellLocations": an nx3 matrix will the Allen locations for all cells.
cellLocations = cellLocations(:,[3 2 1]); % change the order to [AP,DV,ML]
nCells = size(cellLocations,1);

% "penetrationReference": an nx1 vector indicating the penetration for each cell.
nPens = max(penetrationReference);

% Plot the cell locations in 3D
figure; clf
h = scatter3(SurfXYZ(1,:),SurfXYZ(2,:),SurfXYZ(3,:),20, SurfAreas, 'filled'); hold on;
set(h, 'MarkerEdgeAlpha', 0.1, 'MarkerFaceAlpha', 0.1)
hold on
axis equal


axis off
colormap(PipColorMap);
for iPen = 1:nPens
    Locs3D = cellLocations( penetrationReference==iPen, : )/DecFac;
    plot3( Locs3D(:,1), Locs3D(:,2), Locs3D(:,3), 'o' ); hold on
    % as background would be nice to have a wire frame of the brain 
    % or even better one where the wires in cortex are the areas
end

% (hopefully I got the order of the coordinates right)





% this is what I should use instead of DnGradX, etc
NormGrad = cat(4,DnGradX,DnGradY,DnGradZ);
NormGrad = NormGrad ./ repmat( vecnorm(NormGrad,2,4), [1 1 1 3] );

% Plot the cell locations in 2D
% this one does not work yet... I need to write AllenCCF2CortexFlatMap

figure; clf
% plot(BoundaryPoly); hold on
title(['Penetration ' num2str(iPen)]);
axis equal
for iA = 1:nCorticalAreas
    plot( CorticalAreas.Boundary2D{iA}(:,2), CorticalAreas.Boundary2D{iA}(:,1), 'b', 'LineWidth', 1);
    hold on
end
axis equal

for iPen = 1:nPens % iPen = 86
    % check if they are in the right hemisphere, and in cortex
    Locs3D = round( cellLocations( penetrationReference==iPen, : )/DecFac );
  
    % Locs2D = AllenCCF2CortexFlatMap(Locs3D); % I need to write this one!
    % the code below should become that function
    

    
    for iUnit = 1:size(Locs3D,1) % iUnit = 68
        p = Locs3D(iUnit,:);
        % check that it is in cortex!
        UnitLayer = DecLayerVolume(p(1),p(2),p(3));
        if UnitLayer > 0 && UnitLayer < 6 && p(3) < dnn(3)/2
            disp([iPen iUnit]);
%             % plot3( p(1), p(2), p(3), 'ko' ); hold on
%             while UnitLayer ~= 1 && UnitLayer<6 && UnitLayer >0 % bring it to Layer 1
%                 % plot3( p(1), p(2), p(3), '.' ); hold on
%                 % delta = 50*[DnGradX(p(1),p(2),p(3)),DnGradY(p(1),p(2),p(3)),DnGradZ(p(1),p(2),p(3))];
%                 delta = squeeze(NormGrad(roundp(1),roundp(2),roundp(3),:))';
%                 if isnan(delta), error('I got out of the cortex! help!'); end
%                 p = p - delta;
%                 roundp = round(p);
%                 UnitLayer = DecLayerVolume(roundp(1),roundp(2),roundp(3));
%             end
            % plot3( p(1), p(2), p(3), 'ko' ); hold on
            % drawnow
            % I should check that I am still in the cortex <------------------
            
            % now p is a point in layer 1. What triangle is it on
            % pTriang = pointLocation(TriangIn3D,p);
            % DARN THIS DOES NOT WORK NEED TO FIND ANOTHER WAY
            
            % Find the closest Cart3D point
            [~,iCart] = min( vecnorm(Cart3D - repmat(DecFac*p,[size(Cart3D,1), 1]),2,2) );
            
            % now I can plot it in 2D!
            
            plot( Cart2D(iCart,1), Cart2D(iCart,2), 'r.', 'MarkerFaceC', 'r' );
            drawnow
        end
    end
    
end




% 
% 
%             % Go up to the surface
% 
%             
%             plot3( Locs3D(iUnit,1), Locs3D(iUnit,2), Locs3D(iUnit,3), 'o' ); hold on
% 
%         % this might be totally buggy! I should check it first
%         
%     
%     end
%     
%     plot3( Locs3D(:,3), Locs3D(:,1), Locs3D(:,2), 'o' ); hold on
%     % as background would be nice to have a wire frame of the cortex where
%     % the wires are the areas 
% end





%% and reversing the direction: 

% points in 2D should appear in 3D

% also, use that thing about superpixels to find smooth contours? and then plot
% them in 3D



