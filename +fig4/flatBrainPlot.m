function flatBrainPlot
%% Flattening the mouse cortex
load matteoBrainData.mat
MatteoColormapCombinedAreas; % loads cMap
PipColorMap = [ cell2mat(cMap(:,2)); [255 255 255]]/255;
clear cMap
%% Find the 2D boundaries of the areas (slightly annoying that borders are duplicated)

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
     
    Boundary2D{iA} = [];
    [Boundary2D{iA}(:,1), Boundary2D{iA}(:,2)] = ...
        Ref2D.intrinsicToWorld(boundary(:,2),boundary(:,1));
    
end
CorticalAreas = addvars(CorticalAreas,Boundary2D);

figure; clf; 
image(Ref2D.XWorldLimits,Ref2D.YWorldLimits, CleanImg); 
axis equal; axis tight; axis xy
colormap(PipColorMap); 

for iA = 1:nCorticalAreas
    text( CorticalAreas.Medoid2D(iA,1),CorticalAreas.Medoid2D(iA,2), CorticalAreas.acronym(iA),'hori','center','vert','middle', 'color', 'k' ); hold on
    plot( CorticalAreas.Boundary2D{iA}(:,1), CorticalAreas.Boundary2D{iA}(:,2), 'b', 'LineWidth', 3);
end
export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\4a_ColouredBrain', '-pdf', '-painters');
%% Load 3D coordinates of Pip's recordings

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
end

%% Go from 3D points to 2D points 

% one day this should be a function "AllenCCF2CortexFlatMap"

figure; clf
% plot(BoundaryPoly); hold on
axis equal


for iPen = 1:nPens % iPen = 86
    disp(['Penetration ' num2str(iPen)]);
    % check if they are in the right hemisphere, and in cortex
    Locs3D = round( cellLocations( penetrationReference==iPen, : )/DecFac );
  
    % Locs2D = AllenCCF2CortexFlatMap(Locs3D); % I need to write this one!
    % the code below should become that function
   
    nLocs = size(Locs3D,1);
    iiCart = zeros(nLocs,1); % we will assign a Cart point to each (zero if we can't)
    hemis  = zeros(nLocs,1); % will be 1 for left and 2 for right
    for iUnit = 1:nLocs % iUnit = 68
        Loc = Locs3D(iUnit,:);
        UnitLayer = DecLayerVolume(Loc(1),Loc(2),Loc(3));
        if UnitLayer <1 || UnitLayer >5, continue; end % skip it: it's not in cortex
        
        if Loc(3) < dnn(3)/2
            hemis(iUnit) = 1; % it's in the left hemisphere
        else 
            hemis(iUnit) = 2; % it's in the right hemisphere
            Loc(3) = dnn(3) - Loc(3); % bring it to the left hemisphere
        end
        
        % Find the closest Cart3D point
        [~,iCart] = min( vecnorm(Cart3D - repmat(DecFac*Loc,[size(Cart3D,1), 1]),2,2) );
        iiCart(iUnit) = iCart;
    end
        
    % that's it. Each Cart3D point corresponds to a Cart2D point. Plot it.
    rand1 = (rand(length(Cart2D(iiCart(hemis==1),1)),2))*-0.5;
    rand2 = (rand(length(Cart2D(iiCart(hemis==2),1)),2))*-0.5;
    scatter(Cart2D(iiCart(hemis==1),1)+rand1(:,1), Cart2D(iiCart(hemis==1),2)+rand1(:,2), 50, 'r', 'filled', 'MarkerFaceAlpha', 0.2 ); hold on
    scatter(Cart2D(iiCart(hemis==2),1)+rand2(:,1), Cart2D(iiCart(hemis==2),2)+rand2(:,2), 50, 'g', 'filled', 'MarkerFaceAlpha', 0.2 )
    drawnow
    
end

for iA = 1:nCorticalAreas
    plot(CorticalAreas.Boundary2D{iA}(:,1), CorticalAreas.Boundary2D{iA}(:,2), 'k', 'LineWidth', 1);
    hold on
end
axis equal

export_fig('D:\OneDrive\Papers\Coen_2021\FigureParts\4a_FlatBrainWithCells', '-pdf', '-painters');

end
