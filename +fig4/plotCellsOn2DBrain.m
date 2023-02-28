function plotCellsOn2DBrain(hemisphere)
load('C:\Users\Pip\Dropbox (Neuropixels)\MouseData\cellData4Matteo.mat')
load('C:\Users\Pip\Dropbox (Neuropixels)\MouseData\matteoBrainData.mat')
MatteoColormapCombinedAreas;
PipColorMap = [ cell2mat(cMap(:,2)); [255 255 255]]/255;

if ~exist('hemisphere', 'var'); hemisphere = 'r'; end
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

%% Load 3D coordinates of Pip's recordings
% "cellLocations": an nx3 matrix will the Allen locations for all cells.
cellLocations = cellLocations(:,[3 2 1]); % change the order to [AP,DV,ML]
nCells = size(cellLocations,1);

% "penetrationReference": an nx1 vector indicating the penetration for each cell.
nPens = max(penetrationReference);
%% Go from 3D points to 2D points 
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
    if contains(lower(hemisphere), 'l'); title('LEFT HEMISPHERE'); hIdx = 1; end
    if contains(lower(hemisphere), 'r'); title('RIGHT HEMISPHERE'); hIdx = 2; set(gca, 'XDir','reverse'); end
    scatter( Cart2D(iiCart(hemis==hIdx),1), Cart2D(iiCart(hemis==hIdx),2), 50, 'k', 'filled', 'MarkerFaceAlpha', 0.2 ); hold on
    drawnow
    
end
axis equal




