%% Load Pip's data

% a mat file with the location of Pip's cells.  
load cellData4Matteo.mat

% "cellLocations": an nx3 matrix will the Allen locations for all cells.
nCells = size(cellLocations,1);

% "penetrationReference": an nx1 vector indicating the penetration for each cell.
nPens = max(penetrationReference);

%% Plot the cell locations in 3D

figure;
for iPen = 1:nPens
    Locs3D = cellLocations( penetrationReference==iPen, : );
    plot3( Locs3D(:,3), Locs3D(:,1), Locs3D(:,2), 'o' ); hold on
    % as background would be nice to have a wire frame of the brain 
    % or even better one where the wires in cortex are the areas
end

% Pip: when I plot these cells in Matlab with jitter3d, I plot columns
% (3,1,2) for (x,y,z).

%% Plot the cell locations in 2D

% this one does not work yet... I need to write AllenCCF2CortexFlatMap
figure;
for iPen = 1:nPens
    Locs3D = cellLocations( penetrationReference==iPen, : );
    Locs2D = AllenCCF2CortexFlatMap(Locs3D); % I need to write this one!
    plot3( Locs3D(:,3), Locs3D(:,1), Locs3D(:,2), 'o' ); hold on
    % as background would be nice to have a wire frame of the cortex where
    % the wires are the areas 
end

