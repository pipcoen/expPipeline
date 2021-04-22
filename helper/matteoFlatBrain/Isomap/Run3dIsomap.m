function [Xstrain, XGraph] = Run3dIsomap(X, v)
% RunIsomap
%
% Xstrain = Run3dIsomap(X) with X a matrix 3 x n
%
% Xstrain = Run3dIsomap(X, v) shows graphics during the process, assigning
% to each point in X a value v, where v is 1 x n (or 'k' if you want them
% to all have the same color.
%
% Xstrain = Run3dIsomap (without parameters) runs on 3D "Swiss Roll" data
%
% [Xstrain, XGraph] = Run3dIsomap() also returns a graph of 6-nearest
% neighbor adjacencies
%
% 2021 Matteo Carandini, largely copied from
% https://www.numerical-tours.com/matlab/shapes_7_isomap/
%
% I think it only works for 3D data (haven't tried with 2D)

if nargin<2
    DoGraphics = false;
else
    DoGraphics = true;
end

if nargin < 1
    % make a 3-D swiss roll
    n = 5000;
    X  = zeros(3,n);
    v = 3*pi/2 * (.1 + 2* rand(1,n)); % the values at each point
    X(2,:) = 20 * rand(1,n);
    X(1,:) = - cos( v ) .* v;
    X(3,:) = sin( v ) .* v;
else
    n = size(X,2);
end

if size(X,1)~=3
    error('X must be 3 x n');
end

if DoGraphics
    figure; clf; ax = [];
    
    % Parameter for display.
    psize = round(300/(log10(n))^2); % smaller dots if there are lots of points
    view1 = -15;
    view2 = 20;
    
    ax(1) = subplot(1,3,1);
    scatter3(X(1,:),X(2,:),X(3,:),psize,v, 'filled');
    colormap jet(256);
    view(view1,view2); axis('equal'); axis('off');
    drawnow
end

%% Pre-calcs 

% pairwise Euclidean distance matrix
EuclDist = repmat(sum(X.^2,1),n,1);
EuclDist = real(sqrt(EuclDist + EuclDist' - 2*(X'*X)));

% k-nearest neighbors (NN) connectivity
k = 6; % Number of nearest neighbors
[DNN,NN] = sort(EuclDist);
NN = NN(2:k+1,:);
DNN = DNN(2:k+1,:); % this one is k by n

% Weighted adjacency matrix (the metric on the graph).
B = repmat(1:n, [k 1]);
W = sparse(B(:), NN(:), DNN(:));

% make a graph of it
WeightedAdj = full(W);
WeightedAdj = (WeightedAdj+WeightedAdj')/2; % symmetrize it
XGraph = graph(WeightedAdj);

if DoGraphics  

    ax(2) = subplot(1,3,2);
    p = plot(XGraph); hold on;
    p.XData = X(1,:);
    p.YData = X(2,:);
    p.ZData = X(3,:);
    p.EdgeColor = 'k';
    
    scatter3(X(1,:),X(2,:),X(3,:),psize,v, 'filled');
    view(view1,view2); axis('equal'); axis('off');
    colormap jet(256);
    drawnow
end

%% Compute Pairwise Geodesic Distances

GeodesicDists = distances(XGraph,'Method','positive');

%% Remove any disconnected components

% Find index of vertices that are not connected to the main manifold.
Iremove = (GeodesicDists(:,1)==Inf);

% Remove Inf remaining values (disconnected components).
GeodesicDists(GeodesicDists==Inf) = 0;

%% Isomap with Classical Multidimensional Scaling

% Isomap performs the dimensionality reduction by applying multidimensional
% scaling.

% centered kernel
J = eye(n) - ones(n)/n;
K = -1/2 * J*(GeodesicDists.^2)*J;

% diagonalization
opt.disp = 0; 
[Xstrain, val] = eigs(K, 2, 'LR', opt);
Xstrain = Xstrain .* repmat(sqrt(diag(val))', [n 1]);
Xstrain = Xstrain';

%% Reorient the result

% Reorient the points using the two leading eigenvectors of the covariance
% matrix (PCA correction).

[U,~] = eig(Xstrain*Xstrain' / n);
Xstrain = U'*Xstrain;

% Remove problematic points.
Xstrain(:,Iremove) = Inf;

%% Display the final result 

if DoGraphics
    
    ax(3) = subplot(1,3,3);
    p = plot(XGraph); hold on;
    p.XData = Xstrain(1,:);
    p.YData = Xstrain(2,:);
    p.EdgeColor = 'k';
    
    scatter(Xstrain(1,:),Xstrain(2,:),psize,v, 'filled');
    colormap jet(256);
    axis('equal'); axis tight 
    drawnow
end

