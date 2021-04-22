%% Tutorial on isomap

% From https://www.numerical-tours.com/matlab/shapes_7_isomap/


%% A 3-D swiss roll

n = 10000; % works for n = 1000 but chokes for 10000

% Random position on the parameteric domain.
x = rand(2,n);

% Mapping on the manifold.
v = 3*pi/2 * (.1 + 2*x(1,:));
X  = zeros(3,n);
X(2,:) = 20 * x(2,:);
X(1,:) = - cos( v ) .* v;
X(3,:) = sin( v ) .* v;


%% Display the point cloud.

% Parameter for display.
ms = 5;
v1 = -15; v2 = 20;

figure; clf;
scatter3(X(1,:),X(2,:),X(3,:),ms,v, 'filled');
colormap jet(256);
view(v1,v2); axis('equal'); axis('off');

%% compute pairwise Euclidean distance matrix

D1 = repmat(sum(X.^2,1),n,1);
D1 = real(sqrt(D1 + D1' - 2*(X'*X)));

%% compute the k-nearest neighbors (NN) connectivity

% Number of NN for the graph.

k = 6;

% Compute the k-NN connectivity.

[DNN,NN] = sort(D1);
NN = NN(2:k+1,:);
DNN = DNN(2:k+1,:); % this one is k by n

%% Adjacency matrix

% Adjacency matrix, and weighted adjacency.

B = repmat(1:n, [k 1]);
A = sparse(B(:), NN(:), ones(k*n,1));

% Weighted adjacency (the metric on the graph).

W = sparse(B(:), NN(:), DNN(:));

%% plot the graph

% clf; hold on;
% scatter3(X(1,:),X(2,:),X(3,:),ms,v, 'filled');
% 
% for iNode = 1:n
%     for iEdge = 1:k
%         foo = [X(:,iNode),X(:,NN(iEdge,iNode))]';
%         plot3(foo(:,1),foo(:,2),foo(:,3));
%     end
% end
% colormap jet(256);
% view(v1,v2); axis('equal'); axis('off');
% zoom(.8);

%% Alternative using Matlab's graph functions

% from:
% uk.mathworks.com/help/matlab/math/graph-plotting-and-customization.html

D = full(W);
D = (D+D')/2;

MyGraph = graph(D);


figure; clf
p = plot(MyGraph);
p.XData = X(1,:);
p.YData = X(2,:);
p.ZData = X(3,:);
view(v1,v2); axis('equal'); axis('off');
zoom(.8);

p.EdgeColor = 'k';

hold on;
scatter3(X(1,:),X(2,:),X(3,:),ms,v, 'filled');
colormap jet(256);
view(v1,v2); axis('equal'); axis('off');
zoom(.8);

%% Floyd Algorithm to Compute Pairwise Geodesic Distances

% % A simple algorithm to compute the geodesic distances between all pairs of
% % points on a graph is Floyd iterative algorithm. Its complexity is O(n^3)
% % where n is the number of points. It is thus quite slow for sparse graph,
% % where Dijkstra runs in O(log(n)*n^2).
% 
% % Make the graph symmetric.
% D = full(W);
% D = (D+D')/2;
% 
% % Initialize the matrix.
% D(D==0) = Inf;
% 
% % Add connexion between a point and itself.
% D = D - diag(diag(D));
% 
% % Implement the Floyd algorithm to compute the full distance matrix D,
% % where D(i,j) is the geodesic distance between i and j
% for i=1:n
%     D = min(D,repmat(D(:,i),[1 n])+repmat(D(i,:),[n 1])); 
% end

%% ALTERNATIVE using Matlab

D = distances(MyGraph,'Method','positive');


%% Remove any disconnected components

% Find index of vertices that are not connected to the main manifold.
Iremove = find(D(:,1)==Inf);

% Remove Inf remaining values (disconnected components).
D(D==Inf) = 0;

%% Isomap with Classical Multidimensional Scaling

% Isomap perform the dimensionality reduction by applying multidimensional
% scaling.

% centered kernel
J = eye(n) - ones(n)/n;
K = -1/2 * J*(D.^2)*J;

% diagonalization
opt.disp = 0; 
[Xstrain, val] = eigs(K, 2, 'LR', opt);
Xstrain = Xstrain .* repmat(sqrt(diag(val))', [n 1]);
Xstrain = Xstrain';

%% plot graph

% clf; hold on;
% scatter(Xstrain(1,:),Xstrain(2,:),ms,v, 'filled'); 
% % plot_graph(A, Xstrain, options);
% for iNode = 1:n
%     for iEdge = 1:k
%         foo = [Xstrain(:,iNode),Xstrain(:,NN(iEdge,iNode))]';
%         plot(foo(:,1),foo(:,2),'k-');
%     end
% end
% colormap jet(256);
% axis('equal'); axis('off'); 

%% plot graph using Matlab graph functions

figure; clf
p = plot(MyGraph);
p.EdgeColor = 'k';
p.XData = Xstrain(1,:);
p.YData = Xstrain(2,:);

hold on;
scatter(Xstrain(1,:),Xstrain(2,:),ms,v, 'filled');
colormap jet(256);
axis('equal'); axis('off');


%% Reorient the result

% Redess the points using the two leading eigenvectors of the covariance
% matrix (PCA correction).

[U,L] = eig(Xstrain*Xstrain' / n);
Xstrain1 = U'*Xstrain;

% Remove problematic points.
Xstrain1(:,Iremove) = Inf;

%% Display the final result 

% figure; clf; hold on;
% scatter(Xstrain1(1,:),Xstrain1(2,:),ms,v, 'filled');
% for iNode = 1:n
%     for iEdge = 1:k
%         foo = [Xstrain1(:,iNode),Xstrain1(:,NN(iEdge,iNode))]';
%         plot(foo(:,1),foo(:,2),'k-');
%     end
% end
% colormap jet(256);
% axis('equal'); axis('off');

%% plot graph using Matlab graph functions

figure; clf
p = plot(MyGraph);
p.EdgeColor = 'k';
p.XData = Xstrain1(1,:);
p.YData = Xstrain1(2,:);

hold on;
scatter(Xstrain1(1,:),Xstrain1(2,:),ms,v, 'filled');
colormap jet(256);
axis('equal'); axis('off');

