function clustersInBrain(clusterLocations, colorLables, sameAxis)
if ~exist('colorLables', 'var'); colorLables = repmat([0 0 1], length(clusterLocations.areaIdx),1); end
if ~exist('clusterLocations', 'var'); error('clusterLocations must be provided'); end
if ~exist('sameAxis', 'var'); sameAxis = 0; end

if sameAxis; axH = gca; else, axH = axes; end
plotBrainGrid([], axH);
hold on;
axis vis3d equal off manual
h = rotate3d(gca);
h.Enable = 'on';
view([-30,25]);
caxis([0 300]);
apMax = 1320; dvMax = 800; mlMax = 1140;
zoomLev = 50;
xlim([zoomLev,apMax-zoomLev]);
ylim([zoomLev,mlMax-zoomLev]);
zlim([zoomLev,dvMax-zoomLev]);

jitteredLocations = round(clusterLocations.position+randn(size(clusterLocations.position))*5);
for currColor = unique(colorLables, 'rows')'
    idx = ismember(colorLables, currColor', 'rows');
    scatter3(jitteredLocations(idx,3), jitteredLocations(idx,1), jitteredLocations(idx,2), 20, 'filled', 'MarkerFaceAlpha', 1, ...
        'MarkerFaceColor', currColor','MarkerEdgeColor', 'none');
end
end