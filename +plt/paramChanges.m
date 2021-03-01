function paramChanges(subject)
if ~exist('subject', 'var'); error('Must supply subject'); end
%%
s = spatialAnalysis(subject, 'all',0,0);
%%
figure;
allAVCombos = unique(cell2mat(s.blks.exp.conditionParametersAV(:)), 'rows');
conditionGrid = zeros(size(allAVCombos,1)+2, s.blks.tot.experiments);

ephys = strcmp(s.blks.exp.expType, 'ephys');
inactivations = strcmp(s.blks.exp.expType, 'inactivation');

numOfInactivationSites = 1-mat2gray(cellfun(@(x) length(find(~isnan(x(:,1)))), s.blks.exp.inactivationSites)'.*inactivations');
variableInactivationDelays = arrayfun(@(x) s.blks.tri.inactivation.laserOnsetDelay(s.blks.tri.expRef==x), 1:s.blks.tot.experiments, 'uni', 0)';
variableInactivationDelays = ~(cellfun(@(x) length(unique(x(~isnan(x))))>1, variableInactivationDelays)'.*inactivations');

changePoints = conditionGrid(1,:);
for i = 1:s.blks.tot.experiments
    conditionGrid(:,i) = [~ismember(allAVCombos, s.blks.exp.conditionParametersAV{i}, 'rows'); numOfInactivationSites(i); variableInactivationDelays(i)];
    if i > 1 && ~all(conditionGrid(:,i) == conditionGrid(:,i-1))
        changePoints(i) = 1;
    end
end
imagesc(mat2gray(conditionGrid));
colormap gray
axis square
box off;
hold on;
arrayfun(@(x) plot([x x]-0.5, ylim, 'r', 'lineWidth', 2.5), find(changePoints));
yLim = ylim;
patches = arrayfun(@(x) patch([[x x]-0.5 [x x]+0.5 x-0.5], [yLim(1) yLim(2) yLim(2) yLim(1) yLim(1)], 'c'), find(inactivations));
arrayfun(@(x) set(x, 'EdgeColor', 'none', 'FaceAlpha', 0.2) ,patches);
patches = arrayfun(@(x) patch([[x x]-0.5 [x x]+0.5 x-0.5], [yLim(1) yLim(2) yLim(2) yLim(1) yLim(1)], 'y'), find(ephys));
arrayfun(@(x) set(x, 'EdgeColor', 'none', 'FaceAlpha', 0.2) ,patches);

xTickLocations = unique([1, find(changePoints), s.blks.tot.experiments find(diff(inactivations'))+1 s.blks.tot.experiments find(diff(ephys'))+1]);
set(gca, 'xTick', xTickLocations-0.5, 'xTickLabels', s.blks.exp.expDate(xTickLocations),'XTickLabelRotation', 45);
xlabel('Experiment Date')

yTickLabels = arrayfun(@(x) ['A: ' num2str(allAVCombos(x,1)) ', V: ' num2str(allAVCombos(x,2))], 1:size(allAVCombos,1), 'uni', 0);
yTickLabels = [yTickLabels, {'lasSites'}, {'variLasDelay'}];
set(gca, 'yTick', 1:size(allAVCombos,1)+2, 'yTickLabels', yTickLabels);

end