function psychoFits(blk, splitValues, opt)
if ~exist('splitValues', 'var') || isempty(splitValues); splitValues = blk.audValues; end
colorChoices = plt.selectRedBlueColors(splitValues);
numTrials = prc.makeGrid(blk, blk.responseMade, @length, 1);
numRightTurns = prc.makeGrid(blk, blk.responseMade==2, @sum, 1);

if all(ismember(splitValues, blk.audValues)) 
    gridOfConditions = blk.grids.visValues; 
    gridOfSplits = blk.grids.audValues;
elseif all(ismember(splitValues, blk.visValues))
    gridOfConditions = blk.grids.audValues; 
    gridOfSplits = blk.grids.visValues;
else, error('Selected split values could not be found');
end

hold on; box off;
StimLevelsFineGrain = min(gridOfConditions(:)):(range(gridOfConditions(:))/1000):max(gridOfConditions(:));
opt.Linewidth = 1.2;
for j = 1:length(splitValues)
    opt.Color = colorChoices(j,:);
    idx = gridOfSplits==splitValues(j) & numTrials>0;
    [paramsValues, fittingFunction] = fit.psychoCurve(gridOfConditions(idx), numRightTurns(idx), numTrials(idx));
    plot(StimLevelsFineGrain, fittingFunction(paramsValues, StimLevelsFineGrain),opt);
end
xL = xlim; yL = ylim;
plot([0,0], yL, '--k', 'linewidth', 1.5);
plot(xL, [0.5,0.5], '--k', 'linewidth', 1.5);