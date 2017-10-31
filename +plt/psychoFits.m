function psychoFits(blk, audValues)
if ~exist('audValues', 'var'); audValues = blk.audValues; end
colorChoices = plt.selectRedBlueColors(audValues);
hold on; box off;
numTrials = prc.makeGrid(blk, blk.response, @length, 1);
numRightTurns = prc.makeGrid(blk, blk.response==2, @sum, 1);
StimLevelsFineGrain = min(blk.grids.visValues(:)):(range(blk.grids.visValues(:))/1000):max(blk.grids.visValues(:));
plotOpts.Linewidth = 1.2;
for j = 1:length(audValues)
    plotOpts.Color = colorChoices(j,:);
    idx = blk.grids.audValues==audValues(j) & numTrials>0;
    [paramsValues, fittingFunction] = fit.psychoCurve(blk.grids.visValues(idx)', numRightTurns(idx), numTrials(idx));
    plot(StimLevelsFineGrain, fittingFunction(paramsValues, StimLevelsFineGrain),plotOpts);
end
xL = xlim; yL = ylim;
plot([0,0], yL, '--k', 'linewidth', 1.5);
plot(xL, [0.5,0.5], '--k', 'linewidth', 1.5);