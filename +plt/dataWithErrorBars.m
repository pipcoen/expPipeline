function dataWithErrorBars(block, errorOn, audValues)
if ~exist('audValues', 'var'); audValues = block.audValues; end
if ~exist('errorOn', 'var'); errorOn = 1; end
colorChoices = plt.selectRedBlueColors(audValues);
numTrials = prc.makeGrid(block, block.responseMade, @length, 1);
numRightTurns = prc.makeGrid(block, block.responseMade==2, @sum, 1);
[prob,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), numRightTurns, numTrials, 'uni', 0);
centerPoints = cell2mat(cellfun(@(x) permute(x, [3,1,2]), prob, 'uni', 0));
lowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
highBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));

for audVal = audValues(:)'
    idx = find(block.grids.audValues==audVal & numTrials>0);
    err = [centerPoints(idx)-lowBound(idx), highBound(idx) - centerPoints(idx)];
    if errorOn == 1
        errorbar(block.grids.visValues(idx),centerPoints(idx),err(:,1),err(:,2),'.','MarkerSize',20, 'Color', colorChoices(audValues==audVal,:));
    else
        plot(block.grids.visValues(idx),centerPoints(idx),'.','MarkerSize',20, 'Color', colorChoices(audValues==audVal,:));
    end
    hold on;
end