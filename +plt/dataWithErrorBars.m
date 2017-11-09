function dataWithErrorBars(blk, errorOn, splitValues)
if ~exist('splitValues', 'var') || isempty(splitValues); splitValues = blk.audValues; end
if ~exist('errorOn', 'var'); errorOn = 1; end

colorChoices = plt.selectRedBlueColors(splitValues);
numTrials = prc.makeGrid(blk, blk.responseMade, @length, 1);
numRightTurns = prc.makeGrid(blk, blk.responseMade==2, @sum, 1);
[prob,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), numRightTurns, numTrials, 'uni', 0);
centerPoints = cell2mat(cellfun(@(x) permute(x, [3,1,2]), prob, 'uni', 0));
lowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
highBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));

if all(ismember(splitValues, blk.audValues)) 
    gridOfConditions = blk.grids.visValues; 
    gridOfSplits = blk.grids.audValues;
elseif all(ismember(splitValues, blk.visValues))
    gridOfConditions = blk.grids.audValues; 
    gridOfSplits = blk.grids.visValues;
else, error('Selected split values could not be found');
end

for splitVal = splitValues(:)'
    idx = find(gridOfSplits==splitVal & numTrials>0);
    err = [centerPoints(idx)-lowBound(idx), highBound(idx) - centerPoints(idx)];
    if errorOn == 1
        errorbar(gridOfConditions(idx),centerPoints(idx),err(:,1),err(:,2),'.','MarkerSize',20, 'Color', colorChoices(splitValues==splitVal,:));
    else
        plot(gridOfConditions(idx),centerPoints(idx),'.','MarkerSize',20, 'Color', colorChoices(splitValues==splitVal,:));
    end
    hold on;
end