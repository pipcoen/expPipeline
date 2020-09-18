function subsampledData = makeFreqUniform(vectorOfIndicies, numberOfShuffles, outputData, maxFreq)
%%
if ~exist('numberOfShuffles', 'var'); numberOfShuffles = 1; end
if ~exist('maxFreq', 'var'); maxFreq = inf; end
uniValues = unique(vectorOfIndicies, 'stable');
frqValues = arrayfun(@(x) sum(vectorOfIndicies==x), uniValues);
minThresh = min([frqValues(:); maxFreq]);
if isempty(frqValues); minThresh = []; end
frqCells = arrayfun(@(x) [ones(minThresh,1); zeros(x-minThresh,1)], frqValues, 'uni', 0);
frqCells = repmat(frqCells,1,numberOfShuffles);
subsampledData = cell2mat(cellfun(@(x) x(randperm(length(x))), frqCells, 'uni', 0))>0;
if exist('outputData', 'var') && ~isempty(outputData)
    separatedShuffles = num2cell(subsampledData,1);
    subsampledData = cell2mat(cellfun(@(x) outputData(x>0), separatedShuffles, 'uni', 0));
end
end