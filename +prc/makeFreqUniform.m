function indices2Keep = makeFreqUniform(vectorOfIndicies)
%%
uniValues = unique(vectorOfIndicies);
frqValues = arrayfun(@(x) sum(vectorOfIndicies==x), uniValues);
idxValues = arrayfun(@(x) find(ismember(vectorOfIndicies, x)), uniValues,'uni', 0);
indices2Keep = sort(cell2mat(cellfun(@(x,y) randsample(x,min(frqValues),0), idxValues, 'uni', 0)));
indices2Keep = ismember(1:length(vectorOfIndicies), indices2Keep);
end