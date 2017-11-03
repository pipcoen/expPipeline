function combinedParams = combineParams(params)
%% Function to combine and filter block files.
%Inputs(default)

fieldNames = fields(params);
combinedParams = struct;
for i = 1:length(fieldNames)
    tDat = {params(:).(fieldNames{i})}';
    sizeOftDat = cell2mat(cellfun(@(x) [size(x,1) size(x,2)], tDat, 'uni', 0));
    if length(uniquecell(tDat)) == 1
        combinedParams.(fieldNames{i})=tDat{1};
    elseif size(unique(sizeOftDat,'rows'),1) == 1 && sizeOftDat(1,2) == 1;
        combinedParams.(fieldNames{i}) = cell2mat(tDat);
    else
        combinedParams.(fieldNames{i}) = tDat;
    end
end
end
