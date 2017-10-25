function combinedParams = combineParams(params)
%% Function to combine and filter block files.
%Inputs(default)

fieldNames = fields(params);
combinedParams = struct;
for i = 1:length(fieldNames)
    tDat = eval(['{params(:).' fieldNames{i} '}'';']);
    if length(uniquecell(tDat)) == 1
        eval(['combinedParams.' fieldNames{i} ' = tDat{1};']);
    elseif length(unique(cellfun(@(x) size(x,2), tDat))) == 1
        eval(['combinedParams.' fieldNames{i} ' = cell2mat(tDat);']);
    else
        eval(['combinedParams.' fieldNames{i} ' = tDat;']);
    end
end
