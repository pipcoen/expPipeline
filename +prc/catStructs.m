function combinedStruct = catStructs(structArray)
%% Function to fulter structure.
structFields = fields(structArray);
for i = 1:length(structFields)
    if ~any(cellfun(@isstruct, {structArray.(structFields{i})}'))
        if ischar(structArray(1).(structFields{i}))
            tDat = {structArray(:).(structFields{i})};
            combinedStruct.(structFields{i}) = tDat(:);
        else, combinedStruct.(structFields{i}) = vertcat(structArray(:).(structFields{i}));
        end
    else
        combinedStruct.(structFields{i}) = prc.catStructs(vertcat(structArray.(structFields{i})));
    end
end
end
