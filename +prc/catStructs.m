function combinedStruct = catStructs(structArray)
%% Function to fulter structure.
structFields = fields(structArray);
for i = 1:length(structFields)
    if ~any(cellfun(@isstruct, {structArray.(structFields{i})}'))
        combinedStruct.(structFields{i}) = vertcat(structArray(:).(structFields{i}));
        if ischar(combinedStruct.(structFields{i}))
            combinedStruct.(structFields{i}) = num2cell(combinedStruct.(structFields{i}),2); 
        end
    else
        combinedStruct.(structFields{i}) = prc.catStructs(vertcat(structArray.(structFields{i})));
    end
end
end
