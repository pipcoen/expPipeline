function filteredStructure = filtStruct(inputStructure, criterion)
%% Function to fulter structure.
criterion = criterion>0;
refLength = length(criterion);
if ~any(criterion); filteredStructure = []; return; end

fieldNames = fields(inputStructure);
for fieldName = fieldNames'
    currField = fieldName{1};
    fieldSize = size(inputStructure.(currField));
    if ~any(fieldSize == refLength); continue; end
    if find(fieldSize==refLength)~=1; error('Fitlering only takes place along rows'); end
    inputStructure.(currField) = inputStructure.(currField)(criterion,:);
end

filteredStructure = inputStructure;
