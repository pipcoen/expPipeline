function combinedBlocks = combineBlocks(blocks, criterion, rawData)
%% Function to combine and filter block files.
if exist('rawData', 'var'); blocks = prc.catstruct(blocks, rawData); end
if ~exist('criterion', 'var') || isempty(criterion); criterion = vertcat(blocks(:).conditionLabel)*0+1; end
nTrials = sum(arrayfun(@(x) size(x.trialStartEnd(:,1),1), blocks));

if length(blocks) == 1 
    if ~isfield(blocks, 'sessionIdx'); combinedBlocks.sessionIdx = ones(nTrials, 1);
    else, combinedBlocks.sessionIdx = blocks.sessionIdx;
    end
else, combinedBlocks.sessionIdx = cell2mat(arrayfun(@(x) blocks(x).trialStartEnd(:,1)*0+x, 1:length(blocks), 'uni', 0)');
end
if length(unique(arrayfun(@(x) num2str(x.uniqueConditions(:)'),blocks,'uni',0))) > 1
    unableToMerge = 1;
    warning('Can only semi-concatenate blocks of same parameter sets. Some fields will be empty')
else, unableToMerge = 0;
end
%test
tkIdx = criterion>0;
if ~any(tkIdx); combinedBlocks = []; return; end

fieldNames = fields(blocks);
combinedBlocks.nSessions = length(unique(combinedBlocks.sessionIdx));
for i = fieldNames'
    field = i{1};
    if strcmp(field, 'sessionNum'); continue; end
    if strcmp(field, 'nSessions'); continue; end
%     if strcmp(field, 'rigName'); continue; end
    tDat = vertcat(blocks(:).(field));
    if size(tDat,1) == nTrials; tDat = tDat(tkIdx,:); end

    if iscell(tDat) && all(cellfun(@ischar, tDat)) && length(unique(tDat))==1
        combinedBlocks.(field) = unique(tDat);
    elseif iscell(tDat) || size(tDat,1) == sum(tkIdx)
       combinedBlocks.(field) = tDat;
    elseif ~unableToMerge
        if strcmp(field, 'expDate') 
            combinedBlocks.(field) = tDat;
        else, combinedBlocks.(field) = blocks(1).(field);
        end
    else, combinedBlocks.(field) = [];
    end
end
