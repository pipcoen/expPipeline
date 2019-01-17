function combinedBlocks = combineBlocks(blocks)
%% Function to combine and filter block files.
nTrials = sum(arrayfun(@(x) size(x.trialStartEnd(:,1),1), blocks));
blnkMat = cellfun(@(x) x*0, {blocks.feedback}', 'uni', 0);
if length(blocks) == 1 
    if ~isfield(blocks, 'sessionIdx'); combinedBlocks = struct('sessionIdx',  ones(nTrials, 1), 'subjectIdx', ones(nTrials, 1));
    else; combinedBlocks = struct('sessionIdx',  blocks.sessionIdx, 'subjectIdx', blocks.subjectIdx);
    end
else
    combinedBlocks.sessionIdx = cell2mat(arrayfun(@(x) blnkMat{x}+x, 1:length(blocks), 'uni', 0)');
    [~, subjectReference] = ismember({blocks.subject}', unique({blocks.subject}'));
    combinedBlocks.subjectIdx = cell2mat(arrayfun(@(x) blnkMat{x}+subjectReference(x), 1:length(blocks), 'uni', 0)');
end


if length(unique(arrayfun(@(x) num2str(x.uniqueConditions(:)'),blocks,'uni',0))) > 1
    unableToMerge = 1;
    warning('Can only semi-concatenate blocks of same parameter sets. Some fields will be empty')
else, unableToMerge = 0;
end
%test
fieldNames = fields(blocks);
combinedBlocks.nSessions = length(unique(combinedBlocks.sessionIdx));
for fieldName = fieldNames'
    currField = fieldName{1};
    if strcmp(currField, 'sessionNum'); continue; end
    if strcmp(currField, 'nSessions'); continue; end
    if strcmp(currField, 'subjectIdx'); continue; end
    if contains(currField, 'rig'); continue; end
    tDat = vertcat(blocks(:).(currField));

    if iscell(tDat) && all(cellfun(@ischar, tDat)) && length(unique(tDat))==1
        combinedBlocks.(currField) = unique(tDat);
    elseif iscell(tDat) || size(tDat,1) == nTrials
       combinedBlocks.(currField) = tDat;
    elseif ~unableToMerge
        if strcmp(currField, 'expDate') 
            combinedBlocks.(currField) = tDat;
        else, combinedBlocks.(currField) = blocks(1).(currField);
        end
    else, combinedBlocks.(currField) = [];
    end
end
