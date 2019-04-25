function combinedBlocks = combineBlocks(blocks)
%% Function to combine and filter block files.
nTrials = sum(arrayfun(@(x) size(x.trialStartEnd(:,1),1), blocks));
blnkMat = cellfun(@(x) x(:,1)*0, {blocks.trialStartEnd}', 'uni', 0);
nSpikes = 0;
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

combinedParms = vertcat(blocks(:).params);
fieldNames = fields(blocks);
combinedBlocks.nSessions = length(unique(combinedBlocks.sessionIdx));

if isfield(blocks, 'ephRecordingSiteIdx')
    nSpikes = sum(arrayfun(@(x) size(x.ephRecordingSiteIdx(:,1),1), blocks));
    blnkMat = cellfun(@(x) x(:,1)*0, {blocks.ephRecordingSiteIdx}', 'uni', 0);
    combinedBlocks.ephSessionIdx = uint8(cell2mat(arrayfun(@(x) blnkMat{x}+x, 1:length(blocks), 'uni', 0)'));
    [~, subjectReference] = ismember({blocks.subject}', unique({blocks.subject}'));
    combinedBlocks.ephSubjectIdx = cell2mat(arrayfun(@(x) blnkMat{x}+subjectReference(x), 1:length(blocks), 'uni', 0)');
end

for fieldName = fieldNames'
    currField = fieldName{1};
    if contains(currField, {'sessionNum';'nSessions';'subjectIdx'}); continue; end
    if contains(currField, {'expType'; 'rigName'}); tDat = {blocks(:).(currField)}';
    else, tDat = vertcat(blocks(:).(currField));
    end

    if iscell(tDat) && all(cellfun(@ischar, tDat)) && length(unique(tDat))==1
        combinedBlocks.(currField) = unique(tDat);
    elseif iscell(tDat) || size(tDat,1) == nTrials || size(tDat,1) == nSpikes
       combinedBlocks.(currField) = tDat;
    elseif ~unableToMerge
        if strcmp(currField, 'expDate') 
            combinedBlocks.(currField) = tDat;
        else, combinedBlocks.(currField) = blocks(1).(currField);
        end
    else, combinedBlocks.(currField) = [];
    end
end
combinedBlocks.params = combinedParms;
