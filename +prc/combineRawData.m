function combinedBlocks = combineRawData(raw, blocks, criterion)
%% Function to combine and filter block files.
%Inputs(default)
r.visAzimuthTimeValue = prc.indexByTrial(n, e.visAzimuthTimes', [e.visAzimuthTimes' e.visAzimuthValues'], [1 0]);

if ~exist('criterion', 'var'); criterion = vertcat(blocks(:).conditions)*0+1; end
nTrials = sum(arrayfun(@(x) size(x.trialStart,1), blocks));

if length(blocks) == 1 
    if ~isfield(blocks, 'sessionIdx'); combinedBlocks.sessionIdx = ones(nTrials, 1);
    else, combinedBlocks.sessionIdx = blocks.sessionIdx;
    end
else, combinedBlocks.sessionIdx = cell2mat(arrayfun(@(x) blocks(x).trialStart*0+x, 1:length(blocks), 'uni', 0)');
end
if length(unique(arrayfun(@(x) num2str(x.uniqueConditions(:)'),blocks,'uni',0))) > 1
    unableToMerge = 1;
    warning('Can only semi-concatenate blocks of same parameter sets. Some fields will be empty')
else, unableToMerge = 0;
end

tkIdx = criterion>0;
if ~any(tkIdx); combinedBlocks = []; return; end

fieldNames = fields(blocks);
combinedBlocks.sessionIdx = combinedBlocks.sessionIdx;
combinedBlocks.nSessions = length(unique(combinedBlocks.sessionIdx));
for i = 1:length(fieldNames)
    if strcmp(fieldNames{i}, 'sessionNum'); continue; end
    if strcmp(fieldNames{i}, 'nSessions'); continue; end
    tDat = eval(['vertcat(blocks(:).' fieldNames{i} ');']);
    if size(tDat,1) == nTrials; tDat = tDat(tkIdx,:); end

    if iscell(tDat) && length(unique(tDat))==1
        eval(['combinedBlocks.' fieldNames{i} ' = unique(tDat);']);
    elseif iscell(tDat) || size(tDat,1) == sum(tkIdx)
        eval(['combinedBlocks.' fieldNames{i} ' = tDat;']);
    elseif ~unableToMerge
        eval(['combinedBlocks.' fieldNames{i} '=blocks(1).' fieldNames{i} ';']);
    else, eval(['combinedBlocks.' fieldNames{i} '=[];']);
    end
end
