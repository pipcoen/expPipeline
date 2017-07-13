function concatenatedBlocks = concatinateBlocks(blocks, criterion)
%%
if ~exist('criterion', 'var'); criterion = vertcat(blocks(:).conditions)*0+1; end
if isempty(criterion); criterion = 'none'; end

if length(unique(arrayfun(@(x) num2str(x.uniqueConditions(:)'),blocks,'uni',0))) > 1
    warning('Can only semi-concatenate blocks of same parameter sets')
end

fieldNames = fields(blocks);
nTrials = sum(arrayfun(@(x) size(x.trialStart,1), blocks));
sessionIdx = cell2mat(arrayfun(@(x) blocks(x).trialStart*0+x, 1:length(blocks), 'uni', 0)');

tkIdx = criterion>0;
sessionIdx = sessionIdx(tkIdx);
concatenatedBlocks.sessionIdx = sessionIdx;
concatenatedBlocks.nSessions = length(unique(sessionIdx));
for i = 1:length(fieldNames)
    tDat = eval(['vertcat(blocks(:).' fieldNames{i} ');']);
    if size(tDat,1) == nTrials; tDat = tDat(tkIdx,:); end
        
    if ischar(tDat); tDat = num2cell(tDat,2); end
    if iscell(tDat) && length(unique(tDat))==1
        eval(['concatenatedBlocks.' fieldNames{i} ' = unique(tDat);']);
    elseif iscell(tDat) || size(tDat,1) == sum(tkIdx)
        eval(['concatenatedBlocks.' fieldNames{i} ' = tDat;']);
    elseif contains(fieldNames{i}, {'uniqueConditions'; 'Grid'})
        eval(['concatenatedBlocks.' fieldNames{i} '=blocks(1).' fieldNames{i} ';']);
    else, error('Unregonized block field');
    end
end