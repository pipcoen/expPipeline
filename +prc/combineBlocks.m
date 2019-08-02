function combinedBlocks = combineBlocks(blocks)
%% Function to combine and filter block files.
nTrials = sum(arrayfun(@(x) size(x.trialStartEnd(:,1),1), blocks));
nEPhys = [0 0];
[~, combinedBlocks.subExpPenLink] = ismember({blocks.subject}', unique({blocks.subject}'));
combinedBlocks.subExpPenLink = num2cell([combinedBlocks.subExpPenLink (1:length(blocks))']);
if length(unique(arrayfun(@(x) num2str(x.uniqueConditions(:)'),blocks,'uni',0))) > 1
    noMerge = 1;
    fprintf('WARNING: Can only semi-concatenate blocks. Some fields will be turned into cells \n')
else, noMerge = 0;
end
blnkMat = cellfun(@(x) x(:,1)*0, {blocks.trialStartEnd}', 'uni', 0);
combinedBlocks.trialExperimentIdx = cell2mat(arrayfun(@(x) blnkMat{x}+x, 1:length(blocks), 'uni', 0)');
combinedParms = vertcat(blocks(:).params);
fieldNames = fields(blocks);

if contains('eph_spikeTimes',fieldNames)
    nEPhys = sum([cellfun(@length, {blocks.eph_spikeTimes}') cellfun(@length, {blocks.eph_clusterDepths}')]);
    combinedBlocks.subExpPenLink = [combinedBlocks.subExpPenLink, cellfun(@unique, {blocks.eph_clusterPenetrationIdx}', 'uni', 0)];
end

for fieldName = fieldNames'
    currField = fieldName{1};
    if contains(currField, [fields(combinedBlocks);'experimentNum']); continue; end
    mergeCriticalList = {'visVal'; 'audVal';'uniq';'Label';'grid'};
    alwaysCellsList = {'expT'; 'rigN';'subject';'expD';'expN';'channelMap';'clusterTemp';'lfpPower'};
    if contains(currField, alwaysCellsList); combinedBlocks.(currField) = {blocks(:).(currField)}'; continue; end
    if contains(currField, mergeCriticalList) && noMerge; combinedBlocks.(currField) = {blocks(:).(currField)}'; continue; end
    if contains(currField, {'visVal'; 'audVal'}) && ~noMerge; combinedBlocks.(currField) = blocks(1).(currField); continue; end
    
    tDat = vertcat(blocks(:).(currField));
    if iscell(tDat) && all(cellfun(@ischar, tDat)) && length(unique(tDat))==1
        combinedBlocks.(currField) = unique(tDat);
    elseif iscell(tDat) || ismember(size(tDat,1), [nTrials nEPhys])
        combinedBlocks.(currField) = tDat;
    elseif ~isstruct(tDat) && size(unique(tDat, 'rows'),1) ==1
        combinedBlocks.(currField) = blocks(1).(currField);
    elseif ~noMerge
        combinedBlocks.(currField) = blocks(1).(currField);
    else, combinedBlocks.(currField) = [];
    end
end
combinedBlocks = rmfield(combinedBlocks, 'params');
%%
allfields = fields(combinedBlocks);
firstFields = {'subject';'expDate';'expNum';'rigName';'expType';'expDef';'subExpPenLink'; 'trialExperimentIdx'; 'validTrials';'trialStartEnd';'stimPeriodStart';...
'closedLoopStart';'correctResponse';'feedback';'responseTime';'timeToWheelMove';'responseMade';'trialType';'sequentialTimeOuts';...
'timeOutsBeforeResponse';'repeatsAfterResponse';'audAmplitude';'audInitialAzimuth';'audDiff';'audValues';'visContrast';'visInitialAzimuth';...
'visDiff';'visValues';'uniqueConditions';'uniqueDiff';'uniqueConditionLabels';'conditionLabelRow';'laserType';'laserPower';'galvoType';...
'laserOnsetDelay';'galvoPosition';'laserOnOff';'grids'};
firstFields = firstFields(contains(firstFields,allfields));
nextFields = sort(allfields(~contains(allfields, '_') & ~cellfun(@(x) any(strcmp(x, firstFields)), allfields)));
extraFields = allfields(~cellfun(@(x) any(strcmp(x, firstFields)), allfields) & ~cellfun(@(x) any(strcmp(x, nextFields)), allfields));
combinedBlocks = orderfields(combinedBlocks, [firstFields;nextFields;extraFields]);
combinedBlocks.params = combinedParms;
