function [expList] = scanForNewFiles(rebuildList)
if ~exist('rebuildList', 'var'); rebuildList = false; end
%#ok<*AGROW>
expInfo = pathFinder('expInfo');
includedMice = {'PC006'; 'PC010'; 'PC011'; 'PC012'; 'PC013'; 'PC015'; 'PC016'; 'PC017'};
if rebuildList; expList = []; else; load(pathFinder('expList')); end
lastWeek = datestr(floor(now)-6:datenum(floor(now)), 'yyyy-mm-dd');
%%
lastBlockFiles = [];
for i = 1:length(includedMice)
    detectedFiles = arrayfun(@(x) dir([expInfo '/' includedMice{i} '/' lastWeek(x,:) '/**/' '*block*']), 1:7, 'uni', 0);
    lastBlockFiles = cat(1, lastBlockFiles, detectedFiles{:});
end
rigList = {'zym1', 'training'; 'zym2', 'training'; 'zym3', 'training'; ...
    'zredone', 'training'; 'zredtwo', 'training'; 'zoolatry', 'training';
    'zredthree', 'training'; 'zgreyfour', 'training'; ...
    'zurprise', 'twophoton';
    'zooropa', 'widefield'; 'zgood', 'widefield'};
%%
if rebuildList
    newBlocks = cat(1, cellfun(@(x) dir([expInfo '/' x '/**/*block*']), includedMice, 'uni', 0));
    newBlocks = cat(1, newBlocks{:});
    cellfun(@(x) syncfolder([pathFinder('expinfo') x], [pathFinder('rawbackup') x], 2), includedMice);
else 
    [~, nIdx] = setdiff({lastBlockFiles.folder}',{expList.rawFolder}'); 
    newBlocks = lastBlockFiles(nIdx);
end
tDat = struct;
for i = 1:length(newBlocks)
    splitStr = regexp(newBlocks(i).folder,'[\\/]','split');
    sessionNum = splitStr{end};
    expDate = splitStr{end-1};
    subject = splitStr{end-2};
    
    tDat(i,1).subject = subject;
    tDat(i,1).expDate = expDate;
    tDat(i,1).sessionNum = sessionNum;
    tDat(i,1).expDef = 'ChoiceWorld';
    tDat(i,1).rawBlock = pathFinder('rawBlock', subject, expDate, sessionNum);
    tDat(i,1).rawTimeline = pathFinder('rawTimeline', subject, expDate, sessionNum);
    
    backUpFolder = pathFinder('rawbackupfolder', subject, expDate, sessionNum);
    syncfolder(newBlocks(i).folder, backUpFolder, 2);

    if exist(pathFinder('galvoLog', subject, expDate, sessionNum), 'file')
        tDat(i,1).gavloLog = pathFinder('galvoLog', subject, expDate, sessionNum);
    else
        tDat(i,1).gavloLog = 0;
    end
    
    tDat(i,1).rawParams = pathFinder('rawParameters', subject, expDate, sessionNum);
    tDat(i,1).processedData = pathFinder('processedData', subject, expDate, sessionNum);
    tDat(i,1).sharedData = pathFinder('sharedData', subject, expDate, sessionNum);
    tDat(i,1).suite2POutput = pathFinder('suite2POutput', subject, expDate, sessionNum);
        
    txtF = [dir([fileparts(newBlocks(i).folder) '\*Exclude*.txt']);...
        dir([newBlocks(i).folder '\*Exclude*.txt'])];
    if ~isempty(txtF); tDat(i,1).excluded = 1;
        tDat(i,1).rawFolder = newBlocks(i).folder;
        tDat(i,1).rigName = 'nan';
        tDat(i,1).rigType = 'nan';
        continue;
    end
    
    load([backUpFolder '/' newBlocks(i).name]); b = block;
    
    fprintf('Adding recording %s %s %s\n', expDate, subject, sessionNum);
    if ~contains(b.rigName, rigList(:,1)); error([b.rigName ' not recognized']); end
    if isfield(b, 'expDef'); [~, tDat(i,1).expDef] = fileparts(b.expDef); end
    tDat(i,1).rigName = b.rigName;
    tDat(i,1).rigType = rigList{strcmp(b.rigName, rigList(:,1)),2};
    if isfield(b, 'duration'); tDat(i,1).expDuration = b.duration; else, tDat(i,1).expDuration = 0; end
    tDat(i,1).blockHelper = [tDat(i,1).expDef '_Blk_Proc.m'];
    tDat(i,1).blockFunction = str2func(tDat(i,1).blockHelper);
    tDat(i,1).rawFolder = newBlocks(i).folder;
    tDat(i,1).validRepeat = 0;
    tDat(i,1).excluded = 0;
%%
    if strcmp(tDat(i,1).expDef, 'ChoiceWorld'); tDat(i,1).excluded = 1; continue; end
    if tDat(i,1).expDuration < 60; tDat(i,1).excluded = 1; continue; end

   
    txtF = dir([tDat(i).rawFolder '\*NewFOV*.txt']);
    if ~isempty(txtF); tDat(i,1).suite2POutput = [fileparts(expList(i,1).suite2POutput) '\' txtF.name]; end
    
    txtF = dir([tDat(i).rawFolder '\*ValidRepeat*.txt']); 
    if ~isempty(txtF); tDat(i,1).validRepeat = 1; end
end
%%
if ~isempty(newBlocks); expList = [expList;tDat]; else; fprintf('No new files found\n'); end
%%
folderList = cellfun(@fileparts, {expList.rawFolder}', 'uni', 0);
folderList = cellfun(@(x,y,z) [x, y, z], folderList, {expList.rigName}', {expList.expDef}', 'uni', 0);
[~, uIdx] = unique(folderList);
duplicates = unique(folderList(setdiff(1:length(folderList),uIdx)));
for i = 1:length(duplicates)
    tDat = expList(strcmp(folderList,duplicates{i}));
    if strcmp(tDat(1).rigType, 'training') || sum([tDat.excluded]>0)>=length(tDat)-1 || ...
            sum([tDat.excluded]>0)+sum([tDat.validRepeat])== length(tDat)
        continue;
    end
    
    [selF] = listdlg('ListString', cellfun(@num2str, {tDat.expDuration}', 'uni', 0), ...
        'PromptString', [{'Select durations to REMOVE for:'} ...
        {[tDat(1).subject ' ' tDat(1).expDate ' ' tDat(1).eDef]} {''}]);
    [tDat(:).validRepeat] = deal(1);
    [tDat(selF).excluded] = deal(1);
    [tDat(selF).validRepeat] = deal(0);
    expList(strcmp(folderList,duplicates{i})) = tDat; 
end
%%
folderList = cellfun(@fileparts, {expList.rawFolder}', 'uni', 0);
folderList = cellfun(@(x,y) [x, y], folderList, {expList.rigType}', 'uni', 0);
[~, uIdx] = unique(folderList);
duplicates = unique(folderList(setdiff(1:length(folderList),uIdx)));
for i = 1:length(duplicates)
    tDat = expList((strcmp(folderList,duplicates{i}).*([expList.excluded]'==0))>0);
    if length(tDat) < 2; continue; end
    if strcmp(tDat(1).rigType, 'training')
        fprintf('Multiple training files for %s %s. Keeping largest\n', tDat(1).subject, tDat(1).expDate);
        maxDurationIdx = num2cell([tDat.expDuration] ~= max([tDat.expDuration]));
        [tDat.excluded] = maxDurationIdx{:};
        expList((strcmp(folderList,duplicates{i}).*([expList.excluded]'==0))>0) = tDat;
    elseif strcmp(tDat(1).rigType, 'twophoton') && length(unique({tDat.out2})) > 1
        newFOV = cellfun(@isempty, (strfind({tDat.suite2POutput}', 'NewFOV')));
        [containingFolder, fileExtension] = cellfun(@fileparts, {tDat.suite2POutput}', 'uni', 0);
        sExt = fileExtension(newFOV);
        sExt(1:end-1) = cellfun(@(x) [x '_'], fileExtension(1:end-1), 'uni', 0);
        [tDat(newFOV).out2] = deal([containingFolder{1} sep cell2mat(sExt')]);
        [tDat(newFOV).sessionNum] = deal({tDat(newFOV).sessionNum}');
        expList((strcmp(folderList,duplicates{i}).*([expList.excluded]'==0))>0) = tDat;
    end
end
expList = nestedSortStruct(expList, {'subject', 'expDate'});
availableDirectories = pathFinder('directoryCheck'); 
if availableDirectories(1); save(pathFinder('dropboxlist'), 'expList'); end 
if availableDirectories(1); save(pathFinder('sharedlist'), 'expList'); end 
end