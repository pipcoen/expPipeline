function expList = scanForNewFiles(rebuildList)
%% A funciton to check for new files for specified mice, copy them to a local directory, and update a related experimental list object 

% Inputs(default values)
% rebuildList(0)----A tag to rebuild entire experimental list. Should be used if new files need to be excluded or new fields populated.

% Outputs
% expList-----------A structure containing fields for every file detected for the selected mice in the selected daterange
    %.subject-----------Name of the mouse
    %.expDate-----------Date that the experiment was recorded
    %.sessionNum--------Session number for experiment
    %.expDef------------The experimental definition file used
    %.sharedData--------Location to uploaded copy of processed file on zserver
    %.excluded----------Tag for excluded files (these may be excluded for a variety of reasons)
    %.rigName-----------Name of the rig where the experiment took place
    %.rigType-----------Type of the rig where the experiment took place
    %.expDuration-------Duration of the experiment (s)
    %.blockFunction-----Name of helper function for processing the experiment
    %.validRepeat-------Tag for a valid repeat (in cases where there were two experiments of same type, on the same day, on the same rig.
    %.rawBlock----------Path to the local copy of the block file (generated in all experiments)
    %.rawTimeline-------Path to the local copy of the timeline file
    %.galvoLog----------Path to the local copy of the galvo log file (generated for inactivation expeirments)
    %.rawParams---------Path to the local copy of the parameters file (generated in all experiments)
    %.processedData-----Path to the processed data file (which will is generated by "convertExpFiles.m" for all experiments)
    %.rawFolder---------Path to the folder were raw data is recorded on zserver
    %.suite2POutput-----Path to the processed suite2p file (which will is generated by "convertExpFiles.m" for 2p experiments)

%#ok<*AGROW>
%% Define which mice should be included and which start/end dates (defaults to all files for that mouse)
if ~exist('rebuildList', 'var'); rebuildList = 0; end
expInfo = prc.pathFinder('expInfo');
includedMice = {'Dylan'; 'PC005'; 'PC006'; 'PC009'; 'PC010'; 'PC011'; 'PC012'; 'PC013'; 'PC014'; 'PC015'; 'PC016'; 'PC017'; 'PC018';'PC019'; ...
    'PC020'; 'PC021'; 'PC022'; 'PC023'; 'PC024'; 'PC025'; 'PC026'; 'PC027'; 'PC028'; 'PC029'; 'DJ001'; 'DJ002'; 'DJ003'; 'DJ004'; 'DJ005'; ...
    'DJ006'; 'DJ007';};
startedDates = {'', ''};
retiredDates = {...
    'PC005' '2017-07-28'; ...
    'PC006' '2017-07-28'; ...
    'PC009' '2017-03-22'; ...
    'PC014' '2017-07-01'};
includedMice(:,2) = {datenum('2015-09-01', 'yyyy-mm-dd')};
includedMice(:,3) = {datenum(floor(now))};

[~, startedIdx] = ismember(includedMice(:,1), startedDates(:,1));
if any(startedIdx); includedMice(startedIdx > 0,3) = num2cell(datenum(startedDates(startedIdx(startedIdx>0), 2), 'yyyy-mm-dd')); end
[~, retiredIdx] = ismember(includedMice(:,1), retiredDates(:,1));
if any(retiredIdx); includedMice(retiredIdx > 0,3) = num2cell(datenum(retiredDates(retiredIdx(retiredIdx>0), 2), 'yyyy-mm-dd')); end

%% List of rig names and the type of rig (possible values are training, twophoton, widefield
rigList = {'zym1', 'training'; 'zym2', 'training'; 'zym3', 'training'; 'zym4', 'training'; 'zredone', 'training'; 'zredtwo', 'training'; 'zoolatry', 'training'; ...
    'zredthree', 'training'; 'zgreyfour', 'training'; 'zurprise', 'twophoton'; 'zooropa', 'widefield'; 'zgood', 'widefield'};

%% Check for all files generated in the past 3 weeks for included mice.
lastBlockFiles = [];
lastWeekOnly = [];
lastWeek = datestr(datenum(datetime('now'))-7:datenum(datetime('now')), 'yyyy-mm-dd');
for i = 1:size(includedMice,1)
    lastWeeks = datestr(datenum(includedMice{i,3})-20:datenum(includedMice{i,3}), 'yyyy-mm-dd');
    detectedFiles = arrayfun(@(x) dir([expInfo '/' includedMice{i} '/' lastWeeks(x,:) '/**/' '*block*']), 1:length(lastWeeks), 'uni', 0);
    lastBlockFiles = cat(1, lastBlockFiles, detectedFiles{:});
    
    detectedFiles = arrayfun(@(x) dir([expInfo '/' includedMice{i} '/' lastWeek(x,:) '/**/' '*block*']), 1:size(lastWeek,1), 'uni', 0);
    lastWeekOnly = cat(1, lastWeekOnly, detectedFiles{:});
end

%% Build list of unregistered files from the last week (or find all files if rebuilding list)
if rebuildList == 1
    expList = struct;
    newBlocks = cat(1, cellfun(@(x) dir([expInfo '/' x '/**/*block*']), includedMice(:,1), 'uni', 0));
    newBlocks = cat(1, newBlocks{:});
else 
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
    if rebuildList ~= 2
        [~, nIdx] = setdiff({lastBlockFiles.folder}',{expList.rawFolder}');
        newBlocks = lastBlockFiles(nIdx);
    else, newBlocks = lastWeekOnly;
    end
end


%% Loop to created struture for new files and add them to the expList
addedFiles  = 0;
for i = 1:length(newBlocks)
    %Identify details about the experiment from the folder name.
    splitStr = regexp(newBlocks(i).folder,'[\\/]','split');
    sessionNum = splitStr{end};
    expDate = splitStr{end-1};
    subject = splitStr{end-2};
    mouseIdx = strcmp(subject, includedMice(:,1));
    
    %Ignore files outside of the specified date range for that mouse. Otherwise, sync the raw folder with the local folder.
    if datenum(expDate, 'yyyy-mm-dd') < includedMice{mouseIdx,2} || datenum(expDate, 'yyyy-mm-dd') > includedMice{mouseIdx,3}
        if exist([prc.pathFinder('rawbackup') subject '\' expDate], 'dir'); rmdir([prc.pathFinder('rawbackup') subject '\' expDate], 's'); end
        continue;
    else,  prc.syncfolder([prc.pathFinder('expinfo') subject '\' expDate], [prc.pathFinder('rawbackup') subject '\' expDate], 2);
    end

    fprintf('Adding recording %s %s %s\n', expDate, subject, sessionNum);
    addedFiles = 1;
    %Poplate fields for addition to expList
    expList(end+(length(fields(expList))>0*1),1).subject = subject;
    expList(end).expDate = expDate;
    expList(end).sessionNum = sessionNum;
    expList(end).expDef = 'ChoiceWorld';
    
    backUpFolder = prc.pathFinder('rawbackupfolder', subject, expDate, sessionNum);
    prc.syncfolder(newBlocks(i).folder, backUpFolder, 2);

    if exist(prc.pathFinder('galvoLog', subject, expDate, sessionNum), 'file')
        expList(end).galvoLog = prc.pathFinder('galvoLog', subject, expDate, sessionNum);
    else
        expList(end).galvoLog = 0;
    end
      
    %If a text file with "Exclude" in the name is detected in the experiment directory, exclude this experiment.
    %This exclude is placed here to avoid wasting time in loading unwanted block files, and 
    txtF = [dir([fileparts(newBlocks(i).folder) '\*Exclude*.txt']); dir([newBlocks(i).folder '\*Exclude*.txt'])];
    if ~isempty(txtF); expList(end).excluded = 1; else, expList(end).excluded = 0; end
    if expList(end).excluded; expList(end).rigName = 'nan'; expList(end).rigType = 'nan'; continue; end
    
    load([backUpFolder '/' newBlocks(i).name]); b = block;    
    if ~contains(b.rigName, rigList(:,1)); error([b.rigName ' not recognized']); end
    if isfield(b, 'expDef'); [~, expList(end).expDef] = fileparts(b.expDef); end
    expList(end).rigName = b.rigName;
    expList(end).rigType = rigList{strcmp(b.rigName, rigList(:,1)),2};
    if isfield(b, 'duration'); expList(end).expDuration = b.duration; else, expList(end).expDuration = 0; end
    expList(end).blockFunction = str2func(['prc.expDef.' expList(end).expDef]);
    expList(end).validRepeat = 0;
    expList(end).excluded = 0;
    
   
    %If a ChoiceWorld experiment, or experiment lasts less than 60s, then ignore (Pip doesn't use Choiceworld)
    if strcmp(expList(end).expDef, 'ChoiceWorld'); expList(end).excluded = 1; continue; end
    if expList(end).expDuration < 60; expList(end).excluded = 1; continue; end

    %If a text file with "NewFOV" in the name is detected in the experiment directory, alter suite2Poutput to account for multiple fields of view.
    txtF = dir([newBlocks(i).folder '\*NewFOV*.txt']);
    if ~isempty(txtF); newBlocks(i).folder = [fileparts(expList(end).suite2POutput) '\' txtF.name]; end
    
    %If a text file with "ValidRepeat" in the name is detected in the experiment directory, tag as a valid repeat.
    txtF = dir([newBlocks(i).folder '\*ValidRepeat*.txt']); 
    if ~isempty(txtF); expList(end).validRepeat = 1; end
end
% If new files were identified, add them to expList
if ~addedFiles; fprintf('No new files found\n'); end

%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
for i = 1:length(expList)
    expList(i).rawBlock = prc.pathFinder('rawBlock', expList(i).subject, expList(i).expDate, expList(i).sessionNum); 
    expList(i).rawTimeline = prc.pathFinder('rawTimeline', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).rawFolder = prc.pathFinder('origblockfolder', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).rawParams = prc.pathFinder('rawParameters', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).processedData = prc.pathFinder('processedData', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).sharedData = prc.pathFinder('sharedData', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).suite2POutput = prc.pathFinder('suite2POutput', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
end

%% This section section deals with cases of files on non-training rigs when multiple files are detected for a mouse on same day, rig, and expDef
folderList = cellfun(@fileparts, {expList.rawFolder}', 'uni', 0);
folderList = cellfun(@(x,y,z) [x, y, z], folderList, {expList.rigName}', {expList.expDef}', 'uni', 0);
[~, uniqueFileIdx] = unique(folderList);
duplicates = unique(folderList(setdiff(1:length(folderList),uniqueFileIdx)));
for i = 1:length(duplicates)
    tDat = expList(strcmp(folderList,duplicates{i}));
    if strcmp(tDat(1).rigType, 'training') || sum([tDat.excluded]>0)>=length(tDat)-1 || ...
            sum([tDat.excluded]>0)+sum([tDat.validRepeat])== length(tDat)
        continue;
    end
    [selectedDuplicates] = listdlg('ListString', cellfun(@num2str, {tDat.expDuration}', 'uni', 0), ...
        'PromptString', [{'Select durations to REMOVE for:'} ...
        {[tDat(1).subject ' ' tDat(1).expDate ' ' tDat(1).eDef]} {''}]);
    [tDat(:).validRepeat] = deal(1);
    [tDat(selectedDuplicates).excluded] = deal(1);
    [tDat(selectedDuplicates).validRepeat] = deal(0);
    expList(strcmp(folderList,duplicates{i})) = tDat; 
end


%% This section section deals with short files when multiple files are detected for a mouse on same day and rig type
folderList = cellfun(@fileparts, {expList.rawFolder}', 'uni', 0);
folderList = cellfun(@(x,y) [x, y], folderList, {expList.rigType}', 'uni', 0);
[~, uniqueFileIdx] = unique(folderList);
duplicates = unique(folderList(setdiff(1:length(folderList),uniqueFileIdx)));
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

%% Saving new version of expList in the shared and local directory
expList = nestedSortStruct(expList, {'subject', 'expDate'});
availableDirectories = prc.pathFinder('directoryCheck'); 
if availableDirectories(1); save(prc.pathFinder('dropboxlist'), 'expList'); save(strrep(prc.pathFinder('dropboxlist'),'dData','dDataLite'), 'expList'); end 
if availableDirectories(1); save(prc.pathFinder('sharedlist'), 'expList'); end 
end