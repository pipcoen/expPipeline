function expList = scanForNewFiles(rebuildList, checkDirectories, checkExpRigs)
%% A funciton to check for new files for specified mice, copy them to a local directory, and update a related experimental list object

% Inputs(default values)
% rebuildList(0)----A tag to rebuild entire experimental list. Should be used if new files need to be excluded or new fields populated.

% Outputs
% expList-----------A structure containing fields for every file detected for the selected mice in the selected daterange
%.subject-----------Name of the mouse
%.expDate-----------Date that the experiment was recorded
%.expNum--------Session number for experiment
%.expDef------------The experimental definition file used
%.sharedData--------Location to uploaded copy of processed file on zserver
%.excluded----------Tag for excluded files (these may be excluded for a variety of reasons)
%.rigName-----------Name of the rig where the experiment took place
%.expType-----------Type of the rig where the experiment took place
%.expDuration-------Duration of the experiment (s)
%.blockFunction-----Name of helper function for processing the experiment
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
if ~exist('checkDirectories', 'var'); checkDirectories = 0; end
if ~exist('checkExpRigs', 'var'); checkExpRigs = 0; end
expInfo = prc.pathFinder('expInfo');
includedMice = {'PC010'; 'PC011'; 'PC012'; 'PC013'; 'PC015'; 'PC022'; 'PC025'; 'PC027'; 'PC029'; 'PC030'; 'PC031'; 'PC032'; 'PC033'; 'PC034';...
    'PC035'; 'PC036'; 'PC037'; 'PC038'; 'PC039'; 'PC040'; 'PC041'; 'PC042'; 'PC043'; 'PC044'; ...
    'DJ006'; 'DJ007'; 'DJ008'; 'DJ010'; 'CR010'};

startedDates = {...
    'CR010' '2018-01-29'};
retiredDates = {...
    'PC005' '2017-07-28'; ...
    'PC006' '2017-07-28'; ...
    'PC009' '2017-03-22'; ...
    'PC014' '2017-07-01'; ...
    'DJ011' '2018-06-08'; ...
    'CR010' '2019-01-29'};
includedMice(:,2) = {datenum('2015-09-01', 'yyyy-mm-dd')};
includedMice(:,3) = {datenum(floor(now))};

[~, startedIdx] = ismember(includedMice(:,1), startedDates(:,1));
if any(startedIdx); includedMice(startedIdx > 0,2) = num2cell(datenum(startedDates(startedIdx(startedIdx>0), 2), 'yyyy-mm-dd')); end
[~, retiredIdx] = ismember(includedMice(:,1), retiredDates(:,1));
if any(retiredIdx); includedMice(retiredIdx > 0,3) = num2cell(datenum(retiredDates(retiredIdx(retiredIdx>0), 2), 'yyyy-mm-dd')); end

%% Check for all files generated in the past 3 weeks for included mice.
recentDates = cellfun(@(x,y) num2cell([repmat([expInfo y '\'],[15,1]) datestr(datenum(x)-14:datenum(x), 'yyyy-mm-dd')],2),includedMice(:,3), includedMice(:,1), 'uni', 0);

if rebuildList ~= 1
    cycles = 2;
    processList = vertcat(recentDates{:});
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
else
    cycles = 3;
    processList = cellfun(@(x) [expInfo x], includedMice(:,1), 'uni', 0);
    expList = struct;
end

for i = 1:cycles
    fprintf('Detecting folder level %d ... \n', i);
    processList = cellfun(@(x) java.io.File(x), processList, 'uni', 0);
    processList = cellfun(@(x) arrayfun(@char,x.listFiles,'uni',0), processList, 'uni', 0);
    processList = vertcat(processList{:});
    processList = processList(~contains(processList, {'Lightsheet';'ephys';'Backup'},'IgnoreCase',true));
end
processList = processList(~cellfun(@isempty, regexp(processList, '20.*arameters.mat')));
if ~rebuildList && isempty(processList); fprintf('No new files found\n'); return;
elseif ~rebuildList, processList = processList(~contains(processList,{expList.rawFolder}'));
end

if length(unique(cellfun(@fileparts, processList, 'uni', 0))) ~= length(processList); error('Every detected folder should be unique'); end
%% Loop to created struture for new files and add them to the expList
addedFiles  = 0;
for i = 1:length(processList)
    %Identify details about the experiment from the folder name.
    splitStr = regexp(processList{i},'[\\/]','split');
    expNum = splitStr{end-1};
    expDate = splitStr{end-2};
    subject = splitStr{end-3};
    mouseIdx = strcmp(subject, includedMice(:,1));
    
    %Ignore files outside of the specified date range for that mouse. Otherwise, sync the raw folder with the local folder.
    if datenum(expDate, 'yyyy-mm-dd') < includedMice{mouseIdx,2} || datenum(expDate, 'yyyy-mm-dd') > includedMice{mouseIdx,3}
        if exist([prc.pathFinder('rawbackup') subject '\' expDate], 'dir'); rmdir([prc.pathFinder('rawbackup') subject '\' expDate], 's'); end
        continue;
    else, prc.syncfolder([prc.pathFinder('expinfo') subject '\' expDate '\' expNum], [prc.pathFinder('rawbackup') subject '\' expDate '\' expNum], 2);
    end
    
    timeSinceParamFileCreation = dir(prc.pathFinder('origparams', subject, expDate, expNum));
    timeSinceParamFileCreation = -([timeSinceParamFileCreation.datenum]-now)*24*60;
    if ~exist(prc.pathFinder('origblock', subject, expDate, expNum), 'file') && timeSinceParamFileCreation < 80; continue; end
    
    fprintf('Adding recording %s %s %s\n', expDate, subject, expNum);
    addedFiles = 1;
    %Poplate fields for addition to expList
    expList(end+(length(fields(expList))>0*1),1).subject = subject;
    expList(end).expDate = expDate;
    expList(end).expNum = expNum;
    expList(end).expDef = 'ChoiceWorld';
    
    backUpFolder = prc.pathFinder('rawbackupfolder', subject, expDate, expNum);
    prc.syncfolder(fileparts(processList{i}), backUpFolder, 2); %#ok<*NODEF>
    
    if exist(prc.pathFinder('galvoLog', subject, expDate, expNum), 'file')
        expList(end).galvoLog = prc.pathFinder('galvoLog', subject, expDate, expNum);
    else
        expList(end).galvoLog = 0;
    end
    
    tempLoc = prc.updatePaths(expList(end), 0);
    if exist(tempLoc.rawBlock, 'file'); load(tempLoc.rawBlock, 'block'); b = block;
    else, clear b; load(processList{i});
        if exist(tempLoc.rawTimeline, 'file'); load(tempLoc.rawTimeline);
        else, expList(end).excluded = 1; continue; end
        if isfield(parameters, 'experimentType') && strcmp(parameters.experimentType, 'mpep')
            b.expDef = parameters.Protocol.xfile;
            b.rigName = 'lilrig-stim';
            b.duration = Timeline.lastTimestamp;
        else, expList(end).excluded = 1; continue;
        end
    end
    
    if isfield(b, 'expDef'); [~, expList(end).expDef] = fileparts(b.expDef); end
    expList(end).rigName = b.rigName;
    expList(end).expType = 'training';
    
    if contains(expList(end).rigName, {'zym1'; 'zym2'}) && exist(tempLoc.galvoLog, 'file'); expList(end).expType = 'inactivation'; end
    
    %     if strcmp(expList(end).rigName, 'lilrig-stim') && ~isempty(dir([fileparts(tempLoc.rawFolder) '\*hys*'])); expList(end).expType = 'widefield'; end
    if strcmp(expList(end).rigName, 'zatteo') && ~isempty(dir([tempLoc.rawFolder '\*fus.mat*'])); expList(end).expType = 'fusi'; end
    if isfield(b, 'duration'); expList(end).expDuration = b.duration; else, expList(end).expDuration = 0; end
    expList(end).blockFunction = str2func(['prc.expDef.' expList(end).expDef]);
    expList(end).excluded = 0;
    
    %If a ChoiceWorld experiment, or experiment lasts less than 60s, then ignore (Pip doesn't use Choiceworld)
    if strcmp(expList(end).expDef, 'ChoiceWorld'); expList(end).excluded = 1; continue; end
    if expList(end).expDuration < 120; expList(end).excluded = 1; continue; end
    
    txtF = [dir([fileparts(fileparts(processList{i})) '\*Exclude*.txt']); dir([fileparts(processList{i}) '\*Exclude*.txt'])];
    if ~isempty(txtF); expList(end).excluded = 1; end
    
    if contains(expList(end).expType, {'fusi'; 'ephys'}) && expList(end).excluded == 0; expList(end).excluded = -1; end
end
% If new files were identified, add them to expList
if ~addedFiles; fprintf('No new files found\n'); end
[expList(cellfun(@isempty, {expList.rigName}')).rigName] = deal('');

if checkExpRigs
    expRigs = find(contains({expList.rigName}', {'zatteo'; 'lilrig-stim'}));
    for i = expRigs
        if strcmp(expList(i).rigName, 'lilrig-stim') && ~isempty(dir([fileparts(tempLoc.rawFolder) '\*hys*'])); expList(end).expType = 'ephys'; end
    
    end
end

%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
expList = prc.updatePaths(expList,0);
if checkDirectories
    for i = 1:length(expList)
        if expList(i).excluded == 1; continue; end
        txtF = [dir([fileparts(expList(i).rawFolder) '\*Exclude*.txt']); dir([expList(i).rawFolder '\*Exclude*.txt'])];
        if ~isempty(txtF); expList(i).excluded = 1; end
    end
end
%% This section section deals with cases of files on non-training rigs when multiple files are detected for a mouse on same day, rig, and expDef
processList = cellfun(@fileparts, {expList.rawFolder}', 'uni', 0);
processList = cellfun(@(x,y,z) [x, y, z], processList, {expList.expType}', {expList.expDef}', 'uni', 0);
[~, uniqueFileIdx] = unique(processList);
duplicates = unique(processList(setdiff(1:length(processList),uniqueFileIdx)));
for i = 1:length(duplicates)
    duplicatesIdx = strcmp(processList,duplicates{i}) & [expList.excluded]'==0;
    tDat = expList(duplicatesIdx);
    if length(tDat) < 2 || (~isempty(tDat(1).expType) && contains(tDat(1).expType, {'training';'inactivation'})); continue; end
    [selectedDuplicates] = listdlg('ListString', cellfun(@(x,y) [num2str(x),': ' num2str(y)], {tDat.expNum}', {tDat.expDuration}', 'uni', 0), ...
        'PromptString', [{'Select durations to REMOVE for:'} ...
        {[tDat(1).subject ' ' tDat(1).expDate ' ' tDat(1).expDef  ' ' tDat(1).expType]} {''}]);
    [tDat(selectedDuplicates).excluded] = deal(1);
    expList(duplicatesIdx) = tDat;
end


%% This section section deals with short files when multiple files are detected for a mouse on same day and experiment type type
processList = cellfun(@fileparts, {expList.rawFolder}', 'uni', 0);
processList = cellfun(@(x,y) [x, y], processList, {expList.expType}', 'uni', 0);
[~, uniqueFileIdx] = unique(processList);
duplicates = unique(processList(setdiff(1:length(processList),uniqueFileIdx)));
for i = 1:length(duplicates)
    duplicatesIdx = strcmp(processList,duplicates{i}) & ~[expList.excluded]'>0;
    tDat = expList(duplicatesIdx);
    if length(tDat) < 2; continue; end
    if contains(tDat(1).expType, {'training';'inactivation'})
        fprintf('Multiple training files for %s %s %s. Keeping largest\n', tDat(1).subject, tDat(1).expDate, tDat(1).rigName);
        maxDurationIdx = num2cell([tDat.expDuration] ~= max([tDat.expDuration]));
        if all([tDat.expDuration] == max([tDat.expDuration])); maxDurationIdx{1} = 1; end
        [tDat.excluded] = maxDurationIdx{:};
        expList(duplicatesIdx) = tDat;
    end
end

%% Saving new version of expList in the shared and local directory
expList = nestedSortStruct(expList, {'subject', 'expDate'});
availableDirectories = prc.pathFinder('directoryCheck');
if availableDirectories(1); save(prc.pathFinder('dropboxlist'), 'expList'); save(strrep(prc.pathFinder('dropboxlist'),'dData','dDataLite'), 'expList'); end
if availableDirectories(1); save(prc.pathFinder('sharedlist'), 'expList'); end
end