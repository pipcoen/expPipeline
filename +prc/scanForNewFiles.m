function expList = scanForNewFiles(rebuildList, rechkExcludeTxts, chkExpRigs)
%% A funciton to check for new files for specified mice, copy them to a local directory, and update a related experimental list object
% Inputs(default values)
% rebuildList(0)---------A tag to rebuild entire experimental list from scratch.
% rechkExcludeTxts(0)----Tag to checks all directories for text files named "exclude" in all experimental folders
% chkExpRigs(1)----------Tag to checks all directories for text files named "exclude" in all experimental folders

% Outputs
% expList----------------A structure containing fields for every file detected for the selected mice in the selected daterange
%   .subject----------------Name of the mouse
%   .expDate----------------Date that the experiment was recorded
%   .expNum-----------------EXperiment number for session
%   .expDef-----------------The experimental definition file used
%   .excluded---------------Tag for excluded files (these may be excluded for a variety of reasons)
%   .rigName----------------Name of the rig where the experiment took place
%   .expType----------------Type of experiment (e.g. training, ephys, fusi, etc.)
%   .expDuration------------Duration of the experiment (s)
%   .expDets----------------Empty placeholder for extra experimental details to be later loaded

%#ok<*AGROW>
%% Define which mice should be included and which start/end dates (defaults to all files for that mouse)
%Check whether inputs exist. If not, assign default values
if ~exist('rebuildList', 'var') || isempty(rebuildList); rebuildList = 0; end
if ~exist('rechkExcludeTxts', 'var') || isempty(rechkExcludeTxts); rechkExcludeTxts = 0; end
if ~exist('chkExpRigs', 'var') || isempty(chkExpRigs); chkExpRigs = 1; end

<<<<<<< Updated upstream
expInfo = prc.pathFinder('expInfo');
includedMice = [cellfun(@(x) ['PC0' x], ...
    split({'10,11,12,13,15,22,25,27,29,30,31,32,33,34,36,37,38,41,43,45,46,48,50,51'},','), 'uni', 0); ...
        {'DJ006'; 'DJ007'; 'DJ008'; 'DJ010'; 'CR015';'AN002'}];
aliveMice = {'AN002'};

startedDates = {...
    'CR015' '2019-07-30'};
retiredDates = {...
    'PC005' '2017-07-28'; ...
    'PC006' '2017-07-28'; ...
    'PC009' '2017-03-22'; ...
    'PC014' '2017-07-01'; ...
    'DJ011' '2018-06-08'; ...
    'CR010' '2019-01-29'};
=======
expInfo = prc.pathFinder('expInfo'); %Get the top-level of the experimental files

%The full list of subjects to be included (generally, mice that learnt the task and have a decent amount of data). 
includedMice = [... 
    cellfun(@(x) ['PC0' x], split({'10,11,12,13,15,22,25,27,29,30,31,32,33,34,36,37,38,41,43,45,46,48,50,51'},','), 'uni', 0); ...
    cellfun(@(x) ['DJ0' x], split({'06,07,08,10'},','), 'uni', 0)];
aliveMice = {''}; %Mice that are currently alive (i.e. may generate new data)

%Optional "started" and "retired" dates. This could be relevant if the same subject name was used by other people, or if you wanted to exclude a
%swathe of dates from a particular mouse for some reason. Defaults assume all data from a given mouse is included.
startedDates = {'' ''};
retiredDates = {'' ''};
>>>>>>> Stashed changes
includedMice(:,2) = {datenum('2015-09-01', 'yyyy-mm-dd')};
includedMice(:,3) = {datenum(floor(now))};
[~, startedIdx] = ismember(includedMice(:,1), startedDates(:,1));
if any(startedIdx); includedMice(startedIdx > 0,2) = num2cell(datenum(startedDates(startedIdx(startedIdx>0), 2), 'yyyy-mm-dd')); end
[~, retiredIdx] = ismember(includedMice(:,1), retiredDates(:,1));
if any(retiredIdx); includedMice(retiredIdx > 0,3) = num2cell(datenum(retiredDates(retiredIdx(retiredIdx>0), 2), 'yyyy-mm-dd')); end

<<<<<<< Updated upstream
%% Check for all files generated in the past 10 days for alive mice.
nDays2Chk = 90;
=======
%% Create the search tree for experiments
%Depending on rebuildList, we assign a number of folder cycles to look through, folders from the last 10 days (for living mice), etc.
nDays2Chk = 10;
>>>>>>> Stashed changes
if rebuildList ~= 1
    cycles = 2;
    mice2Update = includedMice(contains(includedMice(:,1), aliveMice),:);
    recentDates = cell(size(mice2Update,1),1);
    for i = 1:size(mice2Update,1)
        dateRange = num2cell(datestr(datenum(mice2Update{i,3})-nDays2Chk:datenum(mice2Update{i,3}), 'yyyy-mm-dd'),2);
        recentDates{i,1} = cellfun(@(x) fileparts(fileparts(prc.pathFinder('serverfolder',mice2Update{i,1},x,'1'))),dateRange, 'uni', 0);
    end
    processList = vertcat(recentDates{:});
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
else
    cycles = 3;
    processList = [];
    for i = 1:length(expInfo)
        processList = [processList;cellfun(@(x) [expInfo{i} x], includedMice(:,1), 'uni', 0)];
    end
    expList = struct;
end

%% Look for all files in the search tree
%Using this "java.io.File" stuff because I found it to be much faster for large lists than the MATLAB alternatives.
for i = 1:cycles
    if isempty(processList); continue; end
    fprintf('Detecting folder level %d ... \n', i);
    processList = cellfun(@(x) java.io.File(x), processList, 'uni', 0);
    processList = cellfun(@(x) arrayfun(@char,x.listFiles,'uni',0), processList, 'uni', 0);
    processList = vertcat(processList{:});
    processList = processList(~contains(processList, {'Lightsheet';'ephys';'Backup'},'IgnoreCase',true));
end

%Check whether there are new files, the whole list is being rebuilt etc. If the whole list is being rebuilt, the create an empty "epxList" to append
%experiments to. Otherwise, only process identified files that are not already part of the expList.
if ~isempty(processList); blockList = processList(~cellfun(@isempty, regexp(processList, '20.*_Block.mat'))); end
if ~isempty(processList); excludeTxtFiles = processList(~cellfun(@isempty, regexp(processList, 'Exclude'))); end
dir2Exclude
filteredList = cellfun(@(x) max(blockList

if ~rebuildList && isempty(processList); fprintf('No new files found\n'); return;
elseif ~rebuildList, tLoc = prc.updatePaths(expList); processList = processList(~contains(processList,{tLoc.serverFolder}'));
end
if length(unique(cellfun(@fileparts, processList, 'uni', 0))) ~= length(processList); error('Every detected folder should be unique'); end

%% Loop to create struture for new files and add them to the expList
addedFiles  = 0; %track whether any experiments were actually added

for i = 1:length(processList)    
    %Identify details about the experiment from the folder name. Ignore the file if it is outside the data range for the subject (as defined by the
    %"started" and "retired" dates 
    splitStr = regexp(processList{i},'[\\/]','split');
    datNum = datenum(splitStr{end-2}, 'yyyy-mm-dd');
    mouseIdx = strcmp(includedMice(:,1), splitStr{end-3});
    if datNum < includedMice{mouseIdx,2} || datNum > includedMice{mouseIdx,3}; continue; end
    
    tempLoc = prc.updatePaths(expList(end));
    
    %Populate the new entry with basic subject data
    if isempty(expList); idx = 1; else, idx = length(expList)+1; end
    expList(idx).subject = splitStr{end-3};
    expList(idx).expDate = splitStr{end-2};
    expList(idx).expNum = splitStr{end-1};
    expList(idx).expDef = 'Ignored';
        
    fprintf('Adding recording %s %s %s\n', expList(end).expDate, expList(end).subject, expList(end).expNum);
    addedFiles = 1;
    %Poplate fields for addition to expList
    
    if exist(prc.pathFinder('serverblock',expList(end)), 'file'); load(prc.pathFinder('serverblock',expList(end)), 'block'); b = block;
    else, clear b; load(processList{i});
        if exist(tempLoc.rawTimeline, 'file'); load(tempLoc.rawTimeline);
        else, expList(end).excluded = 1; continue; end
        if isfield(parameters, 'experimentType') && strcmp(parameters.experimentType, 'mpep')
            b.expDef = parameters.Protocol.xfile;
            if ~isempty(dir([tempLoc.serverFolder '\*fus.mat*'])); b.rigName = 'zatteo';
            else, b.rigName = 'lilrig-stim';
            end
            b.duration = Timeline.lastTimestamp;
        else, expList(end).excluded = 1; continue;
        end
    end
        
    if isfield(b, 'expDef'); [~, expList(end).expDef] = fileparts(b.expDef); end
    expList(end).rigName = b.rigName;
    expList(end).expType = 'training';
    
    %Check and label sessions that aren't training based on the rig names
    if contains(expList(end).rigName, {'zym1'; 'zym2'}) && exist(tempLoc.galvoLog, 'file'); expList(end).expType = 'inactivation'; end
    if strcmp(expList(end).rigName, 'zatteo') && ~isempty(dir([tempLoc.serverFolder '\*fus.mat*'])); expList(end).expType = 'fusi'; end
    if strcmp(expList(end).rigName, 'zurprise') && ~isempty(dir([tempLoc.serverFolder '\*2P_00*.tif*'])); expList(end).expType = 'twophoton'; end
    if isfield(b, 'duration'); expList(end).expDuration = b.duration; else, expList(end).expDuration = 0; end
%     expList(end).blockFunction = str2func(['prc.expDef.' expList(end).expDef]);
    expList(end).excluded = 0;
    
    %If a ChoiceWorld experiment, or experiment lasts less than 60s, then ignore (Pip doesn't use Choiceworld)
    if contains(expList(end).expDef, {'Choiceworld' 'Ignored'}); expList(end).excluded = 1; continue; end
    if expList(end).expDuration < 300; expList(end).excluded = 1; continue; end
    
    txtF = [dir([fileparts(fileparts(processList{i})) '\*Exclude*.txt']); dir([fileparts(processList{i}) '\*Exclude*.txt'])];
    if ~isempty(txtF); expList(end).excluded = 1; end
    if contains(expList(end).expType, {'fusi'; 'ephys'}) && expList(end).excluded == 0; expList(end).excluded = -1; end
end
% If new files were identified, add them to expList
if ~addedFiles; fprintf('No new files found\n'); end
if ~isfield(expList, 'expDets'); expList(1).expDets = []; end
[expList(cellfun(@isempty, {expList.rigName}')).rigName] = deal('');
[expList(cellfun(@isempty, {expList.expType}')).expType] = deal('');

%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
tLoc = prc.updatePaths(expList);
if rechkExcludeTxts
    for i = 1:length(expList)
        if expList(i).excluded == 1; continue; end
        txtF = [dir([fileparts(tLoc(i).serverBlockfolder) '\*Exclude*.txt']); dir([tLoc(i).serverBlockfolder '\*Exclude*.txt'])];
        if ~isempty(txtF); expList(i).excluded = 1; end
    end
end

if chkExpRigs
    potentialEphys = find(contains({expList.rigName}', {'lilrig-stim'; 'zrig1'}) & ~contains({expList.expType}', 'ephys'));
    foundEphys = potentialEphys(arrayfun(@(x) ~isempty(dir([fileparts(x.serverFolder(1:end-1)) '\*hys*'])), tLoc(potentialEphys)));
    [expList(foundEphys).expType] = deal('ephys');
    
    potentialFusi = find(contains({expList.rigName}', 'zatteo') & ~contains({expList.expType}', 'fusi'));
    foundFusi = potentialFusi(arrayfun(@(x) exist(x.rawFusiData, 'file'), tLoc(potentialFusi))>0);
    [expList(foundFusi).expType] = deal('fusi');
end
%% Add cleanup fucntion HERE %%

%% This section section deals with cases of files on non-training rigs when multiple files are detected for a mouse on same day, rig, and expDef
processList = cellfun(@(x) fileparts(x(1:end-1)), {tLoc.serverFolder}', 'uni', 0);
processList = cellfun(@(x,y,z) [x, y, z], processList, {expList.expType}', {expList.expDef}', 'uni', 0);
processList([expList.excluded]>0) = num2cell(num2str(rand(sum([expList.excluded]>0),1),15),2);
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


%% This section section deals with short files when multiple files are detected for a mouse on same day and experiment type
processList = cellfun(@(x) fileparts(x(1:end-1)), {tLoc.serverFolder}', 'uni', 0);
processList = cellfun(@(x,y,z) [x, y, z], processList, {expList.expType}', {expList.expDef}', 'uni', 0);
processList([expList.excluded]>0) = num2cell(num2str(rand(sum([expList.excluded]>0),1),15),2);
noExcludedList = processList(~[expList.excluded]');
[~, uniqueFileIdx] = unique(noExcludedList);
duplicates = unique(noExcludedList(setdiff(1:length(noExcludedList),uniqueFileIdx)));
for i = 1:length(duplicates)
    duplicatesIdx = strcmp(processList,duplicates{i}) & [expList.excluded]'==0;
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
expList = prc.nestedSortStruct(expList, {'subject', 'expDate'});
save(prc.pathFinder('expList'), 'expList');
end