function expList = scanForNewFiles(rebuildList, chkExpRigs)
%% A funciton to check for new files for specified mice, copy them to a local directory, and update a related experimental list object
% NOTE: I am using this version of this code to only process signals files. mpep experiments are ignored.
% Inputs(default values)
% rebuildList(0)---------A tag to rebuild entire experimental list from scratch (if =2, also deletes all autoExcluded text files).
% chkExpRigs(1)----------Tag to checks all directories for text files named "exclude" in all experimental folders

% Outputs
% expList----------------A structure containing fields for every file detected for the selected mice in the selected daterange
%   .subject----------------Name of the mouse
%   .expDate----------------Date that the experiment was recorded
%   .expNum-----------------Experiment number for session
%   .expDef-----------------The experimental definition file used
%   .excluded---------------Tag for excluded files (these may be excluded for a variety of reasons)
%   .rigName----------------Name of the rig where the experiment took placeexpList
%   .expType----------------Type of experiment (e.g. training, ephys, fusi, etc.)
%   .expDuration------------Duration of the experiment (s)
%   .expDets----------------Empty placeholder for extra experimental details to be later loaded

%#ok<*AGROW>
%% Define which mice should be included and which start/end dates (defaults to all files for that mouse)
%Check whether inputs exist. If not, assign default values
if ~exist('rebuildList', 'var') || isempty(rebuildList); rebuildList = 0; end
if ~exist('chkExpRigs', 'var') || isempty(chkExpRigs); chkExpRigs = 1; end

expInfo = prc.pathFinder('expInfo');

%The full list of subjects to be included (generally, mice that learnt the task and have a decent amount of data). 
includedMice = [... 
    cellfun(@(x) ['PC0' x], split({'10,11,12,13,15,22,27,29,30,31,32,33,34,43,45,46,48,50,51,52,53,54,55'},','), 'uni', 0); ...
    cellfun(@(x) ['DJ0' x], split({'06,07,08,10'},','), 'uni', 0)];
aliveMice = {'X'}; %Mice that are currently alive (i.e. may generate new data)

%Optional "started" and "retired" dates. This could be relevant if the same subject name was used by other people, or if you wanted to exclude a
%swathe of dates from a particular mouse for some reason. Defaults assume all data from a given mouse is included.
startedDates = {'' ''};
retiredDates = {'' ''};
includedMice(:,2) = {datenum('2015-09-01', 'yyyy-mm-dd')};
includedMice(:,3) = {datenum(floor(now))};
[~, startedIdx] = ismember(includedMice(:,1), startedDates(:,1));
if any(startedIdx); includedMice(startedIdx > 0,2) = num2cell(datenum(startedDates(startedIdx(startedIdx>0), 2), 'yyyy-mm-dd')); end
[~, retiredIdx] = ismember(includedMice(:,1), retiredDates(:,1));
if any(retiredIdx); includedMice(retiredIdx > 0,3) = num2cell(datenum(retiredDates(retiredIdx(retiredIdx>0), 2), 'yyyy-mm-dd')); end

%% Create the search tree for experiments
%Depending on rebuildList, we assign a number of folder cycles to look through, folders from the last 10 days (for living mice), etc.
nDays2Chk = 2000;
if rebuildList == 0
    cycles = 2;
    mice2Update = includedMice(contains(includedMice(:,1), aliveMice),:);
    recentDates = cell(size(mice2Update,1),1);
    for i = 1:size(mice2Update,1)
        pathInfo.expDate = num2cell(datestr(datenum(mice2Update{i,3})-nDays2Chk:datenum(mice2Update{i,3}), 'yyyy-mm-dd'),2);
        pathInfo.datNum = cellfun(@datenum, pathInfo.expDate, 'uni', 0);
        pathInfo.subject = repmat(mice2Update(i,1), length(pathInfo.expDate),1);
        pathInfo.expNum = repmat({'1'}, length(pathInfo.expDate),1);
        recentDates{i,1} = cellfun(@(x) fileparts(fileparts(x)),prc.pathFinder('serverfolder', pathInfo), 'uni', 0);
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
excludeTxtFiles = [];
for i = 1:cycles
    if isempty(processList); continue; end
    fprintf('Detecting folder level %d ... \n', i);
    processList = cellfun(@(x) java.io.File(x), processList, 'uni', 0);
    processList = cellfun(@(x) arrayfun(@char,x.listFiles,'uni',0), processList, 'uni', 0);
    processList = vertcat(processList{:});
    processList = processList(~contains(processList, {'Lightsheet';'ephys';'Backup'},'IgnoreCase',true));
    excludeTxtFiles = [excludeTxtFiles; processList(~cellfun(@isempty, regexp(processList, 'Exclude.txt')))];
end

%CAREFUL--commented so not accidentally used. Delete the autoExcluded experiments if completely rebuilding the expList (basically deletes a bunch 
%of text files)
% if rebuildList == 2
%     cellfun(@delete, excludeTxtFiles(contains(excludeTxtFiles, 'auto2020PaperExclude.txt')));
%     excludeTxtFiles = excludeTxtFiles(~contains(excludeTxtFiles, 'auto2020PaperExclude.txt'));
% end

if ~isempty(processList)
    dir2Exclude = cellfun(@fileparts, excludeTxtFiles, 'uni', 0);
    processList = processList(~cellfun(@isempty, regexp(processList, '20.*_Block.mat'))); 
    processList(contains(processList, dir2Exclude)) = [];
end

%Check whether there are new files, the whole list is being rebuilt etc. If the whole list is being rebuilt, the create an empty "epxList" to append
%experiments to. Otherwise, only process identified files that are not already part of the expList.
oddNameCorrection = strcmp({expList.expDef}', 'multiSpaceWorldNewNames');
if any(oddNameCorrection); [expList(oddNameCorrection).expDef] = deal('multiSpaceWorld'); save(prc.pathFinder('expList'), 'expList'); end
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
    
    %Populate the new entry with basic subject data
    addedFiles = 1;
    newExp = struct;
    newExp.subject = splitStr{end-3};
    newExp.expDate = splitStr{end-2};
    newExp.expNum = splitStr{end-1};
    tempLoc = prc.updatePaths(newExp);
    load(prc.pathFinder('serverblock',newExp), 'block'); 
    b = block;
    
    backUpDir = fileparts(prc.pathFinder('backupblock', newExp));
    if strcmp(hostname, 'zippy') && ~exist(prc.pathFinder('backupblock', newExp), 'file')
        if ~exist(backUpDir, 'dir'); mkdir(backUpDir); end
        copyfile(prc.pathFinder('serverblock', newExp), backUpDir)
        copyfile(prc.pathFinder('serverparams', newExp), backUpDir)
    end
    
    %Notify which files are bing added and load block file
    fprintf('Adding recording %s %s %s: %d of %d\n', newExp.expDate, newExp.subject, newExp.expNum, i, length(processList));
    
    %Add rig name and expType
    [~, newExp.expDef] = fileparts(b.expDef);
    newExp.rigName = b.rigName;
    newExp.expType = 'training';

    %Check the session is an inactivation session based on the rig names and existence of a galvolog file
    if contains(newExp.rigName, {'zym1'; 'zym2'}) && exist(tempLoc.galvoLog, 'file')
        %If galvoLog isn't a struct, or only has one field ("stereotaxCalib") then inactivation didn't run, so it is a "training" session
        galvoLog = load(tempLoc.galvoLog);
        if isstruct(galvoLog) && length(fields(galvoLog))>1; newExp.expType = 'inactivation';; end 
    end
    
    %Check whether to autoExclude for paper (based on whether it's the "wrong" exp type, or whether it is too short etc.) If so, created a text file
    %with the name 'auto2020PaperExclude.txt' that contains the reason for exclusion.
    if ~isfield(b, 'duration') || b.duration < 300
        fid = fopen([tempLoc.serverFolder 'auto2020PaperExclude.txt'],'wt');
        fprintf(fid, 'Duration undetected or less than 5 minutes, so assumed an erroneous');
        fclose(fid);
        continue;
    else, newExp.expDuration = b.duration;
    end
    
    if ~contains(newExp.expDef, 'SpaceWorld')
        fid = fopen([tempLoc.serverFolder 'auto2020PaperExclude.txt'],'wt');
        fprintf(fid, 'Detected non-SpaceWorldExperiment so not relevant for 2020 paper');
        fclose(fid);
        continue;
    end
    if strcmp(newExp.rigName, 'zatteo') && ~isempty(dir([tempLoc.serverFolder '\*fus.mat*']))
        fid = fopen([tempLoc.serverFolder 'auto2020PaperExclude.txt'],'wt');
        fprintf(fid, 'Detected as fusi exp so not relevant for 2020 paper');
        fclose(fid);
        continue;
    end
    if strcmp(newExp.rigName, 'zurprise') && ~isempty(dir([tempLoc.serverFolder '\*2P_00*.tif*'])) 
        fid = fopen([tempLoc.serverFolder 'auto2020PaperExclude.txt'],'wt');
        fprintf(fid, 'Detected as 2P exp so not relevant for 2020 paper');
        fclose(fid);
        continue;
    end
    if contains(newExp.expDef, {'Choiceworld'}) 
        fid = fopen([tempLoc.serverFolder 'auto2020PaperExclude.txt'],'wt');
        fprintf(fid, 'Detected as choiceworld experiment so not relevant for 2020 paper');
        fclose(fid);
        continue;
    end

    newExp.expDets = [];
    if isempty(fields(expList)); expList = newExp; else, expList(end+1) = newExp; end
end

% If new files were identified, add them to expList. Also add empty field "expDets" which is used for future details (i.e. probe info)
if ~addedFiles; fprintf('No new files found\n'); end
if ~isfield(expList, 'expDets'); expList(1).expDets = []; end


%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
tLoc = prc.updatePaths(expList);
if chkExpRigs
    potentialEphys = find(contains({expList.rigName}', {'lilrig-stim'; 'zrig1'}) & ~contains({expList.expType}', 'ephys'));
    foundEphys = potentialEphys(arrayfun(@(x) ~isempty(dir([fileparts(x.serverFolder(1:end-1)) '\*hys*'])), tLoc(potentialEphys)));
    [expList(foundEphys).expType] = deal('ephys');
end

%% This section section deals with short files when multiple files are detected for a mouse on same day and experiment type. Only keeps longest file.
processList = cellfun(@(x) fileparts(x(1:end-1)), {tLoc.serverFolder}', 'uni', 0);
processList = cellfun(@(x,y,z) [x, y, z], processList, {expList.expType}', {expList.expDef}', 'uni', 0);
[~, uniqueFileIdx] = unique(processList);
duplicates = unique(processList(setdiff(1:length(processList),uniqueFileIdx)));
idx2Remove = zeros(length(processList),1, 'logical');
for i = 1:length(duplicates)
    duplicatesIdx = find(strcmp(processList,duplicates{i}));
    tDat = expList(duplicatesIdx);
    tDat = prc.updatePaths(tDat);
    if contains(tDat(1).expType, {'training';'inactivation'})
        fprintf('Multiple training files for %s %s %s. Keeping largest\n', tDat(1).subject, tDat(1).expDate, tDat(1).rigName);
        maxDurationIdx = num2cell([tDat.expDuration] == max([tDat.expDuration]));
        if all([tDat.expDuration] == max([tDat.expDuration])); maxDurationIdx{1} = 1; end
        
        for j = 1:length(maxDurationIdx)
            if maxDurationIdx{j}==1; continue; end
            fid = fopen([tDat(j).serverFolder 'auto2020PaperExclude.txt'],'wt');
            fprintf(fid, 'Multiple training files for same day, so exluded shorter ones for simplicity');
            fclose(fid);
            idx2Remove(duplicatesIdx(j)) = 1;
        end
    end
end
expList = expList(~idx2Remove);

%% Saving new version of expList in the shared and local directory
expList = prc.nestedSortStruct(expList, {'subject', 'expDate'});
save(prc.pathFinder('expList'), 'expList');
end