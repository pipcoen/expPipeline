function [exptList] = New_scanForNewFiles(varargin)
optInput = {0}; 
optInput(1:length(varargin)) = varargin;
[rebldAll] = optInput{:};
tDat = struct;

exptInfo = '\\zserver.cortexlab.net\Data\expInfo';
syncMatFiles;
if rebldAll; exptList = []; else, load(New_savePath('exptList')); end

lastWeek = datestr(floor(now)-6:datenum(floor(now)), 'yyyy-mm-dd');
allFiles = arrayfun(@(x) dir([exptInfo '/*PC*/' lastWeek(x,:) '/**/*block*']), 1:7, 'uni', 0)';
allFiles = cat(1,cellfun(@(x) {x.folder}', allFiles, 'uni', 0));
allFiles(cellfun(@isempty, allFiles)) = [];

excluded = {'PC001';'PC002'; 'PC003'; 'PC004'; 'PC007'; 'PC009'};
rigsList = {'zym1', 'training'; 'zym2', 'training'; 'zym3', 'training'; ...
    'zredone', 'training'; 'zredtwo', 'training'; 'zoolatry', 'training'; ...
    'zredthree', 'training'; 'zgreyfour', 'training'; ...
    'zurprise', 'twophoto';
    'zooropa', 'widfield'; 'zgood', 'widfield'};

if isempty(exptList); newFiles = dir([exptInfo '/PC*/**/*block*']);
else, [~, newIndex] = setdiff({allFiles.folder}',{exptList.foldrNam}'); 
    newFiles = allFiles(newIndex);
end

for i = 1:length(newFiles)
    load([newFiles(i).folder '/' newFiles(i).name]); 
    b = block;
    
    [mouseNam, exptDate, sessnNum] = mouseDetFromFolder(newFiles(i).folder);
    fprintf('Adding recording %s %s %s\n', exptDate, mouseNam, sessnNum);

    tDat(i,1).mouseNam = mouseNam;
    tDat(i,1).exptDate = exptDate;
    tDat(i,1).sessnNum = sessnNum;
    tDat(i,1).foldrNam = newFiles(i).folder;
    
    tDat(i,1).exRigNam = b.rigName;
    tDat(i,1).exRigTyp = rigsList{cellfun(@(x) strcmp(x, b.rigName), rigsList(:,1)),2};
    tDat(i,1).validRep = 0;
    tDat(i,1).excluded = 0;
    if isfield(b, 'duration'); tDat(i,1).duration = b.duration; end
    tDat(i,1).sessnNum = sessnNum;
    
    if isfield(b, 'expDef'); [~, tDat(i,1).exptDefn] = fileparts(b.expDef);
    else, tDat(i,1).exptDefn = 'ChoiceWorld'; tDat(i,1).excluded = 1; continue;
    end
    
    tDat(i,1).rawBlock = New_savePath('rawblock', mouseNam, exptDate, sessnNum);
    tDat(i,1).timeline = New_savePath('timeline', mouseNam, exptDate, sessnNum);
    tDat(i,1).prmsFile = New_savePath('prmsfile', mouseNam, exptDate, sessnNum);
    tDat(i,1).procData = New_savePath('procdata', mouseNam, exptDate, sessnNum);
    tDat(i,1).backData = New_savePath('backdata', mouseNam, exptDate, sessnNum);
    tDat(i,1).out2pLoc = New_savePath('out2ploc', mouseNam, exptDate, sessnNum);

    if isempty(tDat(i,1).exRigTyp); error([b.rigName ' is not a recognized rig']); end
    if contains(mouseNam, excluded); tDat(i,1).excluded = 1; continue; end
    if tDat(i,1).duration < 60; tDat(i,1).excluded = 1; continue; end

    textFile = [dir([exptInfo '/*' mouseNam '/' exptDate '/*Exclude*.txt']);...
        dir([exptInfo '/*' mouseNam '/' exptDate '/' sessnNum '/*Exclude*.txt'])];
    if ~isempty(textFile); tDat(i,1).excluded = 1; continue; end
    
    textFile = dir([exptInfo '/*' mouseNam '/' exptDate '/' sessnNum '/NewFOV*.txt']);
    if ~isempty(textFile); tDat(i,1).out2pLoc = [fileparts(exptList(i,1).out2pLoc) '/' textFile.name];
    end
    
    textFile = dir([exptInfo '/*' mouseNam '/' exptDate '/' sessnNum '/' 'ValidRepeat*.txt']);
    if ~isempty(textFile); tDat(i,1).validRep = 1; end
end
if ~isempty(newFiles); exptList = [exptList;tDat];
else; fprintf('No new files found\n');
end
%% Eliminate duplicate files for non-training rigs
allFoldr = cellfun(@fileparts, {exptList.foldrNam}', 'uni', 0);
allFoldr = cellfun(@(x,y,z) [x, y, z], allFoldr, {exptList.exRigNam}', {exptList.exptDefn}', 'uni', 0);
[~, uniqeIdx] = unique(allFoldr);
dupFoldr = allFoldr; dupFoldr(uniqeIdx) = [];
dupFoldr = unique(dupFoldr);
for i = 1:length(dupFoldr)
    dupliFil = exptList(strcmp(allFoldr,dupFoldr{i}));
    if strcmp(dupliFil(1).exRigTyp, 'training') || ...
            sum([dupliFil.excluded]>0)>=length(dupliFil)-1 || ...
            sum([dupliFil.excluded]>0)+sum([dupliFil.validRep])== length(dupliFil)
        continue;
    end
    [remvFile] = listdlg('ListString', cellfun(@num2str, {dupliFil.duration}', 'uni', 0), ...
        'PromptString', [{'Select durations to REMOVE for:'} ...
        {[dupliFil(1).mouseNam ' ' dupliFil(1).exptDate ' ' dupliFil(1).exptDefn]} {''}]);
    [dupliFil(:).validRep] = deal(1);
    [dupliFil(remvFile).excluded] = deal(1);
    [dupliFil(remvFile).validRep] = deal(0);
    exptList(strcmp(allFoldr,dupFoldr{i})) = dupliFil;
end

%% Eliminate something else
allFoldr = cellfun(@fileparts, {exptList.foldrNam}', 'uni', 0);
allFoldr = cellfun(@(x,y) [x, y], allFoldr, {exptList.exRigNam}', 'uni', 0);
[~, uniqeIdx] = unique(allFoldr);
dupFoldr = allFoldr; dupFoldr(uniqeIdx) = [];
dupFoldr = unique(dupFoldr);
for i = 1:length(dupFoldr)
    tIdx = (strcmp(allFoldr,dupFoldr{i}).*([exptList.excluded]'==0))>0;
    dupliFil = exptList(tIdx);
    if length(dupliFil) < 2; continue; end
    if strcmp(dupliFil(1).exRigTyp, 'training')
        fprintf('Multiple training files for %s %s. Keeping largest\n', ...
            dupliFil(1).mouseNam, dupliFil(1).exptDate);
        maxLenIdx = num2cell([dupliFil.duration] ~= max([dupliFil.duration]));
        [dupliFil.excluded] = maxLenIdx{:};
        exptList(tIdx) = dupliFil;
    elseif strcmp(dupliFil(1).exRigTyp, 'twophoto') && length(unique({dupliFil.out2pLoc})) > 1
        nFOV = cellfun(@isempty, (strfind({dupliFil.out2pLoc}', 'NewFOV')));
        [fileBase, fileExtn] = cellfun(@fileparts, {dupliFil.out2pLoc}', 'uni', 0);
        sExt = fileExtn(nFOV);
        sExt(1:end-1) = cellfun(@(x) [x '_'], fileExtn(1:end-1), 'uni', 0);
        [dupliFil(nFOV).out2pLoc] = deal([fileBase{1} filesep cell2mat(sExt')]);
        [dupliFil(nFOV).sessnNum] = deal({dupliFil(nFOV).sessnNum}');
        exptList(tIdx) = dupliFil;
    end
end
exptList = nestedSortStruct(exptList, {'mouseNam', 'exptDate'});
saveTarg = New_savePath('allLists'); 
if exist(saveTarg{1}, 'file'); save(saveTarg{1}, 'exptList'); end 
if exist(saveTarg{2}, 'file'); save(saveTarg{2}, 'exptList'); end 
end










