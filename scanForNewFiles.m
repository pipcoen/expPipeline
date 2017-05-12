function [expList] = scanForNewFiles(rebuildList)
if ~exist('rebuildList', 'var'); rebuildList = false; end

expInfo = pathFinder('expInfo');
if rebuildList; expList = []; else; load(pathFinder('expList')); end
lastWeek = datestr(floor(now)-6:datenum(floor(now)), 'yyyy-mm-dd');
lastBlockFiles = arrayfun(@(x) dir([expInfo '/*PC*/' lastWeek(x,:) '/**/' '*block*']), 1:7, 'uni', 0);
lastBlockFiles = cat(1,lastBlockFiles{:});

excludedMice = {'PC001';'PC002'; 'PC003'; 'PC004'; 'PC007'; 'PC009'};
rigList = {'zym1', 'training'; 'zym2', 'training'; 'zym3', 'training'; ...
    'zredone', 'training'; 'zredtwo', 'training'; 'zoolatry', 'training';
    'zredthree', 'training'; 'zgreyfour', 'training'; ...
    'zurprise', 'twophoton';
    'zooropa', 'widefield'; 'zgood', 'widefield'};
%%
if rebuildList; newBlocks = dir([expInfo '/PC*/**/*block*']);
else; [~, nIdx] = setdiff({lastBlockFiles.folder}',{expList.folderName}'); 
    newBlocks = lastBlockFiles(nIdx);
end
tDat = struct;
for i = 1:length(newBlocks)
    load([newBlocks(i).folder '/' newBlocks(i).name]); b = block;
    [subject, expDate, sessionNum] = mouseDetFromFolder(newBlocks(i).folder);
    fprintf('Adding recording %s %s %s\n', expDate, subject, sessionNum);

    tDat(i,1).subject = subject;
    tDat(i,1).expDate = expDate;
    tDat(i,1).sessionNum = sessionNum;
    tDat(i,1).eDef = 'ChoiceWorld';
    if isfield(b, 'expDef'); [~, tDat(i,1).expDef] = fileparts(b.expDef); end
    tDat(i,1).rigName = b.rigName;
    tDat(i,1).rigType = rigList{cellfun(@(x) strcmp(x, b.rigName), rigList(:,1)),2};
    if isfield(b, 'duration'); tDat(i,1).eDur = b.duration; end
    tDat(i,1).idxC = sessionNum;
    tDat(i,1).folderName = newBlocks(i).folder;
    tDat(i,1).rawBlock = pathFinder('block', subject, expDate, sessionNum);
    tDat(i,1).timeline = pathFinder('timeline', subject, expDate, sessionNum);
    tDat(i,1).paramsFile = pathFinder('parameters', subject, expDate, sessionNum);
    tDat(i,1).processedData = pathFinder('processed', subject, expDate, sessionNum);
    tDat(i,1).sharedData = pathFinder('shared', subject, expDate, sessionNum);
    tDat(i,1).suite2POutput = pathFinder('suite2POutput', subject, expDate, sessionNum);
    tDat(i,1).validRepeat = 0;
    tDat(i,1).excluded = 0;
%%
    if isempty(tDat(i,1).rTyp)
        error([b.rigName ' is not a recognized rig']);
    end
    if strcmp(tDat(i,1).eDef, 'ChoiceWorld')
        warning('Exlcuding ChoiceWorld File');
        tDat(i,1).excF = 1; continue;
    end
    if any(strcmp(subject, excM)); tDat(i,1).excF = 1; continue; end
    if tDat(i,1).eDur < 60; tDat(i,1).excF = 1; continue; end

    txtF = [dirP([eDir sep '*' subject sep expDate sep '*Exclude*.txt']);...
        dirP([eDir sep '*' subject sep expDate sep sessionNum sep '*Exclude*.txt'])];
    if ~isempty(txtF); tDat(i,1).excF = 1; continue; end
    
    txtF = dirP([eDir sep '*' subject sep expDate sep sessionNum sep 'NewFOV*.txt']);
    if ~isempty(txtF); tDat(i,1).out2 = [fileparts(expList(i,1).out2) sep txtF.name];
    end
    
    txtF = dirP([eDir sep '*' subject sep expDate sep sessionNum sep 'ValidRepeat*.txt']);
    if ~isempty(txtF); tDat(i,1).valR = 1;
    end
end
if ~isempty(newF); expList = [expList;tDat];
else; fprintf('No new files found\n');
end
%%
fLst = cellfun(@fileparts, {expList.fNam}', 'uni', 0);
fLst = cellfun(@(x,y,z) [x, y, z], fLst, {expList.rigN}', {expList.eDef}', 'uni', 0);
[~, uIdx] = unique(fLst);
dupF = fLst; dupF(uIdx) = [];
dupF = unique(dupF);
for i = 1:length(dupF)
    tIdx = strcmp(fLst,dupF{i});
    dupT = expList(tIdx);
    if strcmp(dupT(1).rTyp, 'trnR') || sum([dupT.excF]>0)>=length(dupT)-1 || ...
            sum([dupT.excF]>0)+sum([dupT.valR])== length(dupT)
        continue;
    end
    [selF] = listdlg('ListString', cellfun(@num2str, {dupT.eDur}', 'uni', 0), ...
        'PromptString', [{'Select durations to REMOVE for:'} ...
        {[dupT(1).subject ' ' dupT(1).expDate ' ' dupT(1).eDef]} {''}]);
    [dupT(:).valR] = deal(1);
    [dupT(selF).excF] = deal(1);
    [dupT(selF).valR] = deal(0);
    expList(tIdx) = dupT;
end

fLst = cellfun(@fileparts, {expList.fNam}', 'uni', 0);
fLst = cellfun(@(x,y) [x, y], fLst, {expList.rigN}', 'uni', 0);
[~, uIdx] = unique(fLst);
dupF = fLst; dupF(uIdx) = [];
dupF = unique(dupF);
for i = 1:length(dupF)
    tIdx = (strcmp(fLst,dupF{i}).*([expList.excF]'==0))>0;
    dupT = expList(tIdx);
    if length(dupT) < 2; continue; end
    if strcmp(dupT(1).rTyp, 'trnR')
        fprintf('Multiple training files for %s %s. Keeping largest\n', ...
            dupT(1).subject, dupT(1).expDate);
        mIdx = num2cell([dupT.eDur] ~= max([dupT.eDur]));
        [dupT.excF] = mIdx{:};
        expList(tIdx) = dupT;
    elseif strcmp(dupT(1).rTyp, 'twoR') && length(unique({dupT.out2})) > 1
        nFOV = cellfun(@isempty, (strfind({dupT.out2}', 'NewFOV')));
        [fBas, fExt] = cellfun(@fileparts, {dupT.out2}', 'uni', 0);
        sExt = fExt(nFOV);
        sExt(1:end-1) = cellfun(@(x) [x '_'], fExt(1:end-1), 'uni', 0);
        [dupT(nFOV).out2] = deal([fBas{1} sep cell2mat(sExt')]);
        [dupT(nFOV).idxC] = deal({dupT(nFOV).idxC}');
        expList(tIdx) = dupT;
    end
end
expList = nestedSortStruct(expList, {'subject', 'expDate'});
sLoc = pathFinder('bothexpList'); 
if exist(sLoc{1}, 'file'); save(sLoc{1}, 'expList'); end 
if exist(sLoc{2}, 'file'); save(sLoc{2}, 'expList'); end 
end










