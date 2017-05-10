function [eLst] = scanForNewFiles(rBld)
%%
eDir = '\\zserver.cortexlab.net\Data\expInfo';
syncMatFiles;
sep = filesep;
if exist('rBld', 'var') && rBld; eLst = [];
else; load(savePath('eLst'));
end
lasW = datestr(floor(now)-7:datenum(floor(now)), 'yyyy-mm-dd');
allF = arrayfun(@(x) dirP([eDir sep '*PC*' sep lasW(x,:) sep '**' sep '*block*']), ...
    (1:size(lasW,1))', 'uni', 0);
allF(cellfun(@isempty, allF)) = [];
allF = cat(1,allF{:});

excM = {'PC001';'PC002'; 'PC003'; 'PC004'; 'PC007'; 'PC009'};

rLst = {'zym1', 'trnR'; 'zym2', 'trnR'; 'zym3', 'trnR'; ...
    'zredone', 'trnR'; 'zredtwo', 'trnR'; 'zoolatry', 'trnR';
    'zredthree', 'trnR'; 'zgreyfour', 'trnR'; ...
    'zurprise', 'twoR';
    'zooropa', 'widR'; 'zgood', 'widR'};

if isempty(eLst); newF = dirP([eDir sep 'PC*' sep '**' sep '*block*']);
else; [~, nIdx] = setdiff({allF.folder}',{eLst.fNam}'); newF = allF(nIdx);
end
tDat = struct;
for i = 1:length(newF)
    load([newF(i).folder sep newF(i).name]); b = block;
    [mNam, rDat, sNum] = mouseDetFromFolder(newF(i).folder);
    fprintf('Adding recording %s %s %s\n', rDat, mNam, sNum);

    tDat(i,1).mNam = mNam;
    tDat(i,1).rDat = rDat;
    tDat(i,1).sNum = sNum;
    tDat(i,1).fNam = newF(i).folder;
    tDat(i,1).eDef = 'ChoiceWorld';
    if isfield(b, 'expDef'); [~, tDat(i,1).eDef] = fileparts(b.expDef); end
    tDat(i,1).rigN = b.rigName;
    tDat(i,1).rTyp = rLst{cellfun(@(x) strcmp(x, b.rigName), rLst(:,1)),2};
    tDat(i,1).valR = 0;
    tDat(i,1).excF = 0;
    if isfield(b, 'duration'); tDat(i,1).eDur = b.duration; end
    tDat(i,1).idxC = sNum;
    tDat(i,1).blkR = savePath('rawblock', mNam, rDat, sNum);
    tDat(i,1).timR = savePath('rawtime', mNam, rDat, sNum);
    tDat(i,1).parF = savePath('params', mNam, rDat, sNum);
    tDat(i,1).finD = savePath('findata', mNam, rDat, sNum);
    tDat(i,1).labB = savePath('labback', mNam, rDat, sNum);
    tDat(i,1).out2 = savePath('out2psuite', mNam, rDat, sNum);

    if isempty(tDat(i,1).rTyp)
        error([b.rigName ' is not a recognized rig']);
    end
    if strcmp(tDat(i,1).eDef, 'ChoiceWorld')
        warning('Exlcuding ChoiceWorld File');
        tDat(i,1).excF = 1; continue;
    end
    if any(strcmp(mNam, excM)); tDat(i,1).excF = 1; continue; end
    if tDat(i,1).eDur < 60; tDat(i,1).excF = 1; continue; end

    txtF = [dirP([eDir sep '*' mNam sep rDat sep '*Exclude*.txt']);...
        dirP([eDir sep '*' mNam sep rDat sep sNum sep '*Exclude*.txt'])];
    if ~isempty(txtF); tDat(i,1).excF = 1; continue; end
    
    txtF = dirP([eDir sep '*' mNam sep rDat sep sNum sep 'NewFOV*.txt']);
    if ~isempty(txtF); tDat(i,1).out2 = [fileparts(eLst(i,1).out2) sep txtF.name];
    end
    
    txtF = dirP([eDir sep '*' mNam sep rDat sep sNum sep 'ValidRepeat*.txt']);
    if ~isempty(txtF); tDat(i,1).valR = 1;
    end
end
if ~isempty(newF); eLst = [eLst;tDat];
else; fprintf('No new files found\n');
end
%%
fLst = cellfun(@fileparts, {eLst.fNam}', 'uni', 0);
fLst = cellfun(@(x,y,z) [x, y, z], fLst, {eLst.rigN}', {eLst.eDef}', 'uni', 0);
[~, uIdx] = unique(fLst);
dupF = fLst; dupF(uIdx) = [];
dupF = unique(dupF);
for i = 1:length(dupF)
    tIdx = strcmp(fLst,dupF{i});
    dupT = eLst(tIdx);
    if strcmp(dupT(1).rTyp, 'trnR') || sum([dupT.excF]>0)>=length(dupT)-1 || ...
            sum([dupT.excF]>0)+sum([dupT.valR])== length(dupT)
        continue;
    end
    [selF] = listdlg('ListString', cellfun(@num2str, {dupT.eDur}', 'uni', 0), ...
        'PromptString', [{'Select durations to REMOVE for:'} ...
        {[dupT(1).mNam ' ' dupT(1).rDat ' ' dupT(1).eDef]} {''}]);
    [dupT(:).valR] = deal(1);
    [dupT(selF).excF] = deal(1);
    [dupT(selF).valR] = deal(0);
    eLst(tIdx) = dupT;
end

fLst = cellfun(@fileparts, {eLst.fNam}', 'uni', 0);
fLst = cellfun(@(x,y) [x, y], fLst, {eLst.rigN}', 'uni', 0);
[~, uIdx] = unique(fLst);
dupF = fLst; dupF(uIdx) = [];
dupF = unique(dupF);
for i = 1:length(dupF)
    tIdx = (strcmp(fLst,dupF{i}).*([eLst.excF]'==0))>0;
    dupT = eLst(tIdx);
    if length(dupT) < 2; continue; end
    if strcmp(dupT(1).rTyp, 'trnR')
        fprintf('Multiple training files for %s %s. Keeping largest\n', ...
            dupT(1).mNam, dupT(1).rDat);
        mIdx = num2cell([dupT.eDur] ~= max([dupT.eDur]));
        [dupT.excF] = mIdx{:};
        eLst(tIdx) = dupT;
    elseif strcmp(dupT(1).rTyp, 'twoR') && length(unique({dupT.out2})) > 1
        nFOV = cellfun(@isempty, (strfind({dupT.out2}', 'NewFOV')));
        [fBas, fExt] = cellfun(@fileparts, {dupT.out2}', 'uni', 0);
        sExt = fExt(nFOV);
        sExt(1:end-1) = cellfun(@(x) [x '_'], fExt(1:end-1), 'uni', 0);
        [dupT(nFOV).out2] = deal([fBas{1} sep cell2mat(sExt')]);
        [dupT(nFOV).idxC] = deal({dupT(nFOV).idxC}');
        eLst(tIdx) = dupT;
    end
end
eLst = nestedSortStruct(eLst, {'mNam', 'rDat'});
sLoc = savePath('botheLst'); 
if exist(sLoc{1}, 'file'); save(sLoc{1}, 'eLst'); end 
if exist(sLoc{2}, 'file'); save(sLoc{2}, 'eLst'); end 
end










