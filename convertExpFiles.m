function [eLst] = convertExpFiles(RorS, oneM)
if ~exist('RorS', 'var') || isempty(RorS); RorS = [0 0]; end
if ~exist('oneM', 'var') || isempty(oneM); oneM = 'PC'; end

sLoc = cellfun(@(x) exist(x, 'file')>0, savePath('botheLst'));
if all(sLoc); syncfolder(savePath('labL'), savePath('savL'), 0);
elseif ~any(sLoc); error('No access to experimental directories'); 
end

[eLst] = scanForNewFiles;
if all(sLoc==[0,1]); newF = {eLst.labB}'; [eLst.finD] = newF{:}; end

for i = 1:size(eLst,1)
    if eLst(i).excF && (~exist(eLst(i).finD, 'file') || ~all(sLoc)); continue;
    elseif eLst(i).excF; delete(eLst(i).finD); 
        if exist(eLst(i).labB, 'file'); delete(eLst(i).labB); end
        continue
    end
    if isempty(strfind(eLst(i).mNam, oneM)); continue; end
    x = eLst(i); x.eLst = eLst;
    x.fSep = filesep;                     

    warning('off', 'MATLAB:load:variableNotFound');
    if exist(x.finD, 'file'); load(x.finD, 'whoD');
        if ~exist('whoD', 'var'); whoD = who('-file', x.finD);
           save(x.finD, 'whoD', '-append');
        end
    else; whoD = []; 
    end
    warning('on', 'MATLAB:load:variableNotFound');
    
    if (~any(strcmp(whoD, 'blk')) || RorS(1) == 1) && ...
            exist([x.eDef '_Blk_Proc.m'], 'file') && RorS(1)<2
        convBlockFile(x)
    end
    if strcmp(x.rTyp, 'twoR') && ...
            (~any(strcmp(whoD, 'flu')) || RorS(2) == 1) && RorS(2)<2
        conv2PData(x);
    end
    clear whoD
end
if all(sLoc); syncfolder(savePath('labL'), savePath('savL'), 0); end
end

function convBlockFile(x)
fprintf('Converting block file for %s %s\n', x.rDat,x.mNam);
b = load(x.blkR); b = b.block;

if ~strcmp(x.rTyp, 'trnR')
    t = load(x.timR); t=t.Timeline;
    b.be2t = alignBlockTimes(b, t);
end
p = load(x.parF);
[b, prm] = standardBlkNames(b, p.parameters);
x.vTri = b.events.repeatNumValues(1:length(b.events.endTrialTimes))==1;

if isfield(b.events, 'fBckValues')
    rIdx = diff([b.events.repeatNumValues(1:length(x.vTri)) 1])<0;
    x.rNum = int8([b.paramsValues(x.vTri).maxRetryIfIncorrect]>0 & ...
        b.events.fBckValues(x.vTri)<0);
    if x.rNum(end) == 1 && x.vTri(end)==1; x.rNum(end) = 0; end
    x.rNum(x.rNum>0) = b.events.repeatNumValues(rIdx)-1;
end
if isfield(b.events, 'sSrtTimes')
    blk.sSrt = single(b.events.sSrtTimes(x.vTri)');
end

blk.mNam = x.mNam;
blk.rDat = x.rDat;
blk.sNum = x.sNum;
blk.stEn = single([b.events.newTrialTimes(x.vTri)' b.events.endTrialTimes(x.vTri)']);

delP = [strfind(diff([0,b.inputs.wheelValues])~=0, [0 0]) ...
    strfind(abs(diff([0,b.inputs.wheelValues(1:2:end)]))>1, [0 0])*2];
wPos = b.inputs.wheelValues(setdiff(1:end, delP))';
wTim = b.inputs.wheelTimes(setdiff(1:end, delP))';
blk.wrTV = single([wTim wPos]); 

prm.totalTrials = length(b.events.endTrialTimes);
prm.validTrials = sum(x.vTri);
prm.minutesOnRig = round((b.experimentEndedTime-b.experimentInitTime)/60);

[blk, prm] = eval([x.eDef '_Blk_Proc(x, b, blk, prm)']);

if ~exist(fileparts(x.finD), 'dir'); mkdir(fileparts(x.finD)); end
if ~exist(x.finD, 'file'); whoD = {'blk'; 'prm'}; save(x.finD, 'blk', 'prm', 'whoD');
else; whoD = unique([who('-file', x.finD); 'blk'; 'prm']); save(x.finD, 'blk', 'prm', 'whoD', '-append'); 
end
end

function conv2PData(x)
fprintf('Converting 2P data file for %s %s\n', x.rDat,x.mNam);
if  ~exist(x.out2, 'dir') || isempty(dirP([x.out2 '*\\F*Plane*']))
    fprintf('Running Suite2P for %s %s\n', x.rDat,x.mNam);
    sFOV = str2double({x.eLst(strcmp(x.out2, {x.eLst.out2}')).sNum});
    runSuite2P(x, sFOV);
end
fLst = dirP([x.out2 '*\\F*proc.mat']);
if exist(x.out2, 'dir') && isempty(fLst)
    fprintf('NOTE: Correct 2P data for %s %s Skipping...\n', x.rDat,x.mNam);
    return;
end
t = load(savePath('rawtime', x.mNam, x.rDat, x.sNum)); t = t.Timeline;
idx = find(strcmp(x.sNum, x.idxC));
for i = 1:length(fLst)
    fprintf('Processing plane %d of %d...\n', i,length(fLst));
    load([fLst(i).folder x.fSep fLst(i).name]);
    pNum = dat.ops.iplane;
    nPLn = dat.ops.nplanes;

    nFrm = t.hw.inputs(strcmp({t.hw.inputs.name},'neuralFrames')).arrayColumn;
    nFrm = t.rawDAQData(:,nFrm);
    sTim = t.rawDAQTimestamps(diff(nFrm)==1)';
    sTim = [sTim(pNum:nPLn:end-1) sTim(pNum+1:nPLn:end)];
    sTim(size(dat.Fcell{idx},2)+1:end,:) = [];
    %%
    iCel = [dat.stat.iscell]'>0;
    sLen = cellfun(@(x) size(x,2), dat.Fcell);
    fLim = [sum(sLen(1:idx-1))+1 sum(sLen(1:idx))];
    
    [spkT, spkI] = cellfun(@sort, {dat.stat.st}', 'uni', 0);
    spkA = cellfun(@(x,y) x(y), {dat.stat.c}', spkI, 'uni', 0);
    takS = cellfun(@(x) x>=fLim(1) & x<=fLim(2), spkT, 'uni', 0);
    spkT = cellfun(@(x) x-fLim(1)+1, spkT, 'uni', 0);
    
    nCof = [dat.stat.neuropilCoefficient]';
    nCof(nCof<0.5) = 0.5; nCof(nCof > 1) = 1;
    
    celF = dat.Fcell{idx}-bsxfun(@times, (dat.FcellNeu{idx}), nCof);
    padV = 520 - size(dat.mimg(:,:,1));
    
    
    dat.stat = dat.stat(iCel);
    tmp1.frmT{i,1} = single(mean(sTim,2)');
    tmp1.skew{i,1} = [dat.stat.skew]';
    tmp1.cPnt{i,1} = uint16(round([cell2mat({dat.stat.med}') [dat.stat.iplane]']));
    tmp1.mImg{i,1} = permute(padarray(single(dat.mimg(:,:,2)), padV, nan, 'post'), [3,1,2]);
    tmp1.pImg{i,1} = permute(padarray(single(dat.mimg(:,:,5)), padV, nan, 'post'), [3,1,2]);
    tmp1.celF{i,1} = celF(iCel,:);
    tmp1.nCof{i,1} = single(nCof(iCel));
    
    tmp2.spkT{i,1} = cellfun(@(x,y) x(y), spkT(iCel), takS(iCel), 'uni', 0);
    tmp2.spkA{i,1} = cellfun(@(x,y) x(y), spkA(iCel), takS(iCel), 'uni', 0);
    
    clear dat;
end
minF = min(cellfun(@length, tmp1.celF));
redL = @(x) cellfun(@(y) y(:,1:min([size(y,2) minF]),:),x, 'uni',0);
tmp1 = structfun(@(x) redL(x), tmp1, 'uni', 0);

flu = structfun(@cell2mat, tmp1, 'uni', 0);
flu.mImg = permute(flu.mImg, [2,3,1]);
flu.pImg = permute(flu.pImg, [2,3,1]);

flu.spkA = fun.map(@(x,y) x(y<=minF), cellflat(tmp2.spkA), cellflat(tmp2.spkT));
flu.spkT = fun.map(@(x) x(x<=minF), cellflat(tmp2.spkT));
flu.mNam = x.mNam;
flu.rDat = x.rDat;
flu.sNam = x.sNum;

if ~exist(fileparts(x.finD), 'dir'); mkdir(fileparts(x.finD)); end
if ~exist(x.finD, 'file'); whoD = {'flu'}; save(x.finD, 'flu', 'whoD');
else; whoD = [who('-file', x.finD); 'flu']; save(x.finD, 'flu', 'whoD', '-append'); 
end
end

function runSuite2P(x, expts)
db(1).mouse_name    = x.mNam;   db(1).date          = x.rDat;
db(1).expts         = expts;    db(1).nchannels     = 1;
db(1).gchannel      = 1;        db(1).nplanes       = 4;
db(1).expred        = [];       db(1).nchannels_red = [];
db(1).comments      = '';       db(1).planesToProcess = [1,2,3,4];
master_file_Pip;
end

