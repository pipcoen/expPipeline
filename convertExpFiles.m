function convertExpFiles(redoBlocks, redoSuite2P, selectedMice)
if ~exist('redoBlocks', 'var') || isempty(redoBlocks); redoBlocks = 0; end
if ~exist('redoSuite2P', 'var') || isempty(redoSuite2P); redoSuite2P = 0; end
if ~exist('selectedMice', 'var') || isempty(redoSuite2P); selectedMice = {'PC'}; end

existDirectories = pathFinder('directoryCheck'); 
if all(existDirectories); syncfolder(pathFinder('sharedFolder'), pathFinder('processedFolder'), 0); end

expList = scanForNewFiles;
if all(existDirectories==[0,1])
    newPathList = {expList.sharedData}'; 
    [expList.processedData] = newPathList{:}; 
end
processedFiles = [{expList.processedData}', {expList.sharedData}'];
existProcessed = cellfun(@(x) exist(x,'file'), processedFiles)>0;
listNotExcluded = ~[expList.excluded]';
if ~any([redoBlocks redoSuite2P])
    processList = find(listNotExcluded & ~all(existProcessed, 2));
else, processList = find(listNotExcluded | any(existProcessed, 2));
end

for i = processList'
    if ~contains(expList(i).subject, selectedMice); continue; end
    if contains(expList(i).subject, {'PC008'}); continue; end
    x = expList(i); x.expList = expList;
   
    x.existProcessed = existProcessed(i,:); 
    if expList(i).excluded; cellfun(@delete, processedFiles(x.existProcessed)); continue; end
    
    warning('off', 'MATLAB:load:variableNotFound'); 
    if x.existProcessed(1)
        load(x.processedData, 'whoD');
        if ~exist('whoD', 'var')
            whoD = who('-file', x.processedData);
           save(x.processedData, 'whoD', '-append');
        end
    else; whoD = [];
    end
    warning('on', 'MATLAB:load:variableNotFound');
    
    varIdx = contains({'blk', 'flu'}, ['ignore'; whoD]);
    
    if any([~varIdx(2) redoSuite2P])
        switch lower(x.rigNameType{2})
            case 'twophoton'; conv2PData(x);            
            case 'widefield'; continue; %convWidefieldData(x);
        end
    end
    
    if ~exist(x.blockHelper, 'file'); continue; end
    if any([~varIdx(1) redoBlocks]); convBlockFile(x); end
end
if all(existDirectories); syncfolder(pathFinder('sharedFolder'), pathFinder('processedFolder'), 0); end
end

function convBlockFile(x)
fprintf('Converting block file for %s %s\n', x.expDate,x.subject);
x.oldBlock = load(x.rawBlock); x.oldBlock = x.oldBlock.block;
x.oldParams = load(x.rawParams); x.oldParams = x.oldParams.parameters;

if ~strcmp(x.rigNameType{2}, 'training')
    x.timeline = load(x.rawTimeline); x.timeline=x.timeline.Timeline;
    x.oldBlock.blockTimeOffset = alignBlockTimes(x.oldBlock, x.timeline);
end
[x.standardizedBlock, x.standardizedParams] = standardBlkNames(x.oldBlock, x.oldParams);

x.validTrials = x.standardizedBlock.events.repeatNumValues(1:length(x.standardizedBlock.events.endTrialTimes))==1;

if isfield(x.standardizedBlock.events, 'feedbackValues')
    repeatIdx = diff([x.standardizedBlock.events.repeatNumValues(1:length(x.validTrials)) 1])<0;
    x.repeatNum = int8([x.standardizedBlock.paramsValues(x.validTrials).maxRepeatIncorrect]>0 ...
        & x.standardizedBlock.events.feedbackValues(x.validTrials)<0);
    if x.repeatNum(end) == 1 && x.validTrials(end)==1; x.repeatNum(end) = 0; end
    x.repeatNum(x.repeatNum>0) = x.standardizedBlock.events.repeatNumValues(repeatIdx)-1;
end

x.newBlock.subject = x.subject;
x.newBlock.expDate = x.expDate;
x.newBlock.sessionNum = x.sessionNum;
x.newBlock.trialStart = single(x.standardizedBlock.events.newTrialTimes(x.validTrials)');
x.newBlock.trialEnd = x.standardizedBlock.events.endTrialTimes(x.validTrials)';

if isfield(x.standardizedBlock.events, 'stimPeriodOnOffTimes')
    stimPeriodStart = x.standardizedBlock.events.stimPeriodOnOffTimes(x.standardizedBlock.events.stimPeriodOnOffValues == 1)';
    stimPeriodStart = indexByTrial(x.newBlock, stimPeriodStart, stimPeriodStart, 0);
    x.newBlock.stimPeriodStart = single(cell2mat(stimPeriodStart));
end

repeatPoints = [strfind(diff([0,x.standardizedBlock.inputs.wheelValues])~=0, [0 0]) ...
    strfind(abs(diff([0,x.standardizedBlock.inputs.wheelValues(1:2:end)]))>1, [0 0])*2];
wheelValue = x.standardizedBlock.inputs.wheelValues(setdiff(1:end, repeatPoints))';
wheelTime = x.standardizedBlock.inputs.wheelTimes(setdiff(1:end, repeatPoints))';
x.newBlock.rawWheelTimeValue = single([wheelTime wheelValue]); 

x.standardizedParams.totalTrials = length(x.standardizedBlock.events.endTrialTimes);
x.standardizedParams.validTrials = sum(x.validTrials);
x.standardizedParams.minutesOnRig = round((x.standardizedBlock.experimentEndedTime-x.standardizedBlock.experimentInitTime)/60);

x.blockFunciton = str2func(x.blockHelper(1:end-2));

[blk, prm] = x.blockFunciton(x); %#ok
if ~exist(fileparts(x.processedData), 'dir'); mkdir(fileparts(x.processedData)); end
if ~x.existProcessed(1); whoD = {'blk'; 'prm'}; save(x.processedData, 'blk', 'prm', 'whoD'); %#ok
else; whoD = unique([who('-file', x.processedData); 'blk'; 'prm']); save(x.processedData, 'blk', 'prm', 'whoD', '-append'); %#ok
end
end

function conv2PData(x)
fprintf('Converting 2P data file for %s %s\n', x.expDate,x.subject);
if  ~exist(x.out2, 'dir') || isempty(dirP([x.out2 '*\\F*Plane*']))
    fprintf('Running Suite2P for %s %s\n', x.expDate,x.subject);
    sFOV = str2double({x.expList(strcmp(x.out2, {x.expList.out2}')).sessionNum});
    runSuite2P(x, sFOV);
end
fLst = dirP([x.out2 '*\\F*proc.mat']);
if exist(x.out2, 'dir') && isempty(fLst)
    fprintf('NOTE: Correct 2P data for %s %s Skipping...\n', x.expDate,x.subject);
    return;
end
t = load(pathFinder('rawtime', x.subject, x.expDate, x.sessionNum)); t = t.Timeline;
idx = find(strcmp(x.sessionNum, x.idxC));
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
flu.subject = x.subject;
flu.expDate = x.expDate;
flu.sNam = x.sessionNum;

if ~exist(fileparts(x.processedData), 'dir'); mkdir(fileparts(x.processedData)); end
if ~exist(x.processedData, 'file'); whoD = {'flu'}; save(x.processedData, 'flu', 'whoD');
else; whoD = [who('-file', x.processedData); 'flu']; save(x.processedData, 'flu', 'whoD', '-append'); 
end
end

function runSuite2P(x, expts)
db(1).mouse_name    = x.subject;   db(1).date          = x.expDate;
db(1).expts         = expts;    db(1).nchannels     = 1;
db(1).gchannel      = 1;        db(1).nplanes       = 4;
db(1).expred        = [];       db(1).nchannels_red = [];
db(1).comments      = '';       db(1).planesToProcess = [1,2,3,4];
master_file_Pip;
end

