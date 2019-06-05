function convertExpFiles(redoTag, dataType, selectedSubjects, selectedDates)
%% A funciton process experimental recodings into more concise, matching, local copies for future analysis.

%Inputs(default values)
%instructBlk(0)------EITHER, a single RAW block file to process (must be that block's directory) A tag to redo all block files
%instruct2P(0)-----A tag to redo all 2P files
%filterOptions('0')-Run for a specific mouse, or group of mice containing this string

%Outputs
%An output files is generated for each experiment in the form subject_yymmdd_expNumProc.mat These files will contain the following
%"blk" is a structure comprising a reduced, processed block file which contains all essential information. Fields common to all experiments are:
%.subject-----------Name of the mouse
%.expDate-----------Date that the experiment was recorded
%.expNum------------Session number for experiment
%.rigName-----------Name of the rig where the experiment took place
%.expType-----------Type of the rig where the experiment took place
%.trialStart--------nx1 vector of trial start times relative to the start of the experiment (s)
%.trialEnd----------nx1 vector of trial end times relative to the start of the experiment (s)
%.????????----------Additional fields are specified in the helper function for each experimental definition

%"prm" is a structure comprising a reduced, processed parameters file which contains all essential information. Fields common to all experiments are:
%.subject-----------Name of the mouse
%.expDate-----------Date that the experiment was recorded
%.expNum------------Session number for experiment
%.rigName-----------Name of the rig where the experiment took place
%.expType-----------Type of the rig where the experiment took place
%.minutesOnRig------Number of minutes spent on the rig
%.numRepeats--------1xn vector of (max) number of repeats for each parameter conditions.
%.????????----------Additional fields are specified in the helper function for each experimental definition

%"raw" is a structure comprising potentially useful (realtively) raw data (wheel movement etc.) which is not used for a lot of analyses and
%so should only be loaded if necessary (as it is large).

%"whoD" is simply a list of which variables are in the file. It is much faster to load this when processing rather than check the file contents.

%% Set default values, load experimental list, check which processed files already exist, etc.
if ~exist('redoTag', 'var') || isempty(redoTag); redoTag = 0; end
if ~exist('dataType', 'var') || isempty(dataType); dataType = 'all'; end
if ~exist('selectedSubjects', 'var'); selectedSubjects = {'PC'; 'DJ'}; end
if ~exist('selectedDates', 'var'); selectedDates = {'1'}; end
if strcmp(hostname, 'zip'); zipComp = 1; else, zipComp = 0; end

if zipComp; fprintf('Running on Zip so will sync local folder and server... \n');
    prc.syncfolder(prc.pathFinder('processedDirectory'), prc.pathFinder('serverProcessedDirectory'), 0);
end

%Scan for new files. Then, if the only available directory is zserver, change the save destination to be zserver
expList = prc.scanForNewFiles(0,[0 1]);
expList = prc.updatePaths(expList);
expList = expList(cellfun(@(x) contains(x, selectedSubjects), {expList.subject}));
expList = expList(cellfun(@(x) contains(x, selectedDates), {expList.expDate}));

expDefs2Run = unique({expList.expDef}');
expDefs2Remove = expDefs2Run(cellfun(@(x) isempty(which(['prc.expDef.' x])), expDefs2Run));
if ~isempty(expDefs2Remove)
    fprintf('Warning: the following expDefs will be skipped: \n');
    fprintf('-%s \n', expDefs2Remove{:});
    expList = expList(~contains({expList.expDef}', [expDefs2Remove; 'Temporal']));
end

existProcessed = cellfun(@(x) exist(x,'file'), {expList.processedData}')>0;      %Logical of whether processed files exist for those directories
%Create processList and deleteLise. If there are instructions to redo the analysis in the inputs, process list will be all the unexcluded paths,
%otherwise it will only include the paths of files that don't already exist. deleteList is for files that are excluded, but do exist (these were
%probably processed and then excluded at a later date)
if zipComp
    deleteList = find([expList.excluded]'==1 & existProcessed);
    cellfun(@delete, {expList(deleteList).processedData});
    cellfun(@delete, {expList(deleteList).serverProcessedData});
    cellfun(@(x) delete(strrep(x, 'Timeline', 'ProcBlock')), {expList(deleteList).rawTimeline});
    existProcessed = existProcessed([expList.excluded]'~=1);
    expList = expList([expList.excluded]'~=1);

    missBackup = find(~cellfun(@(x) exist(x,'file'), {expList.rawBlock}')>0);
    if ~isempty(missBackup)
        missFolders = cellfun(@fileparts, {expList(missBackup).rawBlock}', 'uni', 0);
        cellfun(@(x) mkdir(x), missFolders(cellfun(@(x) ~exist(x, 'dir'), missFolders)));
        arrayfun(@(x,y) copyfile(prc.pathFinder('serverBlock', expList(x)), y{1}), missBackup, {expList(missBackup).rawBlock}');
        arrayfun(@(x,y) copyfile(prc.pathFinder('serverParams', expList(x)), y{1}), missBackup, {expList(missBackup).rawParams}');
    end
else
    existProcessed = existProcessed([expList.excluded]'~=1);
    expList = expList([expList.excluded]'~=1);
end

fusiIdx = find(contains({expList.expType}', 'fusi'));
kruminBackup = cellfun(@(x) strrep(x, 'Timeline', 'ProcBlock'), {expList.rawTimeline}', 'uni', 0); 
copyKruminBackup = fusiIdx(~cellfun(@(x) exist(x,'file'), kruminBackup(fusiIdx))>0 & existProcessed(fusiIdx));
arrayfun(@(x) copyfile(expList(x).processedData, kruminBackup{x}), copyKruminBackup);

if redoTag==1; files2Run = find(existProcessed*0+1)';
else, files2Run = find(~existProcessed | contains({expList.expType}', {'fusi';'ephys'}))';
end

%% Loop to process the files
%srtIdx can be more than zero if one wants to redo all files, but start in the midded. Useful if MATLAB crashes
srtIdx = 0;
for i = files2Run(files2Run>srtIdx)
    %create x, the processing structure. It is the expList entry for the current index, but also contains the entire list
    x = expList(i); x.expList = expList;
    x.existProcessed = existProcessed(i,:);
    
    %If the processed file exists, try to load the whoD variable (turn off warning if not found) which contains a list of the variable in that
    %processed file. If the variable doesn't exist, check the contents of the file, create the "whoD" variable, and then save it to the file.
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
    %Check whether the file already contains blk (behavior) and flu (imaging) variables. 'ignore' is there to avoid errors if whoD is [];
    
    %     Loop to run conv on any ephys files that are missing the blk variable, or if instructBlk is set to 1. First, check whether the processing
    %     function for the particular experiment definition file exists. If it does, process the block.
    if (~contains({'eph'},['ignore'; whoD]) || redoTag) && strcmpi(x.expType, 'ephys') && contains(dataType, {'all'; 'eph'})
        fprintf('Converting ephys recording data for %s %s idx = %d\n', x.expDate,x.subject,i);
        convEphysFile(x);
        whoD = unique(who('-file', x.processedData));
    end
    
    %Loop to run conv2PData on any 2P files that are missing the flu variable, or if instruct2P is set to 1
    if (~contains({'flu'},['ignore'; whoD]) || redoTag) && contains(lower(x.expType), {'twophoton';'widefield'}) && contains(dataType, {'all'; 'flu'})
        switch lower(x.expType)
            case 'twophoton'
                fprintf('Converting 2P file for %s %s idx = %d\n', x.expDate,x.subject,i);
                conv2PData(x);
            case 'widefield'; continue; %convWidefieldData(x);
        end
    end
    
    if (~contains({'fus'},['ignore'; whoD]) || redoTag) && strcmpi(x.expType(1:3), 'fus') && contains(dataType, {'all'; 'fus'})
        fprintf('Converting fusi recording data for %s %s idx = %d\n', x.expDate,x.subject,i);
        convfUSiFile(x);
        whoD = unique(who('-file', x.processedData));
    end
    
    %If not converted through other processing steps, convert block file.
    if  (~contains({'blk'}, ['ignore'; whoD]) || redoTag) && contains(dataType, {'all'; 'blk'})
        fprintf('Converting block file for %s %s idx = %d\n', x.expDate,x.subject,i);
        convBlockFile(x);
        whoD = unique(who('-file', x.processedData));
    end
end
    prc.syncfolder(prc.pathFinder('processedDirectory'), prc.pathFinder('serverProcessedDirectory'), 0);
end

%% Fucntion to convert block files. Does some basic operations, then passes x to the helper function for that experimental definition
function [x, whoD] = convBlockFile(x)
x.oldBlock = load(x.rawBlock); x.oldBlock = x.oldBlock.block;
if exist(x.galvoLog, 'file'); x.galvoLog = load(x.galvoLog); else, x.galvoLog = 0; end
x.oldParams = load(x.rawParams); x.oldParams = x.oldParams.parameters;
[x.oldBlock,  x.galvoLog] = prc.removeFirstTrialFromBlock(x.oldBlock, x.galvoLog);
timeOffset = x.oldBlock.experimentStartedTime-x.oldBlock.events.expStartTimes;

%%
if timeOffset > 100
    fieldList = fieldnames(x.oldBlock.inputs);
    fieldList = fieldList(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), fieldList));
    for i = 1:length(fieldList); x.oldBlock.inputs.(fieldList{i}) = x.oldBlock.inputs.(fieldList{i}) + timeOffset; end
    
    fieldList = fieldnames(x.oldBlock.outputs);
    fieldList = fieldList(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), fieldList));
    for i = 1:length(fieldList); x.oldBlock.outputs.(fieldList{i}) = x.oldBlock.outputs.(fieldList{i}) + timeOffset; end
    
    fieldList = fieldnames(x.oldBlock.events);
    for i = 2:2:length(fieldList); x.oldBlock.events.(fieldList{i}) = x.oldBlock.events.(fieldList{i}) + timeOffset; end
    %%
end
x.oldBlock.galvoLog = x.galvoLog;
[x.standardizedBlock, x.standardizedParams] = prc.standardBlkNames(x.oldBlock, x.oldParams);
x.validTrials = x.standardizedBlock.events.repeatNumValues(1:length(x.standardizedBlock.events.endTrialTimes))==1;

if contains(x.expType, {'ephys'; 'fusi'})
    x.timeline = load(x.rawTimeline); x.timeline = x.timeline.Timeline;
    [x.standardizedBlock, x.aligned] = prc.alignBlock2Timeline(x.standardizedBlock, x.timeline, x.expDef);
else, x.alignment = 'none';
end
%%

fields2copy = {'subject'; 'expDate'; 'expNum'; 'rigName'; 'expType'; 'expDef'};
for i = 1:length(fields2copy); x.newBlock.(fields2copy{i}) = x.(fields2copy{i}); end
for i = 1:length(fields2copy); x.standardizedParams.(fields2copy{i}) = x.(fields2copy{i}); end

repeatPoints = [strfind(diff([-1000,x.standardizedBlock.inputs.wheelValues])~=0, [0 0]) ...
    strfind(abs(diff([-1000,x.standardizedBlock.inputs.wheelValues(1:2:end)]))>1, [0 0])*2];
wheelValue = x.standardizedBlock.inputs.wheelValues(setdiff(1:end, repeatPoints))';
wheelTime = x.standardizedBlock.inputs.wheelTimes(setdiff(1:end, repeatPoints))';
x.newBlock.rawWheelTimeValue = single([wheelTime wheelValue-wheelValue(1)]);
x.standardizedParams.totalTrials = length(x.standardizedBlock.events.endTrialTimes);
x.standardizedParams.minutesOnRig = round((x.standardizedBlock.experimentEndedTime-x.standardizedBlock.experimentInitTime)/60);

blockFunction = str2func(['prc.expDef.' x.expDef]);
x = blockFunction(x);
blk = x.newBlock;
prm = x.newParams;
raw = x.newRaw;

if ~exist(fileparts(x.processedData), 'dir'); mkdir(fileparts(x.processedData)); end
if ~exist(fileparts(strrep(x.processedData,'dData','dDataLite')), 'dir'); mkdir(fileparts(strrep(x.processedData,'dData','dDataLite'))); end
if ~x.existProcessed(1)
    whoD = {'blk'; 'prm'; 'raw'}; save(x.processedData, 'blk', 'prm', 'whoD', 'raw');
    save(strrep(x.processedData,'dData','dDataLite'), 'blk', 'prm');
else
    whoD = unique([who('-file', x.processedData); 'blk'; 'prm'; 'raw']); save(x.processedData, 'blk', 'prm', 'whoD', 'raw', '-append');  %#ok
    save(strrep(x.processedData,'dData','dDataLite'), 'blk', 'prm');
end
end

%%
function x = convfUSiFile(x)
x = convBlockFile(x);
fus = struct;
fields2copy = {'subject'; 'expDate'; 'expNum'; 'rigName'; 'expType'; 'expDef'};
for i = 1:length(fields2copy); fus.(fields2copy{i}) = x.(fields2copy{i}); end
fields2copy = fields(x.aligned);
for i = 1:length(fields2copy); fus.(fields2copy{i}) = x.aligned.(fields2copy{i}); end
whoD = unique([who('-file', x.processedData); 'fus']);
save(x.processedData, 'fus', 'whoD', '-append');
copyfile(x.processedData, strrep(x.rawTimeline, 'Timeline', 'ProcBlock'));
end

%%
function convEphysFile(x)
siteList = dir([fileparts(x.kilosortOutput) '\*site*']);
if isempty(siteList); siteList = {x.kilosortOutput}; else; siteList = cellfun(@(y) [x.kilosortOutput '\' y], {siteList.name}', 'uni', 0); end
sites2Process = cellfun(@(x) ~exist([x '/spike_templates.npy'], 'file'), siteList);
if ~exist(x.kilosortOutput, 'dir') || any(sites2Process)
    kil.preProcessPhase3(x.subject, x.expDate, sites2Process);
    cellfun(@kil.liklihoodNoise, siteList);
elseif any(cellfun(@(x) ~exist([x '\cluster_pNoise.tsv'], 'file'), siteList))
    cellfun(@kil.liklihoodNoise, siteList);
end

clusterGroups = cell2mat(cellfun(@(x) tdfread([x '\cluster_group.tsv']),siteList,'uni', 0));
clusterpNoise = cell2mat(cellfun(@(x) tdfread([x '\cluster_pNoise.tsv']),siteList,'uni', 0));
%%
if length(vertcat(clusterpNoise.cluster_id))==length(vertcat(clusterGroups.cluster_id))
    fprintf('%s %s has been spike sorted. Loading and aligning data now... \n', x.expDate,x.subject);
    x = convBlockFile(x);
    eph = cell2mat(cellfun(@(y) kil.loadEphysData(x, y), siteList, 'uni', 0));
    numClusters = [0;cellfun(@length, {eph.clusterDepths}')];
    eph(1).spikeSite = uint8(cell2mat(arrayfun(@(x) eph(x).spikeTimes*0+x, 1:length(eph), 'uni', 0)'));
    eph(1).clusterSite = cell2mat(arrayfun(@(x) zeros(numClusters(x+1),1)+x, 1:length(eph), 'uni', 0)');
    eph(1).spikeCluster = cell2mat(arrayfun(@(x) eph(x).spikeCluster+sum(numClusters(1:x)), 1:length(eph), 'uni', 0)');
    combineList = {'spikeTimes';'spikeAmps';'clusterDepths';'clusterDuration';'clusterWaveforms';'clusterAmps'};
    for i = 1:length(combineList); eph(1).(combineList{i}) = cat(1,eph.(combineList{i})); end
    eph(1).clusterTemplates = [eph.clusterTemplates]';
    eph(1).channelMap = [eph.channelMap]';
    eph = eph(1);
    
    whoD = unique([who('-file', x.processedData); 'eph']); save(x.processedData, 'eph', 'whoD', '-append');
    blk = prc.combineBlocks(prc.getFilesFromDates(x.subject, x.expDate, 'bloprmeph', x.expDef));
    eph.clusterSigLevel = kil.findResponsiveCells(blk);
    save(x.processedData, 'eph', 'whoD', '-append');
else
    fprintf('%s %s must be spike sorted before processing further \n', x.expDate,x.subject);
end
if ~exist(x.processedData, 'file'); convBlockFile(x); end
end

%%
function conv2PData(x)
fprintf('Converting 2P data file for %s %s\n', x.expDate,x.subject);
if  ~exist(x.out2, 'dir') || isempty(dirP([x.out2 '*\\F*Plane*']))
    fprintf('Running Suite2P for %s %s\n', x.expDate,x.subject);
    sFOV = str2double({x.expList(strcmp(x.out2, {x.expList.out2}')).expNum});
    runSuite2P(x, sFOV);
end
fLst = dirP([x.out2 '*\\F*proc.mat']);
if exist(x.out2, 'dir') && isempty(fLst)
    fprintf('NOTE: Correct 2P data for %s %s Skipping...\n', x.expDate,x.subject);
    return;
end
t = load(prc.pathFinder('rawtime', x)); t = t.Timeline;
idx = find(strcmp(x.expNum, x.idxC));
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
flu.sNam = x.expNum;

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

