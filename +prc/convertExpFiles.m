function convertExpFiles(redoTag, dataType, selectedSubjects, selectedDates)
%% A funciton process experimental recodings into more concise, matching, local copies for future analysis.

%Inputs(default values)
%redoTag(0)---------------------Value of 1 will redo analyses even if they are already processed
%dataType('all')----------------A tag to process only a specific file type ("flu" for 2P, "fus" for fusi, "eph" for ephys and "all" for all of them)
%selectedSubjects({'PC'; 'DJ'})-Run for a specific mouse, or group of mice containing this string(s)
%selectedSubjects({'1'})-Run for a specific date, or group of dates containing this string(s)

%Outputs
%An output files is generated for each experiment in the form subject_yymmdd_expNumProc.mat These files will contain the following
%"blk" is a structure comprising a reduced, processed block file which contains all essential information.
%"raw" is a structure with some raw, rarely needed information that isn't usually loaded (takes too much time)
%"fus" is a structure with fusi-specific information
%"eph" is a structure with ephys-specific information
%"whoD"is a cell will the names of all variables in the file (loading this is quicker than asking matlab to check)

%% Set default values, load experimental list, check which processed files already exist, etc.

%Set the default values for the inputs and check that they are of the right format
if ~exist('redoTag', 'var') || isempty(redoTag); redoTag = 0; end
if ~exist('dataType', 'var') || isempty(dataType); dataType = 'all'; end
if ~exist('selectedSubjects', 'var'); selectedSubjects = {'PC'; 'DJ'}; end
if ~exist('selectedDates', 'var'); selectedDates = {'1'}; end
if strcmp(hostname, 'zip'); zipComp = 1; else, zipComp = 0; end

%If running on Pip's lab computer, sync the processed data in dropbox with his shared folder on the server.
if zipComp; fprintf('Running on Zip so will sync local folder and server... \n');
    prc.syncfolder(prc.pathFinder('processedDirectory'), prc.pathFinder('serverProcessedDirectory'), 0);
end

%Scan for new files and update the experiment list. Filter the list based on subject and date inputs.
expList = prc.scanForNewFiles(0,[0 1]);
expList = prc.updatePaths(expList);
expList = expList(cellfun(@(x) contains(x, selectedSubjects), {expList.subject}));
expList = expList(cellfun(@(x) contains(x, selectedDates), {expList.expDate}));

%Remove expDefs that don't have a processing helper function and notify the user that these will be skipped.
expDefs2Run = unique({expList.expDef}');
expDefs2Remove = expDefs2Run(cellfun(@(x) isempty(which(['prc.expDef.' x])), expDefs2Run));
if ~isempty(expDefs2Remove)
    fprintf('Warning: the following expDefs will be skipped: \n');
    fprintf('-%s \n', expDefs2Remove{:});
    expList = expList(~contains({expList.expDef}', [expDefs2Remove; 'Temporal']));
end


%Check whether the processed files already exist. If on Pip's computer, then also delete any files that are processed but have since been tagged as
%excluded, along with any fusi-processed blocks on the server. Also, check that any processed fusi files have also been backed up into their
%corresponding server folders
javaHandles = cellfun(@(x) java.io.File(x), {expList.processedData}', 'uni', 0);
processedSize = cellfun(@(x) x.length, javaHandles);
existProcessed = processedSize>100000;
notExcluded = [expList.excluded]'~=1;

fusiIdx = contains({expList.expType}', 'fusi');
fusiList = find(fusiIdx);
kruminBackup = cellfun(@(x) strrep(x, 'Timeline', 'ProcBlock'), {expList.rawTimeline}', 'uni', 0);
javaHandlesKruminServer = cellfun(@(x) java.io.File(x), kruminBackup(fusiList), 'uni', 0);
javaHandlesProcessed = cellfun(@(x) java.io.File(x), {expList(fusiList).processedData}', 'uni', 0);
copyIdx = fusiList(cellfun(@(x) x.lastModified, javaHandlesKruminServer)<cellfun(@(x) x.lastModified, javaHandlesProcessed));
arrayfun(@(x) copyfile(expList(x).processedData, kruminBackup{x}), copyIdx);

if zipComp
    deleteIdx = ([expList.excluded]'==1 & existProcessed) | (processedSize>0 & processedSize<100000);
    cellfun(@delete, {expList(deleteIdx).processedData});
    cellfun(@delete, {expList(deleteIdx).serverProcessedData});
    cellfun(@(x) delete(strrep(x, 'Timeline', 'ProcBlock')), {expList(deleteIdx & fusiIdx).rawTimeline});
    
    javaHandles = cellfun(@(x) java.io.File(x), {expList.rawBlock}', 'uni', 0);
    missBackup = find(cellfun(@(x) x.length, javaHandles)==0 & notExcluded);
    if ~isempty(missBackup)
        missFolders = cellfun(@fileparts, {expList(missBackup).rawBlock}', 'uni', 0);
        cellfun(@(x) mkdir(x), missFolders(cellfun(@(x) ~exist(x, 'dir'), missFolders)));
        arrayfun(@(x,y) copyfile(prc.pathFinder('serverBlock', expList(x)), y{1}), missBackup, {expList(missBackup).rawBlock}');
        arrayfun(@(x,y) copyfile(prc.pathFinder('serverParams', expList(x)), y{1}), missBackup, {expList(missBackup).rawParams}');
    end
end
existProcessed = existProcessed(notExcluded);
expList = expList(notExcluded);

%If redoTag is 1, then process all files of selected animals and dates, otherwise only process ones that don't exist, or are of the experiment type
%"fusi" or "ephys" (since these have extra processing stages).
if redoTag==1; files2Run = find(existProcessed*0+1)';
else, files2Run = find(~existProcessed | contains({expList.expType}', 'ephys'))';
end

%% Loop to process the files
%srtIdx can be more than zero if one wants to redo all files, but start in the middle. Useful if MATLAB crashes
srtIdx = 0;
for i = files2Run(files2Run>srtIdx)
    %create x, the processing structure. It is the expList entry for the current index, but also contains the entire list
    x = expList(i); x.expList = expList;
    x.existProcessed = existProcessed(i,:);
    
    %If the processed file exists, try to load the whoD variable (turn off warning if not found) which contains a list of the variable in that
    %processed file. If the variable doesn't exist, check the contents of the file, create the "whoD" variable, and then save it to the file later.
    warning('off', 'MATLAB:load:variableNotFound');
    if x.existProcessed(1)
        load(x.processedData, 'whoD');
        if ~exist('whoD', 'var')
            whoD = who('-file', x.processedData);
            save(x.processedData, 'whoD', '-append');
        end
    else
        if ~exist(fileparts(x.processedData), 'dir'); mkdir(fileparts(x.processedData)); end
        whoD = {'whoD'}; save(x.processedData, 'whoD');
    end
    warning('on', 'MATLAB:load:variableNotFound');
    %Check whether the file already contains blk (behavior) and flu (imaging) variables. 'ignore' is there to avoid errors if whoD is [];
    
    %Loop to convert ephys files
    if (~contains({'eph'},['ignore'; whoD]) || redoTag) && strcmpi(x.expType, 'ephys') && contains(dataType, {'all'; 'eph'})
        fprintf('Converting ephys recording data for %s %s idx = %d\n', x.expDate,x.subject,i);
        convEphysFile(x, redoTag);
        load(x.processedData, 'whoD');
    end
    
    %If not converted through other processing steps, convert block file.
    if  (~contains({'blk'}, ['ignore'; whoD]) || redoTag == 1) && contains(dataType, {'all'; 'blk'})
        fprintf('Converting block file for %s %s idx = %d\n', x.expDate,x.subject,i);
        convBlockFile(x);
    end
    clear whoD
end
prc.syncfolder(prc.pathFinder('processedDirectory'), prc.pathFinder('serverProcessedDirectory'), 0);
end

%% Fucntion to convert block files. Does some basic operations, then passes x to the helper function for that experimental definition
function convBlockFile(x)
%Load original block, params, and galvo log (if it exsits) Remove first trial because there are often timing issues with this trial
x.oldBlock = load(x.rawBlock); x.oldBlock = x.oldBlock.block;
x.oldParams = load(x.rawParams); x.oldParams = x.oldParams.parameters;
if contains(x.expType, 'inactivation'); x.galvoLog = load(x.galvoLog); else, x.galvoLog = 0; end
[x.oldBlock,  x.galvoLog] = prc.removeFirstTrialFromBlock(x.oldBlock, x.galvoLog);
timeOffset = x.oldBlock.experimentStartedTime-x.oldBlock.events.expStartTimes;

%Due to a silly bug, old experiments had huge offsets between signals event times and the stimwindowupdate times. 
if timeOffset > 100
    fieldList = fieldnames(x.oldBlock.inputs);
    fieldList = fieldList(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), fieldList));
    for i = 1:length(fieldList); x.oldBlock.inputs.(fieldList{i}) = x.oldBlock.inputs.(fieldList{i}) + timeOffset; end
    
    fieldList = fieldnames(x.oldBlock.outputs);
    fieldList = fieldList(cellfun(@(x) ~isempty(strfind(x, 'Times')>0), fieldList));
    for i = 1:length(fieldList); x.oldBlock.outputs.(fieldList{i}) = x.oldBlock.outputs.(fieldList{i}) + timeOffset; end
    
    %This is done differently just in case someone has given an event a name with "Times" in the title
    fieldList = fieldnames(x.oldBlock.events);
    for i = 2:2:length(fieldList); x.oldBlock.events.(fieldList{i}) = x.oldBlock.events.(fieldList{i}) + timeOffset; end
end

%This section standardizes the event names etc. for each experimetns. Events have changed since early instances of multiSpaceWorld
x.oldBlock.galvoLog = x.galvoLog;
[x.standardizedBlock, x.standardizedParams] = prc.standardBlkNames(x.oldBlock, x.oldParams);

%If experiment includes a meaningful timeline, then load timeline, align with the signals block, and also extract visual and auditory timings from the
%timeline, as well as reward and movement onset.
if contains(x.expType, {'ephys'; 'fusi'})
    x.timeline = load(x.rawTimeline); x.timeline = x.timeline.Timeline;
    [x.standardizedBlock, x.aligned] = prc.alignBlock2Timeline(x.standardizedBlock, x.timeline, x.expDef);
else, x.aligned = [];
end

%This copies standard fields over to newBlock and standardizedParams
fields2copy = {'subject'; 'expDate'; 'expNum'; 'rigName'; 'expType'; 'expDef'};
for i = 1:length(fields2copy); x.newBlock.(fields2copy{i}) = x.(fields2copy{i}); end
for i = 1:length(fields2copy); x.standardizedParams.(fields2copy{i}) = x.(fields2copy{i}); end

%Section remove repeats from the wheel inputs (this makes the biggest difference to the size of the saved block file)
repeatPoints = [strfind(diff([-1000,x.standardizedBlock.inputs.wheelValues])~=0, [0 0]) ...
    strfind(abs(diff([-1000,x.standardizedBlock.inputs.wheelValues(1:2:end)]))>1, [0 0])*2];
wheelValue = x.standardizedBlock.inputs.wheelValues(setdiff(1:end, repeatPoints))';
wheelTime = x.standardizedBlock.inputs.wheelTimes(setdiff(1:end, repeatPoints))';
x.newBlock.rawWheelTimeValue = single([wheelTime wheelValue-wheelValue(1)]);

%Run the block funciton for the expDef that was used for the mouse. This adds all relevant fields to "newBlock"
blockFunction = str2func(['prc.expDef.' x.expDef]);
x = blockFunction(x);

blk = x.newBlock;
raw = x.newRaw;
tim = x.aligned;

if ~isempty(tim); whoD = {'blk'; 'raw'; 'tim'; 'whoD'}; save(x.processedData, 'blk', 'whoD', 'raw', 'tim', '-append');
else, whoD = {'blk'; 'raw'; 'whoD'}; save(x.processedData, 'blk', 'whoD', 'raw', '-append');
end
end

%%
function convEphysFile(x, redoTag)
%If not converted through other processing steps, convert block file.
load(x.processedData, 'whoD');
if  ~contains({'blk'}, whoD) || redoTag == 1
    fprintf('Converting block file for %s %s \n', x.expDate,x.subject);
    convBlockFile(x);
    load(x.processedData, 'whoD');
end

siteList = dir([fileparts(x.kilosortOutput) '\*site*']);
if isempty(siteList); siteList = {x.kilosortOutput}; else; siteList = cellfun(@(y) [x.kilosortOutput '\' y], {siteList.name}', 'uni', 0); end
sites2Process = cellfun(@(x) ~exist([x '/spike_templates.npy'], 'file'), siteList);
if ~exist(x.kilosortOutput, 'dir') || any(sites2Process)
    kil.preProcessPhase3(x.subject, x.expDate, sites2Process);
    cellfun(@kil.liklihoodNoise, siteList);
elseif any(cellfun(@(x) ~exist([x '\cluster_pNoise.tsv'], 'file'), siteList))
    cellfun(@kil.liklihoodNoise, siteList);
end
sites2Process = cellfun(@(x) ~exist([x '/lfpPowerSpectra.mat'], 'file'), siteList);
if any(sites2Process); kil.loadAndSpectrogramLFP(x.subject, x.expDate, sites2Process); end

clusterGroups = cell2mat(cellfun(@(x) tdfread([x '\cluster_group.tsv']),siteList,'uni', 0));
clusterpNoise = cell2mat(cellfun(@(x) tdfread([x '\cluster_pNoise.tsv']),siteList,'uni', 0));
%%
spikeSorted = length(vertcat(clusterpNoise.cluster_id))==length(vertcat(clusterGroups.cluster_id));
loadData = isempty(whoD) || ~any(contains({'ephTmp'; 'eph'}, whoD)) || redoTag;
if spikeSorted && loadData; fprintf('%s %s has been spike sorted. Loading and aligning data now... \n', x.expDate,x.subject);
elseif ~spikeSorted && loadData; fprintf('WARNING: %s %s needs to be SORTED. Processing temp file now... \n', x.expDate,x.subject);
else, fprintf('WARNING: %s %s needs to be SORTED. \n', x.expDate,x.subject);
end

if loadData
    eph = cell2mat(cellfun(@(y) kil.loadEphysData(x, y, spikeSorted), siteList, 'uni', 0));
    if spikeSorted; whoD = unique([whoD; 'eph']); save(x.processedData, 'eph', 'whoD', '-append');
    else, ephTmp = eph; whoD = unique([whoD; 'ephTmp']); save(x.processedData, 'ephTmp', 'whoD', '-append');
    end
end
end
