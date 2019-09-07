function [pathOut, directoryCheck] = pathFinder(pathType, subject, expDate, expNum, dateNumber)
%% A funciton to return various paths used in processing and anlysis. Changes in file struture should be reflected here.
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('subject', 'var'); [subject, expDate, expNum, dateNumber] = deal('noDataGiven'); end
if isstruct(subject) && ~exist('expDate', 'var') && ~exist('expNum', 'var')
    if isfield(subject, 'dateNum'); dateNumber = subject.dateNum; end
    expDate = subject.expDate;
    expNum = subject.expNum;
    subject = subject.subject;
end
if exist('expDate', 'var') && ~exist('dateNumber', 'var'); dateNumber = datenum(expDate, 'yyyy-mm-dd'); end
if ~exist('expDate', 'var'); expDate = 'noDataGiven'; end
if ~exist('expNum', 'var'); expNum = 'noDataGiven'; end

if ~iscell(pathType); pathType = {pathType}; end
if ~iscell(subject); subject = {subject}; end
if ~iscell(expDate); expDate = {expDate}; end
if ~iscell(expNum); expNum = {expNum}; end
if ~iscell(dateNumber); dateNumber = {dateNumber}; end

hostName = hostname;
if isnumeric(expNum); expNum = num2str(expNum); end
if isnumeric(expDate); expDate =  datestr(expDate, 'yyyy-mm-dd'); end

if contains(hostName, 'ziptop')
    driveName = 'C:';
else, driveName = 'D:';
end

rawBackup = [driveName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];
serverProcessedDirectory = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
processedDirectory = [driveName '\Dropbox (Neuropixels)\MouseData\ProcessedData\'];
if contains('rawBlock', pathType); pathType{contains(pathType, 'rawBlock')} = 'backupBlock'; end
if contains('rawParams', pathType); pathType{contains(pathType, 'rawParams')} = 'backupParams'; end

if contains(hostName, {'homerig'; 'ziptop'})
    directoryCheck = 'local';
elseif strcmp(hostName, {'zip'})
    directoryCheck = 'all';
else; directoryCheck = 'server';
    processedDirectory = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
    if contains('rawBlock', pathType); pathType{contains(pathType, 'rawBlock')} = 'serverBlock'; end
    if contains('rawParams', pathType); pathType{contains(pathType, 'rawParams')} = 'serverParams'; end
end

pathOut = cell(size(subject,1), length(pathType));
for i = 1:size(subject,1)
    subjectPath = [subject{i} '\' expDate{i} '\' expNum{i} '\'];
    expRef = [expDate{i} '_' expNum{i} '_' subject{i}];
    processedFileName = [subject{i} '\' subject{i} '_' expDate{i}([3:4 6:7 9:10]) '_' expNum{i}  'Proc.mat'];
    
    expInfo = {'\\zubjects.cortexlab.net\Subjects\'; '\\zserver.cortexlab.net\Data\Subjects\'};
    if ~strcmp(expDate{i}, 'noDataGiven')
        if dateNumber{i} > 737589 && strcmp(subject{i}, 'PC037'); expInfo = expInfo{2}; %'2019-06-13'
        elseif dateNumber{i} > 737590 && strcmp(subject{i}, 'PC038'); expInfo = expInfo{2}; %'2019-06-14'
        elseif dateNumber{i} < 737612; expInfo = expInfo{1}; %'2019-07-06'
        else, expInfo = expInfo{2};
        end
    end
    %%
    for j = 1:length(pathType)
        switch lower(pathType{j})
            case 'serverblock'; pathOut{i,j} = [expInfo subjectPath expRef '_Block.mat'];
            case 'serverparams'; pathOut{i,j} = [expInfo subjectPath expRef '_Parameters.mat'];
            case 'serverfolder'; pathOut{i,j} = [expInfo subjectPath];
            case 'servertimeline'; pathOut{i,j} = [expInfo subjectPath expRef '_Timeline.mat'];
            case 'serverprobedata'; pathOut{i,j} = [expInfo subject{i} '\' expDate{i} '\ephys'];
            case 'serverprocesseddata'; pathOut{i,j} = [serverProcessedDirectory processedFileName];
            case 'serverprocessedfolder'; pathOut{i,j} = [serverProcessedDirectory subject{i}];
            case 'serverprocesseddirectory'; pathOut{i,j} = serverProcessedDirectory;
            case 'serverfusi'; pathOut{i,j} = [expInfo subjectPath expRef '_fus.mat'];
            case 'backupblock'; pathOut{i,j} = [rawBackup subjectPath expRef '_Block.mat'];
            case 'backupparams'; pathOut{i,j} = [rawBackup subjectPath expRef '_parameters.mat'];
            case 'backupfolder'; pathOut{i,j} = [rawBackup subjectPath];
            case 'backupdirectory'; pathOut{i,j} = rawBackup;
            case 'processeddata'; pathOut{i,j} = [processedDirectory processedFileName];
            case 'processedfolder'; pathOut{i,j} = [processedDirectory subject{i}];
            case 'processeddirectory'; pathOut{i,j} = processedDirectory;
            case 'galvolog'; pathOut{i,j} = [expInfo subjectPath expRef '_galvoLog.mat'];
            case 'kilosortoutput'; pathOut{i,j} = [expInfo subject{i} '\' expDate{i} '\ephys\kilosort'];
            case 'explist'; pathOut{i,j} = [processedDirectory 'expList.mat'];
            case 'expinfo'; pathOut{i,j} = expInfo;
            case 'ephysrecord'; pathOut{i,j} = [processedDirectory 'ePhysRecord.mat'];
            case 'allenatlas'; pathOut{i,j} = [driveName '\Dropbox (Neuropixels)\MouseData\Atlas\allenCCF\'];
            case 'probepathdata'; pathOut{i,j} = [processedDirectory 'XHistology\' subject{i} '\probe_points.mat'];
        end
    end
end
if length(pathOut) == 1; pathOut = pathOut{1}; end
