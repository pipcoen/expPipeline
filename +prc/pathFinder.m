function [pathOut, directoryCheck] = pathFinder(pathType, subject, expDate, expNum) 
%% A funciton to return various paths used in processing and anlysis. Changes in file struture should be reflected here.
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('subject', 'var'); [subject, expDate, expNum] = deal('noDataGiven'); end
if isstruct(subject) && ~exist('expDate', 'var') && ~exist('expNum', 'var')
    expDate = subject.expDate; 
    expNum = subject.expNum; 
    subject = subject.subject; 
end
hostName = hostname;
if isnumeric(expNum); expNum = num2str(expNum); end

subjectPath = [subject '\' expDate '\' expNum '\'];
expRef = [expDate '_' expNum '_' subject];
processedFileName = [subject '\' subject '_' expDate([3:4 6:7 9:10]) '_' expNum  'Proc.mat'];

if contains(hostName, 'ziptop')
    driveName = 'C:';
else, driveName = 'D:';
end

expInfo = '\\zubjects.cortexlab.net\Subjects\';
rawBackup = [driveName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];

if contains(hostName, {'homerig'; 'ziptop'})
    directoryCheck = 'local';
    if strcmpi('rawBlock', pathType); pathType = 'backupBlock'; end
elseif strcmp(hostName, {'zip'})
    directoryCheck = 'all';
    processedFolder = [driveName '\Dropbox (Neuropixels)\MouseData\ProcessedData\'];
    if strcmpi('rawBlock', pathType); pathType = 'backupBlock'; end
else; directoryCheck = 'server';
    processedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
    if strcmpi('rawBlock', pathType); pathType = 'serverBlock'; end
end
serverFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
%%
switch lower(pathType)
    case 'serverblock'; pathOut = [expInfo subjectPath expRef '_Block.mat'];
    case 'serverparams'; pathOut = [expInfo subjectPath expRef '_Parameters.mat'];
    case 'serverfolder'; pathOut = [expInfo subjectPath(1:end-1)];
    case 'servertimeline'; pathOut = [expInfo subjectPath expRef '_Timeline.mat'];
    case 'serverprobedata'; pathOut = [expInfo subject '\' expDate '\ephys'];
    case 'backupblock'; pathOut = [rawBackup subjectPath expRef '_Block.mat'];
    case 'backupparams'; pathOut = [rawBackup subjectPath expRef '_parameters.mat'];
    case 'backupfolder'; pathOut = [rawBackup subjectPath];
    case 'galvolog'; pathOut = [rawBackup subjectPath expRef '_galvoLog.mat'];
    case 'processeddata'; pathOut = [processedFolder processedFileName];
    case 'serverdata'; pathOut = [processedFolder processedFileName];
    case 'kilosortoutput'; pathOut = [expInfo subject '\' expDate '\ephys\kilosort'];
    case 'explist'; pathOut = [processedFolder 'expList.mat'];
    case 'expinfo'; pathOut = expInfo;
end
