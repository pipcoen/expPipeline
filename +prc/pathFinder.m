function [pathOut, directoryCheck] = pathFinder(pathType, subject, expDate, expNum) 
%% A funciton to return various paths used in processing and anlysis. Changes in file struture should be reflected here.
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('subject', 'var'); subject = '.'; end
if isstruct(subject) && ~exist('expDate', 'var') && ~exist('expNum', 'var');
    expDate = subject.expDate; 
    expNum = subject.expNum; 
    subject = subject.subject; 
end
if ~exist('expDate', 'var'); expDate = datestr(now, 'yyyy-mm-dd'); end
if ~exist('expNum', 'var'); expNum = '1'; end
if isnumeric(expNum); expNum = num2str(expNum); end

subjectPath = [subject '\' expDate '\' expNum '\'];
expRef = [expDate '_' expNum '_' subject];
processedFileName = [subject '\' subject '_' expDate([3:4 6:7 9:10]) '_' expNum  'Proc.mat'];

if contains(hostname, 'ziptop')
    driveName = 'C:';
else, driveName = 'D:';
end

expInfo = '\\zubjects.cortexlab.net\Subjects\';
rawBackup = [driveName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];

if contains(lower(pathType), {'directorycheck';'processedcheck';'explist'})
    if contains(hostname, {'homerig'; 'ziptop'}), directoryCheck = 'local';
    elseif strcmp(hostname, {'zip'}), directoryCheck = 'all';
        processedFolder = [driveName '\Dropbox (Neuropixels)\MouseData\ProcessedData\'];
    else; directoryCheck = 'server';
        processedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
    end
end
%%
switch lower(pathType)
    case 'rawblock'; pathOut = [expInfo subjectPath expRef '_Block.mat'];
    case 'rawparams'; pathOut = [expInfo subjectPath expRef '_Parameters.mat'];
    case 'rawblockfolder'; pathOut = [expInfo subjectPath(1:end-1)];
    case 'rawtimeline'; pathOut = [expInfo subjectPath expRef '_Timeline.mat'];
    case 'rawprobedata'; pathOut = [expInfo subject '\' expDate '\ephys'];
    case 'backupblock'; pathOut = [rawBackup subjectPath expRef '_Block.mat'];
    case 'backupparams'; pathOut = [rawBackup subjectPath expRef '_parameters.mat'];
    case 'backupfolder'; pathOut = [rawBackup subjectPath];
    case 'galvolog'; pathOut = [rawBackup subjectPath expRef '_galvoLog.mat'];
    case 'processeddata'; pathOut = [processedFolder processedFileName];
    case 'shareddata'; pathOut = [sharedFolder processedFileName];
    case 'kilosortoutput'; pathOut = [expInfo subject '\' expDate '\ephys\kilosort'];
    case 'explist'; pathOut = [processedFolder 'expList.mat'];
    case 'expinfo'; pathOut = expInfo;
end
