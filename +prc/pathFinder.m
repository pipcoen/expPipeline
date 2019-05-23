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
processedFolder = [driveName '\Dropbox (Neuropixels)\MouseData\ProcessedData\'];
sharedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
origBackup = [driveName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];
mouseRecord = [driveName '\Dropbox (Neuropixels)\MouseData\MouseRecord.xlsx'];

if contains(lower(pathType), {'directorycheck';'processedcheck';'explist'})
    if contains(hostname, {'homerig'; 'ziptop'}), directoryCheck = 'local';
    elseif strcmp(hostname, {'zip'}), directoryCheck = 'all';
    else; directoryCheck = 'server';
    end
end
%%
switch lower(pathType)
    case 'origblock'; pathOut = [expInfo subjectPath expRef '_Block.mat'];
    case 'origparams'; pathOut = [expInfo subjectPath expRef '_Parameters.mat'];
    case 'origblockfolder'; pathOut = [expInfo subjectPath(1:end-1)];
    case 'origtimeline'; pathOut = [expInfo subjectPath expRef '_Timeline.mat'];
    case 'backupblock'; pathOut = [origBackup subjectPath expRef '_Block.mat'];
    case 'backupparams'; pathOut = [origBackup subjectPath expRef '_parameters.mat'];
    case 'backupfolder'; pathOut = [origBackup subjectPath];
    case 'galvolog'; pathOut = [origBackup subjectPath expRef '_galvoLog.mat'];
    case 'processeddata'; pathOut = [processedFolder processedFileName];
    case 'shareddata'; pathOut = [sharedFolder processedFileName];
    case 'origprobedata'; pathOut = [expInfo subject '\' expDate '\ephys'];
    case 'kilosortoutput'; pathOut = [expInfo subject '\' expDate '\ephys\kilosort'];
    case 'explist'; pathOut = [processedFolder 'expList.mat'];
    case 'sharedfolder'; pathOut = sharedFolder;
    case 'expinfo'; pathOut = expInfo;
    case 'rawbackup'; pathOut = origBackup;
    case 'mouserecord'; pathOut = mouseRecord;
end
