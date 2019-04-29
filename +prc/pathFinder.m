function [pathOut] = pathFinder(pathType, subject, expDate, expNum) 
%% A funciton to return various paths used in processing and anlysis. Changes in file struture should be reflected here.
if nargin < 1; error('pathType required'); end
if nargin < 2; subject = '.'; end
if nargin < 3; expDate = datestr(now, 'yyyy-mm-dd'); end
if nargin < 4; expNum = '1'; end
if isnumeric(expNum); expNum = num2str(expNum); end

subjectPath = [subject '\' expDate '\' expNum '\'];
expRef = [expDate '_' expNum '_' subject];
processedFileName = [subject '\' subject '_' expDate([3:4 6:7 9:10]) '_' expNum  'Proc.mat'];

if contains(hostname, 'ziggurat')
    stubName = 'C:\Users\carandini';
elseif contains(hostname, 'ziptop')
    stubName = 'C:';
else, stubName = 'D:';
end

expInfo = '\\zubjects.cortexlab.net\Subjects\';
processedFolder = [stubName '\Dropbox (Neuropixels)\MouseData\ProcessedData\'];
sharedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
suite2POutput = [stubName '\Dropbox (Neuropixels)\MouseData\Suite2POutput\'];
rawBackup = [stubName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];
mouseRecord = [stubName '\Dropbox (Neuropixels)\MouseData\MouseRecord.xlsx'];

if contains(lower(pathType), {'directorycheck';'processedcheck';'explist'})
    if contains(hostname, {'homerig'; 'ziptop'}), directoryCheck = [1 0];
    else; directoryCheck = cellfun(@(x) exist(x, 'dir'), {processedFolder;sharedFolder})>0;
    end
end
%%
switch lower(pathType)
    case 'origblock'; pathOut = [expInfo subjectPath expRef '_Block.mat'];
    case 'origparams'; pathOut = [expInfo subjectPath expRef '_Parameters.mat'];
    case 'origblockfolder'; pathOut = [expInfo subjectPath(1:end-1)];
    case 'rawblock'; pathOut = [rawBackup subjectPath expRef '_Block.mat'];
    case 'rawtimeline'; pathOut = [expInfo subjectPath expRef '_Timeline.mat'];
    case 'galvolog'; pathOut = [rawBackup subjectPath expRef '_galvoLog.mat'];
    case 'rawparameters'; pathOut = [rawBackup subjectPath expRef '_parameters.mat'];
    case 'processeddata'; pathOut = [processedFolder processedFileName];
    case 'shareddata'; pathOut = [sharedFolder processedFileName];
    case 'suite2poutput'; pathOut = [suite2POutput subjectPath(1:end-1)];
    case 'rawprobedata'; pathOut = [expInfo subject '\' expDate '\ephys'];
    case 'kilosortoutput'; pathOut = [expInfo subject '\' expDate '\ephys\kilosort'];
    case 'explist'
        if directoryCheck(1); pathOut = [processedFolder 'expList.mat'];
        elseif  directoryCheck(2); pathOut = [sharedFolder 'expList.mat'];
        else, error('No access to experimental list')
        end
    case 'dropboxlist'; pathOut = [processedFolder 'expList.mat'];
    case 'sharedlist'; pathOut = [sharedFolder 'expList.mat'];
    case 'directorycheck'; pathOut = directoryCheck;
    case 'processedcheck'; pathOut = directoryCheck;
    case 'processedfolder'; pathOut = processedFolder;
    case 'rawbackupfolder'; pathOut = [rawBackup subjectPath];
    case 'sharedfolder'; pathOut = sharedFolder;
    case 'expinfo'; pathOut = expInfo;
    case 'rawbackup'; pathOut = rawBackup;
    case 'mouserecord'; pathOut = mouseRecord;
end
