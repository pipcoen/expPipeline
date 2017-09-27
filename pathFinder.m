function [pathOut] = pathFinder(pathType, subject, expDate, sessionNum) 
%% A funciton to return various paths used in processing and anlysis. Changes in file struture should be reflected here.
if nargin < 1; error('pathType required'); end
if nargin < 2; subject = '.'; end
if nargin < 3; expDate = datestr(today, 'yyyy-mm-dd'); end
if nargin < 4; sessionNum = '1'; end
if isnumeric(sessionNum); sessionNum = num2str(sessionNum); end

subjectPath = [subject '\' expDate '\' sessionNum '\'];
expRef = [expDate '_' sessionNum '_' subject];
processedFileName = [subject '\' subject '_' expDate([3:4 6:7 9:10]) '_' sessionNum  'Proc.mat'];

expInfo = '\\zserver.cortexlab.net\Data2\Subjects\';
processedFolder = 'E:\Dropbox (Neuropixels)\MouseData\ProcessedData\';
sharedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
suite2POutput = 'E:\Dropbox (Neuropixels)\subjectData\Suite2POutput\';
rawBackup = 'E:\Dropbox (Neuropixels)\MouseData\RawBehavior\';

if contains(lower(pathType), {'directorycheck';'processedcheck';'explist'})
    if strcmp(hostname, 'homerig'), directoryCheck = [1 0];
    else; directoryCheck = cellfun(@(x) exist(x, 'dir'), {processedFolder;sharedFolder})>0;
    end
end
%%
switch lower(pathType)
    case 'origblock'; pathOut = [expInfo subjectPath expRef '_Block.mat'];
    case 'origblockfolder'; pathOut = [expInfo subjectPath(1:end-1)];
    case 'rawblock'; pathOut = [rawBackup subjectPath expRef '_Block.mat'];
    case 'rawtimeline'; pathOut = [rawBackup subjectPath expRef '_Timeline.mat'];
    case 'galvolog'; pathOut = [rawBackup subjectPath expRef '_galvoLog.mat'];
    case 'rawparameters'; pathOut = [rawBackup subjectPath expRef '_parameters.mat'];
    case 'processeddata'; pathOut = [processedFolder processedFileName];
    case 'shareddata'; pathOut = [sharedFolder processedFileName];
    case 'suite2poutput'; pathOut = [suite2POutput subjectPath(1:end-1)];
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
    case 'mouserecord'; pathOut = 'E:\Dropbox (Neuropixels)\MouseData\MouseRecord.xlsx';
end
