function [pathOut] = pathFinder(pathType, subject, expDate, sessionNum) 
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('subject', 'var'); subject = ''; end
if ~exist('expDate', 'var'); expDate = '0000-00-00'; end
if ~exist('sessionNum', 'var'); sessionNum = ''; end

subjectPath = [subject '\' expDate '\' sessionNum '\'];
expRef = [expDate '_' sessionNum '_' subject];
processedFileName = [subject '\' subject '_' expDate([3:4 6:7 9:10]) '_' sessionNum  'Proc.mat'];

expInfo = '\\zserver.cortexlab.net\Data\expInfo\';
processedFolder = 'E:\Dropbox (Neuropixels)\MouseData\ProcessedData\';
sharedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
suite2POutput = 'E:\Dropbox (Neuropixels)\subjectData\Suite2POutput\';
rawBlock = [expInfo subjectPath];

directoryCheck = cellfun(@(x) exist(x, 'dir'), {processedFolder;sharedFolder})>0;
%%
switch lower(pathType)
    case 'rawblock'; pathOut = [expInfo subjectPath expRef '_Block'];
    case 'rawtimeline'; pathOut = [rawBlock expRef '_Timeline'];
    case 'rawparameters'; pathOut = [rawBlock expRef '_parameters'];
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
    case 'sharedfolder'; pathOut = sharedFolder;
    case 'expinfo'; pathOut = expInfo;
end
