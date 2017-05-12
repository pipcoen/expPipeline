function [pathOut] = pathFinder(pathType, subject, expDate, sessionNum) 
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('mouseNam', 'var'); subject = ''; end
if ~exist('expDate', 'var'); expDate = '0000-00-00'; end
if ~exist('sessionNum', 'var'); sessionNum = ''; end

mousePath = [subject '\' expDate '\' sessionNum '\'];
expRef = [expDate '_' sessionNum '_' subject];
processedFileName = [subject '\' subject '_' expDate([3:4 6:7 9:10]) '_' sessionNum  'Proc.mat'];

processedFolder = 'E:\Dropbox (Neuropixels)\MouseData\ProcessedData\';
sharedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
suite2POutput = 'E:\Dropbox (Neuropixels)\MouseData\Suite2POutput\';
rawBlock = ['\\zserver.cortexlab.net\Data\expInfo\' mousePath];

directoryCheck = cellfun(@(x) exist(x, 'dir'), {processedFolder;sharedFolder})>0;

%%

switch lower(pathType)
    case 'block'; pathOut = [rawBlock expRef '_Block'];
    case 'timeline'; pathOut = [rawBlock expRef '_Timeline'];
    case 'paramaters'; pathOut = [rawBlock expRef '_parameters'];
    case 'processed'; pathOut = [processedFolder processedFileName];
    case 'shared'; pathOut = [sharedFolder processedFileName];
    case 'suite2poutput'; pathOut = [suite2POutput mousePath(1:end-1)];
    case 'explist'
        if directoryCheck(1); pathOut = [processedFolder 'expList.mat'];
        elseif  directoryCheck(2); pathOut = [sharedFolder 'expList.mat'];
        else, error('No access to experimental list')
        end
    case 'directorycheck'; pathOut = directoryCheck;
    case 'processedfolder'; pathOut = processedFolder;
    case 'sharedfolder'; pathOut = sharedFolder;
end
