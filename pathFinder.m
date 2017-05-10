function [pathOut] = pathFinder(pathType, mouseName, expDate, sessionNum) 
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('mouseNam', 'var'); mouseName = ''; end
if ~exist('expDate', 'var'); expDate = '0000-00-00'; end
if ~exist('sessionNum', 'var'); sessionNum = ''; end

    
addRequired(inp,'pathType',@ischar);
addOptional(inp,'mouseName', '',@ischar);
addOptional(inp,'expDate', '0000-00-00',@ischar);
addOptional(inp,'sessionNum',@ischar);
inp.parse(pathType, varargin{:})

mousePath = [mouseName '\' expDate '\' sessionNum '\'];
expRef = [expDate '_' sessionNum '_' mouseName];
processedFileName = [mouseName '\' mouseName '_' expDate([3:4 6:7 9:10]) '_' sessionNum  'Proc.mat'];

processedFolder = 'E:\Dropbox (Neuropixels)\MouseData\ProcessedData\';
sharedFolder = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
suite2POutput = 'E:\Dropbox (Neuropixels)\MouseData\Suite2POutput\';
rawBlock = ['\\zserver.cortexlab.net\Data\expInfo\' mousePath];

directoryCheck = cellfun(@(x) exist(x, 'directory'), {processedFolder;sharedFolder});

%%

switch lower(pathType)
    case 'block'; pathOut = [rawBlock expRef '_Block'];
    case 'timeline'; pathOut = [rawBlock expRef '_Timeline'];
    case 'paramaters'; pathOut = [rawBlock expRef '_parameters'];
    case 'processed'; pathOut = [processedFolder processedFileName];
    case 'shared'; pathOut = [sharedFolder processedFileName];
    case 'suite2poutput'; pathOut = [suite2POutput mousePath(1:end-1)];
    case 'explist'
        if a == 1
        elseif  exist([labL 'eLst.mat'], 'file'); pathOut = [labL 'eLst.mat'];
        else, error('No access to experimental list');
        end
    case 'bothelst'
        pathOut{1,1} = [savL 'eLst.mat'];
        pathOut{2,1} = [labL 'eLst.mat'];
    case 'savl' 
       pathOut = savL;
    case 'labl'
       pathOut = labL;
end
