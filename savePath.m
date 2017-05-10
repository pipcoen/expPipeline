function [sPth] = savePath(rTyp, mNam, rDat, sNum) 
if nargin == 1; mNam = ''; rDat = ''; sNum = ''; end
sep  = filesep;
mPth = [mNam sep rDat sep sNum sep];

savL = 'E:\Dropbox (Neuropixels)\MouseData\ProcessedData\';
labL = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
imoL = 'E:\Dropbox (Neuropixels)\MouseData\Suite2POutput\';
blkL = ['\\zserver.cortexlab.net\Data\expInfo\' mPth];
switch lower(rTyp)
    case 'rawblock'
        sPth = [blkL rDat '_' sNum '_' mNam '_Block'];
    case 'rawtime'
        sPth = [blkL rDat '_' sNum '_' mNam '_Timeline'];
    case 'params'
        sPth = [blkL rDat '_' sNum '_' mNam '_parameters'];
    case 'findata'
        sPth = [savL mNam sep mNam '_' rDat([3:4 6:7 9:10]) '_' sNum 'Proc.mat'];
    case 'labback'
        sPth = [labL mNam sep mNam '_' rDat([3:4 6:7 9:10]) '_' sNum 'Proc.mat'];
    case 'out2psuite'
        sPth = [imoL mPth(1:end-1)];
    case 'elst'
        if exist([savL 'eLst.mat'], 'file'); sPth = [savL 'eLst.mat'];
        elseif  exist([labL 'eLst.mat'], 'file'); sPth = [labL 'eLst.mat'];
        else, error('No access to experimental list');
        end
    case 'bothelst'
        sPth{1,1} = [savL 'eLst.mat'];
        sPth{2,1} = [labL 'eLst.mat'];
    case 'savl' 
       sPth = savL;
    case 'labl'
       sPth = labL;
end
