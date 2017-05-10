function [savePath] = New_savePath(dataType, varargin)
optInput = {'' '0000-00-00' ''}; 
optInput(1:length(varargin)) = varargin;
[mouseNam, exptDate, sessnNum] = optInput{:};

mousePth = [mouseNam '\' exptDate '\' sessnNum '\'];
oldFName = [exptDate '_' sessnNum '_' mouseNam];
newFName = [mouseNam '\' mouseNam '_' exptDate([3:4 6:7 9:10]) '_' sessnNum];

savedLoc = 'E:\Dropbox (Neuropixels)\MouseData\ProcessedData\';
shareLoc = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
out2PLoc = 'E:\Dropbox (Neuropixels)\MouseData\Suite2POutput\';
blockLoc = ['\\zserver.cortexlab.net\Data\expInfo\' mousePth];
switch lower(dataType)
    case 'rawblock'; savePath = [blockLoc oldFName '_Block'];
    case 'timeline'; savePath = [blockLoc oldFName '_Timeline'];
    case 'prmsfile'; savePath = [blockLoc oldFName '_parameters'];
    case 'procdata'; savePath = [savedLoc newFName 'Proc.mat'];
    case 'backdata'; savePath = [shareLoc newFName 'Proc.mat'];
    case 'out2ploc'; savePath = [out2PLoc mousePth(1:end-1)];
    case 'exptlist'; savePath = iff(exist([savedLoc 'eLst.mat'], 'file'), [savedLoc 'eLst.mat'], [shareLoc 'eLst.mat']);
    case 'alllists'; savePath{1,1} = [savedLoc 'eLst.mat']; savePath{2,1} = [shareLoc 'eLst.mat'];
    case 'savedloc'; savePath = savedLoc;
    case 'shareloc'; savePath = shareLoc;
end
