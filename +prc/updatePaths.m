function expList = updatePaths(expList, saveList)
if ~exist('expList', 'var')
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
end
if ~exist('saveList', 'var'); saveList = 1; end
%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
for i = 1:length(expList)
    expList(i).rawBlock = prc.pathFinder('rawBlock', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).rawTimeline = prc.pathFinder('rawTimeline', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).rawFolder = prc.pathFinder('origblockfolder', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).rawParams = prc.pathFinder('rawParameters', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).processedData = prc.pathFinder('processedData', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).sharedData = prc.pathFinder('sharedData', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).suite2POutput = prc.pathFinder('suite2POutput', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).rawProbeData = prc.pathFinder('rawProbeData', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    expList(i).kilosortOutput = prc.pathFinder('kilosortOutput', expList(i).subject, expList(i).expDate, expList(i).sessionNum);
    if expList(i).galvoLog ~= 0; expList(i).galvoLog = prc.pathFinder('galvoLog', expList(i).subject, expList(i).expDate, expList(i).sessionNum); end
end
if saveList
    save(prc.pathFinder('dropboxlist'), 'expList'); save(strrep(prc.pathFinder('dropboxlist'),'dData','dDataLite'), 'expList');
end