function expList = updatePaths(expList)
if ~exist('expList', 'var')
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
end
%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
for i = 1:length(expList)
    expList(i).rawBlock = prc.pathFinder('backupBlock', expList(i));
    expList(i).rawParams = prc.pathFinder('backupParams', expList(i));
    expList(i).rawBlock = prc.pathFinder('serverBlock', expList(i));
    expList(i).rawParams = prc.pathFinder('serverParams', expList(i));
    expList(i).rawProbeData = prc.pathFinder('serverProbeData', expList(i));
    expList(i).rawTimeline = prc.pathFinder('serverTimeline', expList(i));
    expList(i).galvoLog = prc.pathFinder('galvoLog', expList(i));
    expList(i).serverFolder = prc.pathFinder('serverFolder', expList(i));
    expList(i).processedData = prc.pathFinder('processedData', expList(i));
    expList(i).kilosortOutput = prc.pathFinder('kilosortOutput', expList(i));
end