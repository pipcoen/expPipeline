function expList = updatePaths(expList)
if ~exist('expList', 'var')
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
end
%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
paths2Add = {'rawBlock'; 'rawParams'; 'serverProbeData';'serverTimeline'; 'serverFusi'; ...
    'galvoLog'; 'serverFolder'; 'processedData'; 'serverProcessedData'; 'kilosortOutput'};
fieldNames = {'rawBlock'; 'rawParams'; 'rawProbeData'; 'rawTimeline'; 'rawFusiData'; ...
    'galvoLog'; 'serverFolder'; 'processedData'; 'serverProcessedData'; 'kilosortOutput'};
for i = 1:length(expList)
    newPaths = prc.pathFinder(paths2Add, expList(i));
    for j = 1:length(fieldNames)
        expList(i).(fieldNames{j}) = newPaths{j};
    end
end