function expList = updatePaths(expList)
if ~exist('expList', 'var')
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
end
%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
paths2Add = {'rawBlock'; 'rawParams'; 'serverProbeData';'serverTimeline'; 'serverFusi'; ...
    'galvoLog'; 'serverFolder'; 'processedData'; 'serverProcessedData'; 'kilosortOutput'};
fieldNames = {'rawBlock'; 'rawParams'; 'rawProbeData'; 'rawTimeline'; 'rawFusiData'; ...
    'galvoLog'; 'serverFolder'; 'processedData'; 'serverProcessedData'; 'kilosortOutput'};
dateNums = num2cell(deal(dtstr2dtnummx({expList.expDate}', 'yyyy-MM-dd')));
newPaths = prc.pathFinder(paths2Add, {expList.subject}', {expList.expDate}', {expList.expNum}', dateNums);
for i = 1:length(expList)
    for j = 1:length(fieldNames)
        expList(i).(fieldNames{j}) = newPaths{i,j};
    end
end
