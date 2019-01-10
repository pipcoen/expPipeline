function expList = updateProcessedContents(expList, saveList)
if ~exist('expList', 'var')
    expList = load(prc.pathFinder('expList'), 'expList'); expList = expList.expList;
end
if ~exist('saveList', 'var'); saveList = 1; end
%% Collect all paths--this also ensures paths are up to date with any changes in prc.pathFinder
warning('off', 'MATLAB:load:variableNotFound');
for i = 1:length(expList)
    if expList.existProcessed(1)
        load(x.processedData, 'whoD');
        if ~exist('whoD', 'var')
            whoD = who('-file', x.processedData);
            save(x.processedData, 'whoD', '-append');
        end
    else; whoD = [];
    end
end
warning('on', 'MATLAB:load:variableNotFound');
if saveList
    save(prc.pathFinder('dropboxlist'), 'expList'); save(strrep(prc.pathFinder('dropboxlist'),'dData','dDataLite'), 'expList');
end