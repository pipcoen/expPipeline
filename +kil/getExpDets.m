function expList = getExpDets(expList)
%%
ephsRec = table2struct(readtable(prc.pathFinder('ephysrecord')));
tDat  = {ephsRec.expDate}'; tDat = cellfun(@(x) datestr(x, 'yyyy-mm-dd'), tDat, 'uni', 0);
[ephsRec.expDate] = deal(tDat{:});
[~,uniIdx,uniLabel] = unique(cell2table([{ephsRec.subject}' {ephsRec.expDate}', {ephsRec.expNum}']),'stable', 'rows');
for i = 1:length(uniIdx)
    expDets = ephsRec(uniLabel == i);
    expIdx = contains({expList.expDate}, {expDets.expDate}) & contains({expList.subject}, {expDets.subject});
    if contains('all', {expDets.expNum})
        [expList(expIdx).expDets] = deal(rmfield(expDets, {'subject', 'expDate', 'expNum', 'estDepth'}));
    end
end
