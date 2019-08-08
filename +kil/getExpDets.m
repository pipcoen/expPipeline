function expList = getExpDets(expList)
%%
ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
uniqueTags = cellfun(@(x,y,z) [x y x], {ephysRecord.subject}', {ephysRecord.expDate}', {ephysRecord.expNum}', 'uni', 0);
[~,uniIdx,uniLabel] = unique(uniqueTags, 'stable');
for i = 1:length(uniIdx)
    expDets = ephysRecord(uniLabel == i);
    expIdx = contains({expList.expDate}, {expDets.expDate}) & contains({expList.subject}, {expDets.subject});
    if contains('all', {expDets.expNum}) && any(expIdx)
        [expList(expIdx).expDets] = deal(rmfield(expDets, {'subject', 'expDate', 'expNum', 'estDepth'}));
    end
end
