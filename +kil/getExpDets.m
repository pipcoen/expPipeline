function expDets = getExpDets(subject, expDate, expNum, folder)
%%
ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
uniqueRecords = cellfun(@(w,x,y,z) [w,x,y,z], {ephysRecord.subject}', {ephysRecord.expDate}', {ephysRecord.expNum}', {ephysRecord.folder}', 'uni', 0);
requestedRecords = cellfun(@(w,x,y,z) [w,x,y,z], subject, expDate, expNum, folder, 'uni', 0);
requestedRecordsAll = cellfun(@(w,x,y,z) [w,x,y,z], subject, expDate, repmat({'all'}, length(folder),1), folder, 'uni', 0);
[~, selectedIdx] = ismember(requestedRecords, uniqueRecords);
[~, selectedIdxAll] = ismember(requestedRecordsAll, uniqueRecords);
selectedIdx = max([selectedIdx selectedIdxAll], [], 2);
expDets = ephysRecord(selectedIdx(selectedIdx>0));


% for i = 1:length(uniIdx)
%     expDets = ephysRecord(uniLabel == i);
%     expIdx = contains({expList.expDate}, {expDets.expDate}) & contains({expList.subject}, {expDets.subject});
%     if isstruct(subject)
%         expList = subject;
%         if contains('all', {expDets.expNum}) && any(expIdx)
%             [expList(expIdx).expDets] = deal(rmfield(expDets, {'subject', 'expDate', 'expNum', 'estDepth'}));
%         end
%     end
% end
