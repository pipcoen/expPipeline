function expDets = getExpDets(subject, expDate, expNum, folder)
%%
if ischar(subject); subject = {subject}; end
if ischar(expDate); expDate = {expDate}; end
if ischar(expNum); expNum = {expNum}; end
if ischar(folder); folder = {folder}; end

ephysRecord = load(prc.pathFinder('ephysrecord')); ephysRecord = ephysRecord.ephysRecord;
allSubjects = {ephysRecord.subject};
allDates = {ephysRecord.expDate};
allExpNums = {ephysRecord.expNum};
allFolders = {ephysRecord.folder};;
uniqueRecords = cellfun(@(w,x,y,z) [w,x,y,z], allSubjects(:), allDates(:), allExpNums(:), allFolders(:), 'uni', 0);
requestedRecords = cellfun(@(w,x,y,z) [w,x,y,z], subject, expDate, expNum, folder, 'uni', 0);
requestedRecordsAll = cellfun(@(w,x,y,z) [w,x,y,z], subject, expDate, repmat({'all'}, length(folder),1), folder, 'uni', 0);
[~, selectedIdx] = ismember(requestedRecords, uniqueRecords);
[~, selectedIdxAll] = ismember(requestedRecordsAll, uniqueRecords);
selectedIdx = max([selectedIdx selectedIdxAll], [], 2);
expDets = ephysRecord(selectedIdx(selectedIdx>0));
for i = 1:length(expDets); expDets(i).ephysRecordIdx = selectedIdx(i); end


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
