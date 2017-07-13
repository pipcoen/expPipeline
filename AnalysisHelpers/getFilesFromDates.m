function [blk, prm, raw] = getFilesFromDates(subject, requestedDates, dataType)
if ~exist('subject', 'var'); error('Must specify subject'); end
if ~exist('requestedDates', 'var') || isempty(requestedDates); requestedDates = {'last'}; end
if ~exist('dataType', 'var'); dataType = 'all'; end
if ~iscell(requestedDates); requestedDates = {requestedDates}; end

expList = load(pathFinder('expList')); expList = expList.expList;
existDirectories = pathFinder('directoryCheck'); 
if all(existDirectories==[0,1])
    newPathList = {expList.sharedData}'; 
    [expList.processedData] = newPathList{:}; 
end

% get list of references and dates for subject
selectedFiles = expList(strcmp({expList.subject}', subject) & [expList.excluded]'==0);
if isempty(expList); warning(['No processed files matching ' subject]); blk = []; prm = []; raw = []; return; end
expDates = datenum(cell2mat({selectedFiles.expDate}'));

switch lower(requestedDates{1}(1:3))
    case 'las'
        if numel(requestedDates{1})>4
            lastDate = str2double(requestedDates{1}(5:end));
            selectedDates = expDates(end-min([lastDate length(expDates)])+1:end);
        else; selectedDates = expDates(end);
        end
    case 'yes'; selectedDates = expDates(end-1);
    case 'all'; selectedDates = expDates;
    case 'rng'; selectedDates = expDates(expDates>=datenum(requestedDates{2}) & expDates<=datenum(requestedDates{3}));
    otherwise; selectedDates = sort(datenum(requestedDates));
end

% find paths to existing block files
selectedFiles = {selectedFiles(ismember(expDates, selectedDates)).processedData}';
blk = []; prm = []; raw = []; 
if isempty(selectedFiles); warning(['No processed files matching ' subject ' for requested dates']); return; end
if contains(lower(dataType), {'blk'; 'blo'; 'all'})
    selectedBlocks = cellfun(@(x) load(x, 'blk'), selectedFiles, 'uni', 0);
    blk = [selectedBlocks{:}]'; blk = [blk(:).blk]';
end
if contains(lower(dataType), {'prm'; 'par'; 'all'})
    selectedParams = cellfun(@(x) load(x, 'prm'), selectedFiles, 'uni', 0);
    prm = [selectedParams{:}]'; prm = [prm(:).prm]';
end
if contains(lower(dataType), {'raw'; 'all'})
    selectedParams = cellfun(@(x) load(x, 'raw'), selectedFiles, 'uni', 0);
    raw = [selectedParams{:}]'; raw = [raw(:).raw]';
end
end