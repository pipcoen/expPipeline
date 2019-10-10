function [data, dates] = getDataFromDates(subject, requestedDates, expDef, extraData)
%% Function to load proessed files from dates. Works with files from the convertExpFiles funciton.
% Inputs(defaults)
% subject(required)------String of subjects name
% requestedDates('all')--Requested dates to load. Can be a single date, an cell array of dates, a range, 'lastXXX' etc.
% extraData('all')--------Indicates types of data to load. A string containing one or more of 'blk', 'prm', and 'raw'

if ~exist('subject', 'var'); error('Must specify subject'); end
if ~exist('requestedDates', 'var') || isempty(requestedDates); requestedDates = {'last'}; end
if ~exist('extraData', 'var'); extraData = 'all'; end
if ~exist('expDef', 'var'); expDef = 'multiSpaceWorld'; end
if ~iscell(requestedDates); requestedDates = {requestedDates}; end
if ~iscell(subject); subject = {subject}; end
load(prc.pathFinder('expList'), 'expList');

% get list of references and dates for subject
availableExps = expList(strcmp({expList.subject}', subject) & [expList.excluded]'~=1 & strcmp({expList.expDef}', expDef));
availableExps = prc.updatePaths(availableExps);
if isempty(availableExps); warning(['No processed files matching ' subject{1}]); return; end
availableDateNums = datenum(cell2mat({availableExps.expDate}'), 'yyyy-mm-dd');

selectedDateNums = cell(size(requestedDates,1),1);
for i = 1:size(requestedDates,1)
    currDat = requestedDates{i};
    if strcmpi(currDat(1:4), 'last')
        if numel(currDat)==4; currDat = [currDat '1']; end %#ok<*AGROW>
        lastDate = str2double(currDat(5:end));
        selectedDateNums{i} = availableDateNums(end-min([lastDate length(availableDateNums)])+1:end);
    elseif strcmpi(currDat(1:5), 'first')
        if numel(currDat)==5; currDat = [currDat '1']; end
        lastDate = str2double(currDat(6:end));
        selectedDateNums{i} = availableDateNums(1:min([length(availableDateNums), lastDate]));
    elseif strcmpi(currDat(1:4), 'yest');  selectedDateNums{i} = availableDateNums(end-1); 
    elseif strcmpi(currDat(1:3), 'all');  selectedDateNums{i} = availableDateNums;
    elseif contains(lower(currDat), ':')
        dateNums = datenum(strsplit(currDat, ':')', 'yyyy-mm-dd');
        selectedDateNums{i} = availableDateNums(availableDateNums>=dateNums(1) & availableDateNums<=dateNums(2));
    else, selectedDateNums = datenum(requestedDates, 'yyyy-mm-dd');
    end
end
if iscell(selectedDateNums); selectedDateNums = unique(cell2mat(selectedDateNums)); end

selectedFiles = {availableExps(ismember(availableDateNums, selectedDateNums)).processedData}';
if isempty(selectedFiles); warning(['No processed files matching ' subject{1} ' for requested dates']); return; end
whoD = cellfun(@(x) load(x, 'whoD'), selectedFiles, 'uni', 0);
whoD = [whoD{:}]'; whoD = {whoD(:).whoD}';
if any(~cellfun(@(x) contains('blk', x), whoD)); error('No blk for one or more requested files'); end
if any(~cellfun(@(x) contains('raw', x), whoD)); error('No raw for one or more requested files'); end

blk = cellfun(@(x) load(x, 'blk'), selectedFiles, 'uni', 0);
blk = [blk{:}]'; blk = [blk(:).blk]';

if contains(lower(extraData), {'raw'; 'all'})
    raw = cellfun(@(x) load(x, 'raw'), selectedFiles, 'uni', 0);
    raw = [raw{:}]'; raw = [raw(:).raw]';
    [blk.raw] = deal(raw);
end

if contains(lower(extraData), {'eph'; 'all'})
    timelineAvailable = find(cellfun(@(x) contains('tim', x), whoD));
    timeline = cellfun(@(x) load(x, 'tim'), selectedFiles(timelineAvailable), 'uni', 0);
    if ~isempty(timeline)
        timeline = [timeline{:}]'; timeline = [timeline(:).tim]';
        for i = 1:length(timelineAvailable); blk(timelineAvailable(i)).timeline = timeline(i); end
    end
    
    ephysAvailable = find(cellfun(@(x) contains('eph', x), whoD));
    ephys = cellfun(@(x) load(x, 'eph'), selectedFiles(ephysAvailable), 'uni', 0);
    if ~isempty(ephys)
        ephys = [ephys{:}]';
        for i = 1:length(ephysAvailable); blk(ephysAvailable(i)).ephys = ephys(i).eph; end
    end
end

data = blk;
dates = 1;
end