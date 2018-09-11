function [blk, prm, raw] = getFilesFromDates(subject, requestedDates, dataType)
%% Function to load proessed files from dates. Works with files from the convertExpFiles funciton.
% Inputs(defaults)
% subject(required)------String of subjects name
% requestedDates('all')--Requested dates to load. Can be a single date, an cell array of dates, a range, 'lastXXX' etc.
% dataType('all')--------Indicates types of data to load. A string containing one or more of 'blk', 'prm', and 'raw'

if ~exist('subject', 'var'); error('Must specify subject'); end
if ~exist('requestedDates', 'var') || isempty(requestedDates); requestedDates = {'last'}; end
if ~exist('dataType', 'var'); dataType = 'all'; end
if ~iscell(requestedDates); requestedDates = {requestedDates}; end
if iscell(requestedDates{1}); requestedDates = requestedDates{1}; end

expList = load(prc.pathFinder('expList')); expList = expList.expList;
existDirectories = prc.pathFinder('directoryCheck'); 
if all(existDirectories==[0,1])
    newPathList = {expList.sharedData}'; 
    [expList.processedData] = newPathList{:}; 
end

% get list of references and dates for subject
selectedFiles = expList(strcmp({expList.subject}', subject) & [expList.excluded]'==0);
excludedFiles = contains({selectedFiles.expDef}', {'Temporal'; 'stimSparseNoiseUncorrAsync'; 'stimKalatsky'});
selectedFiles(excludedFiles)  = [];
if isempty(expList); warning(['No processed files matching ' subject]); blk = {}; prm = {}; raw = {}; return; end
expDates = datenum(cell2mat({selectedFiles.expDate}'));

switch lower(requestedDates{1}(1:3))
    case 'las'
        if numel(requestedDates{1})>4
            lastDate = str2double(requestedDates{1}(5:end));
            selectedDates = expDates(end-min([lastDate length(expDates)])+1:end);
        else; selectedDates = expDates(end);
        end
    case 'fir'
        if numel(requestedDates{1})>5
            lastDate = str2double(requestedDates{1}(6:end));
            selectedDates = expDates(1:min([length(expDates), lastDate]));
        else; selectedDates = expDates(end);
        end
    case 'yes'; selectedDates = expDates(end-1);
    case 'all'; selectedDates = expDates;
    case 'rng'; selectedDates = expDates(expDates>=datenum(requestedDates{2}) & expDates<=datenum(requestedDates{3}));
    otherwise; selectedDates = sort(datenum(requestedDates));
end

% find paths to existing block files
blk = {}; prm = {}; raw = {};
selectedFiles = {selectedFiles(ismember(expDates, selectedDates)).processedData}';
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