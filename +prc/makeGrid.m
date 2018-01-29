function [gridData, gridXY] = makeGrid(blocks, data, operation, gridOpt)
%Make grid is a function for separating trials into a grid of audiovisual conditions. It operates on the output of concatenate blocks.

if ~exist('data', 'var'); error('Need data to sort into grid'); end
if ~exist('blocks', 'var'); error('Need block information to sort into grid'); end
if ~exist('operation', 'var'); operation = @sum; end
if ~exist('girdOpt', 'var') || ~isfield(gridOpt); gridOpt.type = 'condition'; end
if ~exist('girdOpt', 'var') || ~isfield(gridOpt); gridOpt.split = 0; end

sessionIdx = blocks.sessionIdx;
sessions = unique(sessionIdx);


conditions = blocks.conditionLabel;
gridIdx = blocks.grids.conditions;
switch lower(gridOpt.type)
    case 'abscondition'
        conditions = abs(conditions);
    case 'galvouni'
        a = 5;
end

fullGrid = repmat(gridIdx,[1,1,length(sessions)]);
repSessions = arrayfun(@(x) gridIdx*0+x, sessions, 'uni', 0);
repSessions = cat(3,repSessions{:});

if ~exist('selectedSessions', 'var'); selectedSessions = 1:length(sessions); end

switch sortFlag
    case 0; gridData = arrayfun(@(x,y) operation(data(conditions==x & sessionIdx==y)), fullGrid, repSessions, 'uni', 0);
    case 1; gridData = cell2mat(arrayfun(@(x) operation(data(conditions==x)), gridIdx, 'uni', 0));
    case 2; gridData = arrayfun(@(x,y) operation(data(abs(conditions)==x & sessionIdx==y)), fullGrid, repSessions, 'uni', 0);
    case 3; gridData = cell2mat(arrayfun(@(x) operation(data(abs(conditions)==x)), gridIdx, 'uni', 0));
end

if iscell(gridData) && max(cellfun(@length, gridData(:))) == 1
    gridData = cellfun(@squeeze, num2cell(cell2mat(gridData(:,:,selectedSessions)),3), 'uni', 0);
end

gridXY = 5;
end
