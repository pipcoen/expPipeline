function grid = makeGrid(blocks, data, operation, sortFlag, selectedSessions)
%Make grid is a function for separating trials into a grid of audiovisual conditions. It operates on the output of concatenate blocks.

if ~exist('data', 'var'); error('Need data to sort into grid'); end
if ~exist('blocks', 'var'); error('Need block information to sort into grid'); end
if ~exist('operation', 'var'); operation = @sum; end
if ~exist('sortFlag', 'var'); sortFlag = 0; end

sessionIdx = blocks.sessionIdx;
conditions = blocks.conditionLabel;
sessions = unique(sessionIdx);
gridIdx = blocks.grids.conditions;

fullGrid = repmat(gridIdx,[1,1,length(sessions)]);
repSessions = arrayfun(@(x) gridIdx*0+x, sessions, 'uni', 0);
repSessions = cat(3,repSessions{:});

if ~exist('selectedSessions', 'var'); selectedSessions = 1:length(sessions); end

switch sortFlag
    case 0; grid = arrayfun(@(x,y) operation(data(conditions==x & sessionIdx==y)), fullGrid, repSessions, 'uni', 0);
    case 1; grid = cell2mat(arrayfun(@(x) operation(data(conditions==x)), gridIdx, 'uni', 0));
    case 2; grid = arrayfun(@(x,y) operation(data(abs(conditions)==x & sessionIdx==y)), fullGrid, repSessions, 'uni', 0);
    case 3; grid = cell2mat(arrayfun(@(x) operation(data(abs(conditions)==x)), gridIdx, 'uni', 0));
end

if iscell(grid) && max(cellfun(@length, grid(:))) == 1
    grid = cellfun(@squeeze, num2cell(cell2mat(grid(:,:,selectedSessions)),3), 'uni', 0);
end

end
