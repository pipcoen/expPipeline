function [gridData, gridXY] = makeGrid(blocks, data, operation, type, split)
%Make grid is a function for separating trials into a grid of audiovisual conditions. It operates on the output of concatenate blocks.

if ~exist('data', 'var'); error('Need data to sort into grid'); end
if ~exist('blocks', 'var'); error('Need block information to sort into grid'); end
if ~exist('operation', 'var'); operation = @sum; end
if ~exist('type', 'var'); type = 'condition'; end
if ~exist('split', 'var'); split = 0; end

sessions = unique(blocks.sessionIdx);
conditions = blocks.conditionLabel;
gridIdx = num2cell(blocks.grids.conditions);
gridXY = {blocks.audValues; blocks.visValues};
switch lower(type)
    case 'abscondition'
        conditions = abs(conditions);
    case 'galvouni'
        conditions = blocks.galvoPosition;
        LMAxis = unique(blocks.galvoPosition(:,1));
        APAxis = unique(blocks.galvoPosition(:,2));
        [gridXY{1}, gridXY{2}] = meshgrid(LMAxis, APAxis);
        gridIdx = arrayfun(@(x,y) [x y], gridXY{1}, gridXY{2}, 'uni', 0);
end

% fullGrid = repmat(gridIdx,[1,1,length(sessions)]);
% repSessions = arrayfun(@(x) gridIdx*0+x, sessions, 'uni', 0);
% repSessions = cat(3,repSessions{:});

switch split
    case 0; gridData = cell2mat(cellfun(@(x) operation(data(all(conditions==x,2))), gridIdx, 'uni', 0));
    case 1; gridData = arrayfun(@(x,y) operation(data(conditions==x & blocks.sessionIdx==y)), fullGrid, repSessions, 'uni', 0);
end
end
