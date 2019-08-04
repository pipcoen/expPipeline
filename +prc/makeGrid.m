function [gridData, gridXY] = makeGrid(blocks, data, operation, type, split)
%Make grid is a function for separating trials into a grid of audiovisual conditions. It operates on the output of concatenate blocks.

if ~exist('data', 'var'); error('Need data to sort into grid'); end
if ~exist('blocks', 'var'); error('Need block information to sort into grid'); end
if ~exist('operation', 'var'); operation = @sum; end
if ~exist('type', 'var') || isempty(type); type = 'condition'; end
if ~exist('split', 'var') || isempty(split); split = 0; end

experiments = unique(blocks.trialExperimentIdx);
conditions = blocks.conditionLabelRow(:,1);
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

if split == 1
    fullGrid = repmat(gridIdx,[1,1,length(experiments)]);
    repSessions = arrayfun(@(x) zeros(size(gridIdx))+x, experiments, 'uni', 0);
    repSessions = num2cell(cat(3,repSessions{:}));
end

switch split
    case 0; gridData = cell2mat(cellfun(@(x) operation(data(all(conditions==x,2))), gridIdx, 'uni', 0));
    case 1; gridData = cell2mat(cellfun(@(x,y) operation(data(all(conditions==x,2) & blocks.trialExperimentIdx==y)), fullGrid, repSessions, 'uni', 0));
    case 2; gridData = cellfun(@(x) data(all(conditions==x,2),:), gridIdx, 'uni', 0);
    case 3; gridData = cellfun(@(x) prc.filtStruct(data, all(conditions==x,2)), gridIdx, 'uni', 0);
end
end
