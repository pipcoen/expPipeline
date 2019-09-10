function [gridData, gridXY] = makeGrid(blks, data, operation, type, split)
%Make grid is a function for separating trials into a grid of audiovisual conditionLabels. It operates on the output of concatenate blks.

if ~exist('data', 'var'); error('Need data to sort into grid'); end
if ~exist('blks', 'var'); error('Need block information to sort into grid'); end
if ~exist('operation', 'var'); operation = @sum; end
if ~exist('type', 'var') || isempty(type); type = 'condition'; end
if ~exist('split', 'var') || isempty(split); split = 0; end

experiments = blks.tot.experiments;
conditionLabels = blks.tri.stim.conditionLabel;

catBlk = prc.catStructs(blks);
if length(uniquecell(catBlk.exp.conditionParametersAV))~=1; error('Blocks must have same parameter set to make grid'); end
paramsAV = blks.exp.conditionParametersAV{1};
audValues = unique(paramsAV(:,1));
visValues = unique(paramsAV(:,2));

[grids.visValues, grids.audValues] = meshgrid(visValues,audValues);
[~, gridIdx] = ismember(paramsAV, [grids.audValues(:) grids.visValues(:)], 'rows');
grids.conditionLabels = nan*grids.visValues;
grids.conditionLabels(gridIdx) = blks.exp.conditionLabels{1};

%%
gridIdx = num2cell(grids.conditionLabels);
gridXY = {audValues; visValues};
switch lower(type)
    case 'abscondition'
        conditionLabels = abs(conditionLabels);
    case 'galvouni'
        conditionLabels = blks.tri.inactivation.galvoPosition;
        LMAxis = unique(blks.tri.inactivation.galvoPosition(:,1));
        APAxis = unique(blks.tri.inactivation.galvoPosition(:,2));
        [gridXY{1}, gridXY{2}] = meshgrid(LMAxis, APAxis);
        gridIdx = arrayfun(@(x,y) [x y], gridXY{1}, gridXY{2}, 'uni', 0);
end

if split == 1
    fullGrid = repmat(gridIdx,[1,1,length(experiments)]);
    repSessions = arrayfun(@(x) zeros(size(gridIdx))+x, experiments, 'uni', 0);
    repSessions = num2cell(cat(3,repSessions{:}));
end

switch split
    case 0; gridData = cell2mat(cellfun(@(x) operation(data(all(conditionLabels==x,2))), gridIdx, 'uni', 0));
    case 1; gridData = cell2mat(cellfun(@(x,y) operation(data(all(conditionLabels==x,2) & blks.tri.expRef==y)), fullGrid, repSessions, 'uni', 0));
    case 2; gridData = cellfun(@(x) data(all(conditionLabels==x,2),:), gridIdx, 'uni', 0);
    case 3; gridData = cellfun(@(x) prc.filtStruct(data, all(conditionLabels==x,2)), gridIdx, 'uni', 0);
end
end
