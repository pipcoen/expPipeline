function fracitonRightChoiceData(grids, plotOpt)
%% A plotting function that plots the fraction of right choices made by the mouse for each audiovisual condition in a blk
% INPUTS(default values)
% blk(required)-------------This is the block of behavioral data (after being loaded with spatial analysis)
% plotOpt-------------------Stucture with various plotting options
%	.errorBars(0)------------------Whether to include error bars
%	.axesType('normal')------------Can plot fraction of right choices in "normal" or "log" scale (log(pR/pL)
%	.lineType(.)-------------------Can be changed to connect dots with a line
%	.splitType('aud')--------------Whether to split colors by aud or vis values

if ~exist('plotOpt', 'var'); plotOpt = struct; end
if ~isfield(plotOpt, 'errorBars'); plotOpt.errorBars = 0; end
if ~isfield(plotOpt, 'errorType'); plotOpt.errorType = 'patch'; end
if ~isfield(plotOpt, 'plotType'); plotOpt.plotType = 'normal'; end
if ~isfield(plotOpt, 'lineStyle'); plotOpt.lineStyle = 'none'; end
if ~isfield(plotOpt, 'lineWidth'); plotOpt.lineWidth = 2; end
if ~isfield(plotOpt, 'Marker'); plotOpt.marker = '.'; end
if ~isfield(plotOpt, 'MarkerSize'); plotOpt.markerSize = 20; end
if ~isfield(plotOpt, 'contrastPower'); plotOpt.contrastPower = 0.7; end
if ~isfield(grids, 'numTrials'); grids.numTrials = ~isnan(grids.fracRightTurns); end

uniAud = unique(grids.audValues(:));
selectedColors = plt.selectRedBlueColors(uniAud);

if strcmp(plotOpt.plotType, 'log')
    grids.fracRightTurns = log((grids.fracRightTurns./(1-grids.fracRightTurns)));
    grids.visValues = (abs(grids.visValues).^plotOpt.contrastPower).*sign(grids.visValues);
    maxContrast = max(abs(grids.visValues(:)));
    xlim([-maxContrast maxContrast]);
end

for i = uniAud'
    idx = find(grids.audValues==i & grids.numTrials>0);
    
    lineOpt.lineStyle = plotOpt.lineStyle;
    lineOpt.MarkerSize = plotOpt.markerSize;
    lineOpt.Color = selectedColors(uniAud==i,:);
    lineOpt.Marker = plotOpt.marker;
    lineOpt.lineWidth = plotOpt.lineWidth;
    plot(grids.visValues(idx),grids.fracRightTurns(idx),lineOpt);

    if plotOpt.errorBars && strcmpi(plotOpt.errorType, 'patch')
        patchY = [grids.fracRightTurnsHighBound(idx); grids.fracRightTurnsLowBound(flipud(idx(:))); grids.fracRightTurnsHighBound(idx(1))];
        patchX = [grids.visValues(idx); grids.visValues(flipud(idx(:))); grids.visValues(idx(1))];
        patch(patchX, patchY, lineOpt.Color, 'FaceAlpha', 1, 'EdgeColor', 'none')
    elseif plotOpt.errorBars && strcmpi(plotOpt.errorType, 'line')
        lineOpt.Marker = 'none';
        lineOpt.lineStyle = '-';
        lineOpt.lineWidth = 2;
        arrayfun(@(x) plot(grids.visValues(x)*[1 1], [grids.fracRightTurnsLowBound(x) grids.fracRightTurnsHighBound(x)], lineOpt), idx)        
    end
    
    hold on;
end