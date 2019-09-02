function jitter(data, opt)
createOpt = ~exist('opt', 'var');
if createOpt || ~isfield(opt, 'centerOperation'); opt.centerOperation = @nanmean; end
if createOpt || ~isfield(opt, 'errorOperation'); opt.errorOperation = @nanstd; end
if createOpt || ~isfield(opt, 'jitterOffset'); opt.scatterOffset = 0.1; end
if createOpt || ~isfield(opt, 'jitterSpread'); opt.jitterSpread = 0.15; end
if createOpt || ~isfield(opt, 'xTickLocations'); opt.xTickLocations = 1:1:length(data); end
if createOpt || ~isfield(opt, 'xTickLabels'); opt.xTickLabels = []; end
if createOpt || ~isfield(opt, 'showNumbers'); opt.showNumbers = 0; end
if createOpt || ~isfield(opt, 'createAxis'); opt.createAxis = 0; end
if createOpt || ~isfield(opt, 'faceColors'); opt.faceColors = repmat('k', length(data),1); end
if createOpt || ~isfield(opt, 'edgeColors'); opt.edgeColors = opt.faceColors; end
if createOpt || ~isfield(opt, 'pairs2test'); opt.pairs2test = []; end
if createOpt || ~isfield(opt, 'significanceTest'); opt.significanceTest = @ttest; end
if createOpt || ~isfield(opt, 'linkedGroups'); opt.linkedGroups = []; end
if createOpt || ~isfield(opt, 'yLimits'); opt.yLimits = []; end
if createOpt || ~isfield(opt, 'pairLines'); opt.pairLines = 1; end

data = data(:);
data = cellfun(@double, data, 'uni', 0);
if opt.createAxis; figure; end; hold on
centerPoints = cellfun(opt.centerOperation, data);
errorBarLength = cellfun(@opt.errorOperation, data);


%define X axis positions
for i = 1:length(data)
    loopOptions.markerfacecolor = opt.faceColors(i,:);
    loopOptions.color = opt.edgeColors(i,:);
    loopOptions.markersize = 5;
    loopOptions.LineWidth = 2;
    
    lineXLocation = opt.scatterOffset + opt.xTickLocations(i);
    lineYLimits = [centerPoints(i)-errorBarLength(i) centerPoints(i)+errorBarLength(i)];
    line(lineXLocation*[1,1],lineYLimits,loopOptions)
    
    loopOptions.LineWidth = 1.5;
    plot(lineXLocation,centerPoints(i), 'o', loopOptions); 
    
    loopOptions.markersize = 4.5;
    jitterValues = (rand(length(data{i}),1)-0.5)*opt.jitterSpread + lineXLocation-2*opt.scatterOffset;
    plot(jitterValues,data{i}, 'o', loopOptions);
end
if ~isempty(opt.pairs2test)
    [~, significance] = cellfun(@(x) opt.significanceTest(data{x(1)}, data{x(2)}), opt.pairs2test);
    xPositionOfPairs = cellfun(@(x) opt.xTickLocations(x), opt.pairs2test, 'uni',0);
    plt.sigstar(xPositionOfPairs,significance);
end
if ~isempty(opt.yLimits); ylim(opt.yLimits); end
set(gca,'XTick', opt.xTickLocations, 'XTickLabel',opt.xTickLabels);

