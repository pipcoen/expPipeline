function gridSplitByRows(dataGrid, visValues, audValues, plotOpts)
if size(dataGrid,1) ~= length(audValues); error('Must provide an auditory value for each row'); end
colorChoices = plt.selectRedBlueColors(audValues);

for j = 1:length(audValues)
    plotOpts.Color = colorChoices(j,:);
    plot(visValues, dataGrid(j,:), plotOpts);
    hold on;
end