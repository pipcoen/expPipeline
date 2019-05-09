function gridSplitByRows(dataGrid, visValues, audValues, plotOpts)
if size(dataGrid,1) ~= length(audValues); error('Must provide an auditory value for each row'); end
realAudValues = audValues(~isinf(audValues));
colorChoices = plt.selectRedBlueColors(realAudValues);
for j = 1:length(audValues)
    if isinf(audValues(j)); plotOpts.Color = 'k'; 
    else, plotOpts.Color = colorChoices(j,:);
    end
    plot(visValues, dataGrid(j,:), plotOpts);
    hold on;
end