function [fitData] = gridSplitByRows(dataGrid, visValues, audValues, plotOpts)
if size(dataGrid,1) ~= length(audValues); error('Must provide an auditory value for each row'); end
realAudValues = audValues(~isinf(audValues));
colorChoices = plt.selectRedBlueColors(realAudValues);
for j = 1:length(audValues)
    if isinf(audValues(j)); plotOpts.Color = 'k'; 
    else, plotOpts.Color = colorChoices(j,:);
    end
    plot(visValues, dataGrid(j,:), plotOpts);
    hold on;
%     
%     if strcmp(plotOpts.lineStyle, 'none')
%         xDat = visValues(~isnan(dataGrid(j,:)) & ~isinf(dataGrid(j,:)));
%         yDat = dataGrid(j,(~isnan(dataGrid(j,:)) & ~isinf(dataGrid(j,:))));
%         lineFit = polyfit(xDat(:)', yDat(:)',1);
%         plot(visValues,polyval(lineFit, visValues), '-', 'linewidth', 2.5, 'color', colorChoices(j,:));
%         r2 = corrcoef(xDat(:)', yDat(:)').^2;
%         fitData.lineFit(j,:) = lineFit;
%         fitData.r2(j,1) = r2(1,2);
%     end
end