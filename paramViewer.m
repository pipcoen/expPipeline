function updateParamChangeSpreadsheet(subject)



if ~exist('expDate', 'var'); expDate = 'all'; end
if ~exist('prmType', 'var'); prmType = 'all'; end
[blks, prms] = getFilesFromDates(subject, expDate, 'prmblo');

selectedParams = {};
selectedFields = {};
selectedDates = cellfun(@(x) datenum(x, 'yyyy-mm-dd'),{prms.expDate}');
if contains(lower(prmType), {'rewardsize'; 'all'})
    selectedParams = [selectedParams; [prms.rewardSize]];
    selectedFields = [selectedFields; 'rewardSize'];
end
if contains(lower(prmType), {'validtrials'; 'all'})
    selectedParams = [selectedParams; [prms.validTrials]];
    selectedFields = [selectedFields; 'validTrials'];
end
if contains(lower(prmType), {'conditionSetChange'; 'all'})
    [~, ~, setIdx] = uniquecell({blks.uniqueConditions}');
    tDat = [1 diff(setIdx)']~=0;
    tDat = tDat.*cumsum(tDat);
    tDat = (tDat - smooth(tDat,3)').*([1 diff(setIdx)']~=0);
    selectedParams = [selectedParams; tDat];
    selectedFields = [selectedFields; 'conditionSets'];
end

t = find([1 diff(setIdx)']~=0);
t2 = diff([t length(setIdx)]);
for i = 1:length(t)
params{i,8} = blks(t(i)).uniqueConditions;
params{i,4} = size(blks(t(i)).uniqueConditions,1);
params{i,1} = blks(t(i)).expDate;
params{i,5} = t2(i);
params{i,6} = round(unique(blks(t(i)).uniqueConditions(:,1))*100);
params{i,2} = round(unique(blks(t(i)).uniqueConditions(:,2))*100);
params{i,3} = unique(blks(t(i)).trialType);
params{i,7} = unique(blks(t(i)).clickRate);
end

% colours = distinguishable_colors(length(selectedParams));
% aCol = repmat([246 158 219;149 211 79]./255,ceil(length(selectedParams)/2),1);
% figure('Name', [subject '''s parameters'], 'NumberTitle', 'Off', 'Units', 'Normalized');
% %%
% datacursormode on;
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',{@dispDatestr,selectedDates})
% pltH = (0.8/length(selectedParams));
% for i = 1:length(selectedParams)
%     minlim = min(selectedParams{i}(:))-abs((min(selectedParams{i}(:))/100)*length(selectedParams));
%     maxlim = max(selectedParams{i}(:))+abs((max(selectedParams{i}(:))/100)*length(selectedParams));
%     pltB = ((i-1)*pltH)+0.1;
%     axes('Position',[0.1 pltB 0.8 pltH]);
%     hold on
%     cSha = bsxfun(@rdivide, aCol(i,:), repmat((1:size(selectedParams{i},1))',[1,3]));
%     for j = 1:size(selectedParams{i},1)
%         stairs(selectedDates,selectedParams{i}(j,:),'LineWidth',2,'Color',cSha(j,:),'Marker','.');
%         stairs(selectedDates,selectedParams{i}(j,:),'LineWidth',2,'Color',cSha(j,:),'Marker','.');
%     end
%     text(0.01,0.8,selectedFields{i},'Units','normalized');
%     hold off
%     set(gca,'YTick',[],'YTickLabel',[],'YLim',[minlim maxlim],...
%     'XTick',[],'XTickLabel',[],'XLim',[min(selectedDates)-1 max(selectedDates)+1],'box','off');
%     if i ~= 1; continue; end
%     set(gca,'box','off','XTick',min(selectedDates):floor(length(selectedDates)/5):max(selectedDates),...
%         'XTickLabel',datestr(min(selectedDates):floor(length(selectedDates)/5):max(selectedDates)));
% end
% title([subject ': ' datestr(min(selectedDates),'yyyy-mm-dd') '�' datestr(max(selectedDates),'yyyy-mm-dd')],'FontWeight','bold');
% end
% %% Display date string
% function output_txt = dispDatestr(~,event_obj,uniqueUnits)
% % Display the position of the data cursor as date string
% % obj          Currently not used (empty)
% % event_obj    Handle to event object
% % output_txt   Data cursor text string (string or cell array of strings).
% 
% pos = get(event_obj,'Position');
% all_ax = get(gcf,'Children');
% all_ax = sort(all_ax(1:end-1));
% 
% output_txt = {['X: ',datestr(pos(1), 'yyyy-mm-dd')],...
%     ['Y: ',num2str(pos(2)),' ',uniqueUnits(all_ax==gca)]};
% 
% % If there is a Z-coordinate in the position, display it as well
%     if length(pos) > 2
%         output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
%     end
% end