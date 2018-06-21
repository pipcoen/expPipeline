function updateParamChangeSpreadsheet(subject)
%%
recordPath = prc.pathFinder('mouserecord');
load(prc.pathFinder('expList'));
[blks, prms] = prc.getFilesFromDates(subject, 'all', 'prmblo');

%%
clear cache;
[~, fileInfo] = xlsfinfo(recordPath);
sheetValid = any(strcmp(fileInfo, subject));
%%
if sheetValid; [~, ~, xlsData] = xlsread(recordPath, subject);
else; disp(['MUST CREATE SHEET FOR ' subject]); return;
end

headers = xlsData(1:2,:);
xlsData(1:2,:) = [];

xlsData(cellfun(@(x) any(isnan(x)), xlsData(:,1)),:) = [];
laserTypes = cellfun(@(x) unique(double(x(:)))', {blks.laserType}', 'uni', 0);
laserPowers = cellfun(@(x) unique(double(x(:)))', {blks.laserPower}', 'uni', 0);
numberLaserSites = cellfun(@(x, y) size(x, 1)*any(y>0), {prms.galvoCoords}', laserPowers, 'uni', 0);
clickRate = cellfun(@unique, {prms.clickRate}', 'uni', 0);

combineParams = cellfun(@(w,x,y,z,a) [w repmat([x,y,z,a],size(w,1),1)],{blks.uniqueConditions}',laserTypes, laserPowers, numberLaserSites, clickRate, 'uni', 0);
[~, ~, setIdx] = uniquecell(combineParams);
changeIdx = find([1 diff(setIdx)']~=0)';

laserTypes = laserTypes(changeIdx);
laserPowers = laserPowers(changeIdx);
clickRate = clickRate(changeIdx);
numberLaserSites = numberLaserSites(changeIdx);
numberOfSessions = diff([changeIdx; length(setIdx)+1]);
datesInFile = cellfun(@(x) datenum(x, 'dd/mm/yyyy'), xlsData(:,1));
xlsData(:,1) = num2cell(datesInFile);
datesOfChange = cellfun(@(x) datenum(x, 'yyyy-mm-dd'), {blks(changeIdx).expDate}');
visContrasts = cellfun(@(x) unique(x(:,2))*100, {blks(changeIdx).uniqueConditions}', 'uni', 0);
audAmplitudes = cellfun(@(x) unique(x(:,1))*100, {blks(changeIdx).uniqueConditions}', 'uni', 0);
audInitialAzimuths = cellfun(@(x) unique(abs(x(:,3))), {blks(changeIdx).uniqueConditions}', 'uni', 0);
trialTypes = cellfun(@(x) unique(x(:)), {blks(changeIdx).trialType}', 'uni', 0);
numberOfConditions = cellfun(@(x) size(x,1), {blks(changeIdx).uniqueConditions}');

[datesInFile, uniqueIdx] = unique(datesInFile, 'stable');
xlsData = xlsData(uniqueIdx, :);
[~, excessDates] = setdiff(datesInFile, datesOfChange);
xlsData(excessDates, :) = [];
for i = 1:length(datesOfChange)
    idx = find(datesOfChange(i) == datesInFile);
    if isempty(idx); idx = size(xlsData, 1)+1; end
    xlsData{idx, 1} = datesOfChange(i);
    xlsData{idx, 3} = strtrim(sprintf('%g ',visContrasts{i}'));
    xlsData{idx, 4} = strtrim(sprintf('%g ',audAmplitudes{i}'));
    xlsData{idx, 5} = clickRate{i};
    xlsData{idx, 5+1} = strtrim(sprintf('%g ',audInitialAzimuths{i}'));
    xlsData{idx, 6+1} = num2str(trialTypes{i}');
    xlsData{idx, 7+1} = numberOfConditions(i);
    xlsData{idx, 8+1} = num2str(laserTypes{i});
    xlsData{idx, 9+1} = num2str(laserPowers{i}');
    xlsData{idx, 10+1} = num2str(numberLaserSites{i}');
    xlsData{idx, 11+1} = numberOfSessions(i);
end
for i = 3:size(xlsData,2)-1
    [~, ~, changeIdx] = uniquecell(xlsData(:,i));
    changeIdx = [0; diff(changeIdx)==0]>0;
    xlsData(changeIdx,i) = {'Same'};
end
xlsData = sortrows(xlsData,1);
xlsData(:,1) = cellfun(@(x) datestr(x, 'yyyy-mm-dd'), xlsData(:,1), 'uni', 0);
xlsData = [headers; xlsData];
xlsData(end+1:end+50, :) = {[]};
xlswrite(recordPath, xlsData, subject);
copyfile(recordPath, strrep(recordPath,'MouseData\','MouseData\ProcessedDataLite\'), 'f');

