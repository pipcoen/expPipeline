function pathObj = masterPath;
%% Set your main paths here
local = 1; %Keep a local backup of files somewhere to make loading etc faster


hostName = hostname;
if contains(hostName, 'ziptop'); driveName = 'C:';
else, driveName = 'D:';
end

if contains(hostName, {'homerig'; 'ziptop'}); directoryCheck = 'local';
elseif strcmp(hostName, {'zip'}); directoryCheck = 'all';
else; directoryCheck = 'server';
end

%Assign locations for the raw data, processed data etc. depending on the access of the computer being used.
rawBackup = [driveName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];
serverProcessedDirectory = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
if strcmp(directoryCheck, 'server')
    processedDirectory = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
    if contains('rawBlock', pathType); pathType{contains(pathType, 'rawBlock')} = 'serverBlock'; end
    if contains('rawParams', pathType); pathType{contains(pathType, 'rawParams')} = 'serverParams'; end
else
    processedDirectory = [driveName '\Dropbox (Neuropixels)\MouseData\ProcessedData\'];
    if contains('rawBlock', pathType); pathType{contains(pathType, 'rawBlock')} = 'backupBlock'; end
    if contains('rawParams', pathType); pathType{contains(pathType, 'rawParams')} = 'backupParams'; end
end

pathOut = cell(size(subject,1), length(pathType));
for i = 1:size(subject,1)
    %Set up paths based on the subject, date, etc.
    subjectPath = [subject{i} '\' expDate{i} '\' expNum{i} '\'];
    expRef = [expDate{i} '_' expNum{i} '_' subject{i}];
    processedFileName = [subject{i} '\' subject{i} '_' expDate{i}([3:4 6:7 9:10]) '_' expNum{i}  'Proc.mat'];
    
    %Annoying adjustments to account for changes in which server stored the data in lab. This is the expInfo path.
    expInfo = {'\\zubjects.cortexlab.net\Subjects\'; '\\zserver.cortexlab.net\Data\Subjects\'; '\\znas.cortexlab.net\Subjects\'};
    if ~strcmp(expDate{i}, 'noDataGiven')
        if datNum{i} > 737589 && strcmp(subject{i}, 'PC037'); expInfo = expInfo{2}; %'2019-06-13'
        elseif datNum{i} > 737590 && strcmp(subject{i}, 'PC038'); expInfo = expInfo{2}; %'2019-06-14'
        elseif datNum{i} < 737612 && datNum{i} <= 737739; expInfo = expInfo{1}; %'2019-07-06'
        elseif datNum{i} > 737739; expInfo = expInfo{3}; %'2019-11-10'
        end
    end
end
end