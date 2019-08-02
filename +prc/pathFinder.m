function [pathOut, directoryCheck] = pathFinder(pathType, subject, expDate, expNum) 
%% A funciton to return various paths used in processing and anlysis. Changes in file struture should be reflected here.
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('subject', 'var'); [subject, expDate, expNum] = deal('noDataGiven'); end
if isstruct(subject) && ~exist('expDate', 'var') && ~exist('expNum', 'var')
    expDate = subject.expDate; 
    expNum = subject.expNum; 
    subject = subject.subject; 
end
if ~exist('expDate', 'var'); expDate = 'noDataGiven'; end
if ~exist('expNum', 'var'); expNum = 'noDataGiven'; end

if ~iscell(pathType); pathType = {pathType}; end
hostName = hostname;
if isnumeric(expNum); expNum = num2str(expNum); end
if isnumeric(expDate); expDate =  datestr(expDate, 'yyyy-mm-dd'); end

subjectPath = [subject '\' expDate '\' expNum '\'];
expRef = [expDate '_' expNum '_' subject];
processedFileName = [subject '\' subject '_' expDate([3:4 6:7 9:10]) '_' expNum  'Proc.mat'];

if contains(hostName, 'ziptop')
    driveName = 'C:';
else, driveName = 'D:';
end
% expInfo = '\\zubjects.cortexlab.net\Subjects\';
expInfo = {'\\zubjects.cortexlab.net\Subjects\'; '\\zserver.cortexlab.net\Data\Subjects\'};
if ~strcmp(expDate, 'noDataGiven')
    dateNumber = datenum(expDate, 'yyyy-mm-dd' );
    if dateNumber > 737589 && strcmp(subject, 'PC037'); expInfo = expInfo{2}; %'2019-06-13'
    elseif dateNumber > 737590 && strcmp(subject, 'PC038'); expInfo = expInfo{2}; %'2019-06-14'
    elseif dateNumber < 737612; expInfo = expInfo{1}; %'2019-07-06'
    else, expInfo = expInfo{2};
    end
end
rawBackup = [driveName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];
serverProcessedDirectory = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
    processedDirectory = [driveName '\Dropbox (Neuropixels)\MouseData\ProcessedData\'];
if contains('rawBlock', pathType); pathType{contains(pathType, 'rawBlock')} = 'backupBlock'; end
if contains('rawParams', pathType); pathType{contains(pathType, 'rawParams')} = 'backupParams'; end

if contains(hostName, {'homerig'; 'ziptop'})
    directoryCheck = 'local';
elseif strcmp(hostName, {'zip'})
    directoryCheck = 'all';
else; directoryCheck = 'server';
    processedDirectory = '\\zserver.cortexlab.net\lab\Share\Pip\ProcessedData\';
    if contains('rawBlock', pathType); pathType{contains(pathType, 'rawBlock')} = 'serverBlock'; end
    if contains('rawParams', pathType); pathType{contains(pathType, 'rawParams')} = 'serverParams'; end
end
%%
pathOut = pathType;
for i = 1:length(pathType)
switch lower(pathType{i})
    case 'serverblock'; pathOut{i,1} = [expInfo subjectPath expRef '_Block.mat'];
    case 'serverparams'; pathOut{i,1} = [expInfo subjectPath expRef '_Parameters.mat'];
    case 'serverfolder'; pathOut{i,1} = [expInfo subjectPath];
    case 'servertimeline'; pathOut{i,1} = [expInfo subjectPath expRef '_Timeline.mat'];
    case 'serverprobedata'; pathOut{i,1} = [expInfo subject '\' expDate '\ephys'];
    case 'serverprocesseddata'; pathOut{i,1} = [serverProcessedDirectory processedFileName];
    case 'serverprocessedfolder'; pathOut{i,1} = [serverProcessedDirectory subject];
    case 'serverprocesseddirectory'; pathOut{i,1} = serverProcessedDirectory;
    case 'serverfusi'; pathOut{i,1} = [expInfo subjectPath expRef '_fus.mat'];
    case 'backupblock'; pathOut{i,1} = [rawBackup subjectPath expRef '_Block.mat'];
    case 'backupparams'; pathOut{i,1} = [rawBackup subjectPath expRef '_parameters.mat'];
    case 'backupfolder'; pathOut{i,1} = [rawBackup subjectPath];
    case 'backupdirectory'; pathOut{i,1} = rawBackup;
    case 'processeddata'; pathOut{i,1} = [processedDirectory processedFileName];
    case 'processedfolder'; pathOut{i,1} = [processedDirectory subject];
    case 'processeddirectory'; pathOut{i,1} = processedDirectory;
    case 'galvolog'; pathOut{i,1} = [rawBackup subjectPath expRef '_galvoLog.mat'];
    case 'kilosortoutput'; pathOut{i,1} = [expInfo subject '\' expDate '\ephys\kilosort'];
    case 'explist'; pathOut{i,1} = [processedDirectory 'expList.mat'];
    case 'expinfo'; pathOut{i,1} = expInfo;
    case 'ephysrecord'; pathOut{i,1} = [processedDirectory 'ePhysRecord.xlsx'];
    case 'allenatlas'; pathOut{i,1} = [driveName '\Dropbox (Neuropixels)\MouseData\Atlas\allenCCF\'];
    case 'probepathdata'; pathOut{i,1} = [processedDirectory 'XHistology\' subject '\probe_points.mat'];
end
if length(pathOut) == 1; pathOut = pathOut{1}; end
end
