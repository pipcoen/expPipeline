function [pathOut, directoryCheck] = pathFinder(pathType, pathInfo)
%% A funciton to return various paths used in processing and anlysis. Changes in file struture should be reflected here.
% INPUTS(default values)
% pathType(required)-------------A string to indicate the requested path. Can be a cell with multiple strings.
% pathInfo('noDataGiven')--------The subject information which can (but doesn't have to) contain the following fields:
%	.subject('noDataGiven')------------------Name of the subject
%	.expDate('noDataGiven')------------------Date of the experiment
%	.expNum('noDataGiven')-------------------Number of experiment
%	.datNum('noDataGiven')-------------------Number of experiment

% OUTPUTS
% pathOut---------------------------The requested path
% directoryCheck--------------------An indicator of whether the computer has 'server' or only 'local' access. Or 'all'.

%% Check inputs, extract variables from struture, and convert all to cells
if ~exist('pathType', 'var'); error('pathType required'); end
if ~exist('pathInfo', 'var'); [pathInfo.subject, pathInfo.expDate, pathInfo.expNum, pathInfo.datNum] = deal('noDataGiven'); end
if ~isfield(pathInfo, 'subject'); subject = 'noDataGiven'; else, subject = pathInfo.subject;  end
if ~isfield(pathInfo, 'expDate'); expDate = 'noDataGiven'; else, expDate = pathInfo.expDate;  end
if ~isfield(pathInfo, 'expNum'); expNum = repmat({'noDataGiven'}, length(subject),1); else, expNum = pathInfo.expNum;  end
if ~isfield(pathInfo, 'datNum'); datNum = repmat({'noDataGiven'}, length(subject),1); else, datNum = pathInfo.datNum;  end

if isnumeric(expNum); expNum = num2str(expNum); end
if isnumeric(expDate); expDate =  datestr(expDate, 'yyyy-mm-dd'); end

if ~iscell(pathType); pathType = {pathType}; end
if ~iscell(subject); subject = {subject}; end
if ~iscell(expDate); expDate = {expDate}; end
if ~iscell(expNum); expNum = {expNum}; end
if ~iscell(datNum); datNum = {datNum}; end

%% Make initial directory decisions based on dates and the computer that the program is running on.
%Assign the drive name and directoryCheck depending on where Pip keeps his dropbox
hostName = hostname;
if contains(hostName, 'ziptop'); driveName = 'C:';
else, driveName = 'D:';
end

% if contains(hostName, {'homerig'; 'ziptop'}); directoryCheck = 'local';
% elseif strcmp(hostName, {'zip'}); directoryCheck = 'all';
% else; directoryCheck = 'server';
% end
directoryCheck = 'server';

%Assign locations for the raw data, processed data etc. depending on the access of the computer being used.
rawBackup = [driveName '\Dropbox (Neuropixels)\MouseData\RawBehavior\'];
serverProcessedDirectory = '\\zserver.cortexlab.net\lab\Share\Magda\ProcessedMice\';
if strcmp(directoryCheck, 'server')
    processedDirectory = '\\zserver.cortexlab.net\lab\Share\Magda\ProcessedMice\';
    if contains('rawBlock', pathType); pathType{contains(pathType, 'rawBlock')} = 'serverBlock'; end
    if contains('rawParams', pathType); pathType{contains(pathType, 'rawParams')} = 'serverParams'; end
else
    processedDirectory = [driveName '\Dropbox (Neuropixels)\MouseData\ProcessedData4Paper\'];
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
        if strcmp(datNum{i}, 'noDataGiven'); datNum{i} = datenum(expDate{i}, 'yyyy-mm-dd'); end
        if datNum{i} > 737589 && strcmp(subject{i}, 'PC037'); expInfo = expInfo{2}; %'2019-06-13'
        elseif datNum{i} > 737590 && strcmp(subject{i}, 'PC038'); expInfo = expInfo{2}; %'2019-06-14'
        elseif datNum{i} < 737612; expInfo = expInfo{1}; %'2019-07-06'
        elseif datNum{i} > 737612 && datNum{i} <= 737739; expInfo = expInfo{2}; %'2019-07-06'
        elseif datNum{i} > 737739; expInfo = expInfo{3}; %'2019-11-10'
        else, expInfo = expInfo{2};g
        end
    end
    
    for j = 1:length(pathType)
        switch lower(pathType{j})                                                                             %hardcoded location of...

            case 'serverfolder'; pathOut{i,j} = [expInfo subjectPath];                                        %raw data folder on server
            case 'serverblock'; pathOut{i,j} = [expInfo subjectPath expRef '_Block.mat'];                     %raw block file on server
            case 'serverparams'; pathOut{i,j} = [expInfo subjectPath expRef '_Parameters.mat'];               %raw param file on server
            case 'servertimeline'; pathOut{i,j} = [expInfo subjectPath expRef '_Timeline.mat'];               %timeline file on server
            case 'serverprobedata'; pathOut{i,j} = [expInfo subject{i} '\' expDate{i} '\ephys'];              %probe data on server
            case 'serverprocesseddirectory'; pathOut{i,j} = serverProcessedDirectory;                         %directory for processed data server
            case 'serverprocessedfolder'; pathOut{i,j} = [serverProcessedDirectory subject{i}];               %folder for processed data server
            case 'serverprocesseddata'; pathOut{i,j} = [serverProcessedDirectory processedFileName];          %processed data on server

            case 'backupdirectory'; pathOut{i,j} = rawBackup;                                                 %local raw data directory
            case 'backupfolder'; pathOut{i,j} = [rawBackup subjectPath];                                      %local raw data folder
            case 'backupblock'; pathOut{i,j} = [rawBackup subjectPath expRef '_Block.mat'];                   %local raw block file
            case 'backupparams'; pathOut{i,j} = [rawBackup subjectPath expRef '_parameters.mat'];             %local raw param file
            case 'processeddirectory'; pathOut{i,j} = processedDirectory;                                     %local processed data directory
            case 'processedfolder'; pathOut{i,j} = [processedDirectory subject{i}];                           %local processed data folder
            case 'processeddata'; pathOut{i,j} = [processedDirectory processedFileName];                      %local processed data file

            case 'galvolog'; pathOut{i,j} = [expInfo subjectPath expRef '_galvoLog.mat'];                     %galvo records for inactivations
            case 'kilosortoutput'; pathOut{i,j} = [expInfo subject{i} '\' expDate{i} '\ephys\kilosort'];      %kilosort output on server
            case 'explist'; pathOut{i,j} = [processedDirectory 'expList.mat'];                                %the master list of experiments
            case 'expinfo'; pathOut{i,j} = expInfo;                                                           %the expInfo path
            case 'ephysrecord'; pathOut{i,j} = [processedDirectory 'ePhysRecord.mat'];                        %an excel sheet with ephys records
            case 'allenatlas'; pathOut{i,j} = [driveName '\Dropbox (Neuropixels)\MouseData\Atlas\allenCCF\']; %local allan atlas directory
            case 'probepath'; pathOut{i,j} = [processedDirectory 'XHistology\' subject{i} '\probe_histIdx' expNum{i}];  %probe vectors estimated from histology
        end
    end
end
if length(pathOut) == 1; pathOut = pathOut{1}; end
