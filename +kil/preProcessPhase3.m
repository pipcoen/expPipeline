function preProcessPhase3(animal,day,testData)
% AP_preprocess_phase3(animal,day,testData)
if exist('testData','var') && testData; dataPaths =  {'E:\testData'};
else, dataPaths = {['\\zubjects.cortexlab.net\Subjects\' animal '\' day '\ephys']};
end
savePaths = {['\\zubjects.cortexlab.net\Subjects\' animal '\' day '\ephys\kilosort']};

% Check for multiple sites (based on the data path containing only folders)
dataPathDir = dir(dataPaths{1});
if all([dataPathDir(3:end).isdir])
    dataPaths = cellfun(@(x) [dataPaths{1} filesep x],{dataPathDir(3:end).name},'uni',false);
    savePaths = cellfun(@(x) [savePaths{1} filesep x],{dataPathDir(3:end).name},'uni',false);
end

for currSite = 1:length(dataPaths)
    
    currDataPath = dataPaths{currSite};
    currSavePath = savePaths{currSite};
    
    if ~exist(currSavePath,'dir')
        mkdir(currSavePath)
    end
    
    % Filenames are semi-hardcoded in open ephys convention
    apDataDir = dir([currDataPath '\experiment*_10*-0_0.dat']);
    syncDir = dir([currDataPath '\experiment*_all_channels_0.events']);
    settingsDir = dir([currDataPath '\settings*.xml']);
    
    apDataFileName = [currDataPath filesep apDataDir.name];
    syncFilename = [currDataPath filesep syncDir.name];
    settingsFilename = [currDataPath filesep settingsDir.name];
    
    %% Get and save recording parameters
    
    % Get index of electrophysiology channels in recordings
    ephysSettings = kil.xml2struct(settingsFilename);
    
    % Get sample rate, gain, cutoff frequency (separate numbers and suffixes)
    apGain = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.apGainValue,'%d%s');
    lfpGain = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.lfpGainValue,'%d%s');
    filterCut = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.filterCut,'%d%s');
    
    % (0.195x for int16 to uV? how's this change with gain, just another x?)
    
    % Hard-coded parameters
    numChannels = 384;
    apSampleRate = 30000;
    lfpSampleRate = 2500;
    
    params = {'raw_path',['''' currDataPath '''']; ...
        'numChannels',num2str(numChannels); ...
        'apSampleRate',num2str(apSampleRate); ... % this should be 30000 AP, 2500 LFP
        'lfpSampleRate',num2str(lfpSampleRate);
        'apGain',num2str(apGain{1}); ...
        'lfpGain',num2str(lfpGain{1})
        'filterCutoff',num2str(filterCut{1})};
    
    paramFilename = [currSavePath filesep 'dat_params.txt'];
    
    formatSpec = '%s = %s\r\n';
    fid = fopen(paramFilename,'w');
    for curr_param = 1:size(params,1)
        fprintf(fid,formatSpec,params{curr_param,:});
    end
    fclose(fid);
    
    %% Get/save digital input events
    % Get/save digital input event times,
    [syncData, syncTimestamps, syncInfo] = load_open_ephys_data_faster(syncFilename);
    syncChannels = unique(syncData);
    sync = struct('timestamps',cell(size(syncChannels)),'values',cell(size(syncChannels)));
    for currSync = 1:length(syncChannels)
        syncEvents = syncData == (syncChannels(currSync));
        sync(currSync).timestamps = syncTimestamps(syncEvents);
        sync(currSync).values = logical(syncInfo.eventId(syncEvents));
    end
    
    syncSaveFilename = [currSavePath '\sync.mat'];
    save(syncSaveFilename,'sync');
    
    %% Run kilosort
    
    % Set up local directory and clear out
    KilosortPath = 'C:\Temp\kilosort';
    if exist(KilosortPath, 'dir'); rmdir(KilosortPath,'s'); end
    mkdir(KilosortPath);
    
    % Clear out whatever's currently in phy (usually not enough room)
    localPhyPath = 'C:\Temp\phy';
    if exist(localPhyPath, 'dir'); rmdir(localPhyPath,'s'); end
    mkdir(localPhyPath);
    
    % Copy data locally
    disp('Copying data to local drive...')
    apTempFilename = [KilosortPath filesep animal '_' day  '_' 'ephys_apband.dat'];
    if ~exist(KilosortPath,'dir')
        mkdir(KilosortPath)
    end
    copyfile(apDataFileName,apTempFilename);
    disp('Done');
    
    % Subtract common median across AP-band channels (hardcode channels?)
    ops.NchanTOT = 384;
    medianTrace = applyCARtoDat(apTempFilename, ops.NchanTOT);
    apTempCarFilename = [apTempFilename(1:end-4) '_CAR.dat'];
    
    % Get rid of the original non-CAR (usually not enough disk space)
    delete(apTempFilename);
    
    % Run kilosort on CAR data
    kil.runKilosort(apTempCarFilename, apSampleRate);
    
    
    %% Copy kilosort results to server
    
    disp('Copying sorted data to server...');
    resultsPath = [KilosortPath '\results'];
    copyfile(resultsPath,currSavePath);
    
    %% Copy kilosort results and raw data to phy folder for clustering    
    % Clear out whatever's currently in phy
    rmdir(localPhyPath,'s');
    mkdir(localPhyPath);
    
    % Copy the CAR'd data
    [~,apFile,apExt] = fileparts(apTempCarFilename);
    movefile(apTempCarFilename,[localPhyPath filesep apFile apExt])
    
    % Copy the results
    movefile([resultsPath filesep '*'],localPhyPath)
    
    %% Delete all temporarly local data
    rmdir(KilosortPath,'s');
    mkdir(KilosortPath);
    
end

disp('Done processing phase 3 data.');


