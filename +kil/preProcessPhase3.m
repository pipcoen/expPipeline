function preProcessPhase3(animal,day,sites2Process)
% AP_preprocess_phase3(animal,day,testData)
dataPaths = {['\\zubjects.cortexlab.net\Subjects\' animal '\' day '\ephys']};
savePaths = {['\\zubjects.cortexlab.net\Subjects\' animal '\' day '\ephys\kilosort']};

% Check for multiple sites (based on the data path containing only folders)
dataPathDir = dir(dataPaths{1});
dataPathDir = dataPathDir(~contains({dataPathDir.name}', 'kilosort'));
if all([dataPathDir(3:end).isdir])
    dataPaths = cellfun(@(x) [dataPaths{1} filesep x],{dataPathDir(3:end).name},'uni',false);
    savePaths = cellfun(@(x) [savePaths{1} filesep x],{dataPathDir(3:end).name},'uni',false);
end

for currSite = find(sites2Process')
    currDataPath = dataPaths{currSite};
    currSavePath = savePaths{currSite};
    
    disp(['Processing data from ' currDataPath '...'])
    if ~exist(currSavePath,'dir'); mkdir(currSavePath); end
    if exist([currDataPath filesep 'experiment1'],'dir'); newStructure = 1; else, newStructure = 0; end
    
    
    % Filenames are semi-hardcoded in open ephys convention
    if ~newStructure
        apDataDir = dir([currDataPath '\experiment*_10*-0_0.dat']);
        syncDir = dir([currDataPath '\experiment*_all_channels_0.events']);
        settingsDir = dir([currDataPath '\settings*.xml']);
        
        apDataFileName = [currDataPath filesep apDataDir.name];
        syncFilename = [currDataPath filesep syncDir.name];
        settingsFilename = [currDataPath filesep settingsDir.name];
    else
        apDataFileName = [currDataPath '\experiment1\recording1\continuous\Neuropix-3a-100.0\continuous.dat'];
        syncTimestampsFilename = [currDataPath '\experiment1\recording1\events\Neuropix-3a-100.0\TTL_1\timestamps.npy' ];
        syncFilename = [currDataPath '\experiment1\recording1\events\Neuropix-3a-100.0\TTL_1\channel_states.npy'];
        settingsFilename = [currDataPath '\experiment1\recording1\structure.oebin'];
    end
    
    %% Get and save recording parameters
   
    % Get sample rate, gain, cutoff frequency (separate numbers and suffixes)
    if  ~newStructure
        ephysSettings = kil.xml2struct(settingsFilename);
        apGain = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.apGainValue,'%d%s');
        lfpGain = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.lfpGainValue,'%d%s');
        filterCut = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.filterCut,'%d%s');
    else
        apGain = {500};
        lfpGain = {125};
        filterCut = {300};
    end
    
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
    for currParam = 1:size(params,1)
        fprintf(fid,formatSpec,params{currParam,:});
    end
    fclose(fid);
    
    %% Get/save digital input events
    % Get/save digital input event times,
    if ~newStructure
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
    else
        % This way of getting start times is garbage too, assume AP
        % band is listed third (1 = software, 2 = AP, 3 = LFP)
        % (unused at the moment - I guess that assumes start time = 0)
        % Get/save digital input event times,
        syncData = readNPY(syncFilename);
        syncTimestamps = readNPY(syncTimestampsFilename);
        syncTimestamps = double(syncTimestamps)/apSampleRate;
        
        syncChannels = unique(abs(syncData));
        sync = struct('timestamps',cell(size(syncChannels)),'values',cell(size(syncChannels)));
        for currSync = 1:length(syncChannels)
            syncEvents = abs(syncData) == (syncChannels(currSync));
            sync(currSync).timestamps = syncTimestamps(syncEvents);
            sync(currSync).values = sign(syncData(syncEvents)) == 1;
        end
        
        syncSaveFilename = [currSavePath '\sync.mat'];
        save(syncSaveFilename,'sync');
    end
    
    
    %% Run kilosort
    
    % Set up local directory and clear out
    kilosortPath = 'D:\Temp\kilosort'; 
    apTempFilename = [animal '_' day  '_' 'ephys_apband.dat'];
    localPhyPath = 'D:\Temp\phy';
    if exist(kilosortPath, 'dir') 
        fileList = dir(kilosortPath);
        fileList = fileList(~ismember({fileList.name}, {'.', '..',apTempFilename}));
        arrayfun(@(x) delete([fileList(x).folder '/' fileList(x).name]), 1:length(fileList));
    else, mkdir(kilosortPath);
    end
    
    % Clear out whatever's currently in phy (usually not enough room)
    if exist(localPhyPath, 'dir'); rmdir(localPhyPath,'s'); end; mkdir(localPhyPath);
    
    % Copy data locally
    disp('Copying data to local drive...')
    apTempFilename = [kilosortPath filesep animal '_' day  '_' 'ephys_apband.dat'];
    if ~exist(apTempFilename, 'file')
        if exist(kilosortPath, 'dir'); rmdir(kilosortPath,'s'); end; mkdir(kilosortPath);
        copyfile(apDataFileName,apTempFilename); 
    end
    disp('Done'); 
   
    % Subtract common median across AP-band channels (hardcode channels?)
    ops.NchanTOT = 384;
    medianTrace = applyCARtoDat(apTempFilename, ops.NchanTOT);
%     %%
    apTempCarFilenameInit = [apTempFilename(1:end-4) '_CAR.dat'];
    apTempCarFilename = ['D' apTempCarFilenameInit(2:end)];
% 
%     % Get rid of the original non-CAR (usually not enough disk space)
    delete(apTempFilename);
    
%     disp('Moving CAR file before processing');
%     %%
%     movefile(apTempCarFilenameInit, apTempCarFilename);
    
    % Run kilosort on CAR data
    %%
    if ~exist('tRange','var'); tRange = [0,inf]; end
    kil.runKilosort(apTempCarFilename, apSampleRate, tRange);
    
    
    % Copy kilosort results to server
    
    disp('Copying sorted data to server...');
    resultsPath = [kilosortPath '\results'];
    copyfile(resultsPath,currSavePath);
    
    % Copy kilosort results and raw data to phy folder for clustering
    % Clear out whatever's currently in phy
    rmdir(localPhyPath,'s');
    mkdir(localPhyPath);
    
    % Copy the CAR'd data
    [~,apFile,apExt] = fileparts(apTempCarFilename);
    movefile(apTempCarFilename,[localPhyPath filesep apFile apExt])
    
    % Copy the results
    movefile([resultsPath filesep '*'],localPhyPath)
    
    % Delete all temporarly local data
    if exist(kilosortPath, 'dir'); rmdir(kilosortPath,'s'); end
    if exist(kilosortPath, 'dir'); rmdir(kilosortPath,'s'); end 
end

disp('Done processing phase 3 data.');


