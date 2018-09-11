function preProcessPhase3(animal,day)
% AP_preprocess_phase3(animal,day,local_data)
savePaths = {['\\zubjects.cortexlab.net\Subjects\' animal '\' day '\ephys\kilosort']};
dataPaths = {['\\zubjects.cortexlab.net\Subjects\' animal '\' day '\ephys']};

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
    lfpDataDir = dir([currDataPath '\experiment*_10*-1_0.dat']);
    syncDir = dir([currDataPath '\experiment*_all_channels_0.events']);
    messagesDir = dir([currDataPath '\experiment*_messages_0.events']);
    settingsDir = dir([currDataPath '\settings*.xml']);
    
    apDataFileName = [currDataPath filesep apDataDir.name];
    lfpDataFileName = [currDataPath filesep lfpDataDir.name];
    syncFilename = [currDataPath filesep syncDir.name];
    messagesFilename = [currDataPath filesep messagesDir.name];
    settingsFilename = [currDataPath filesep settingsDir.name];
    
    %% Get and save recording parameters
    
    % Get index of electrophysiology channels in recordings
    ephysSettings = xml2struct(settingsFilename);
    
    % Get sample rate, gain, cutoff frequency (separate numbers and suffixes)
    apGain = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.apGainValue,'%d%s');
    lfpGain = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.lfpGainValue,'%d%s');
    filterCut = textscan(ephysSettings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.filterCut,'%d%s');
    
    % (0.195x for int16 to uV? how's this change with gain, just another x?)
    
    % Hard-coded parameters
    numChannels = 384;
    ap_sample_rate = 30000;
    lfp_sample_rate = 2500;
    
    params = {'raw_path',['''' currDataPath '''']; ...
        'numChannels',num2str(numChannels); ...
        'ap_sample_rate',num2str(ap_sample_rate); ... % this should be 30000 AP, 2500 LFP
        'lfp_sample_rate',num2str(lfp_sample_rate);
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
    
    % Get experiment start time (these messages are saved in a super dumb way
    % that are impossible to parse generally, so this is messy)
    messages_id = fopen(messagesFilename);
    messages_text = textscan(messages_id,'%*d %s %s', 'delimiter',{': '});
    fclose(messages_id);
    
    start_time_idx = strcmp(messages_text{1},'start time');
    start_time = str2num(messages_text{2}{start_time_idx}(1:strfind(messages_text{2}{start_time_idx},'@')-1));
    start_time_freq = str2num(messages_text{2}{start_time_idx}(strfind(messages_text{2}{start_time_idx},'@')+1: ...
        strfind(messages_text{2}{start_time_idx},'Hz')-1));
    start_time_sec = start_time/start_time_freq;
    
    % Get/save digital input event times,
    [sync_data, sync_timestamps, sync_info] = load_open_ephys_data_faster(syncFilename);
    sync_channels = unique(sync_data);
    sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
    for curr_sync = 1:length(sync_channels)
        sync_events = sync_data == (sync_channels(curr_sync));
        sync(curr_sync).timestamps = sync_timestamps(sync_events);
        sync(curr_sync).values = logical(sync_info.eventId(sync_events));
        
        % correct for experiment start time (not always necessary??)
        % as far as I can tell it's random whether this is needed or not: if
        % it's not then you get negative numbers at first, so maybe check for
        % those and then it can be automated? it's not a good sign that it's
        % variable though... I should probably just switch to spikeglx
%         if sync(curr_sync).timestamps(1) - start_time_sec > 0
%             sync(curr_sync).timestamps = sync(curr_sync).timestamps - start_time_sec;
%         end
%         sync(curr_sync).timestamps = sync(curr_sync).timestamps - start_time_sec;
    end
    
    sync_saveFilename = [currSavePath filesep 'sync.mat'];
    save(sync_saveFilename,'sync');
    
    %% Run kilosort
    
    % Set up local directory and clear out
    local_kilosort_path = 'E:\data_temp\kilosort';
    rmdir(local_kilosort_path,'s');
    mkdir(local_kilosort_path);
    
    % Clear out whatever's currently in phy (usually not enough room)
    local_phy_path = 'E:\data_temp\phy';
    rmdir(local_phy_path,'s');
    mkdir(local_phy_path);
    
    % Copy data locally
    disp('Copying data to local drive...')
    ap_tempFilename = [local_kilosort_path filesep animal '_' day  '_' 'ephys_apband.dat'];
    if ~exist(local_kilosort_path,'dir')
        mkdir(local_kilosort_path)
    end
    copyfile(apDataFileName,ap_tempFilename);
    disp('Done');
    
    % Subtract common median across AP-band channels (hardcode channels?)
    ops.NchanTOT = 384;
    medianTrace = applyCARtoDat(ap_tempFilename, ops.NchanTOT);
    ap_temp_carFilename = [ap_tempFilename(1:end-4) '_CAR.dat'];
    
    % Get rid of the original non-CAR (usually not enough disk space)
    delete(ap_tempFilename);
    
    % Run kilosort on CAR data
    AP_run_kilosort(ap_temp_carFilename,ap_sample_rate);
    
    
    %% Copy kilosort results to server
    
    disp('Copying sorted data to server...');
    ks_results_path = [local_kilosort_path filesep 'results'];
    copyfile(ks_results_path,currSavePath);
    
    %% Copy kilosort results and raw data to phy folder for clustering
    
    local_phy_path = 'E:\data_temp\phy';
    
    % Clear out whatever's currently in phy
    rmdir(local_phy_path,'s');
    mkdir(local_phy_path);
    
    % Copy the CAR'd data
    [~,ap_file,ap_ext] = fileparts(ap_temp_carFilename);
    movefile(ap_temp_carFilename,[local_phy_path filesep ap_file ap_ext])
    
    % Copy the results
    movefile([ks_results_path filesep '*'],local_phy_path)
    
    %% Delete all temporarly local data
    rmdir(local_kilosort_path,'s');
    mkdir(local_kilosort_path);
    
end

disp('Done processing phase 3 data.');


