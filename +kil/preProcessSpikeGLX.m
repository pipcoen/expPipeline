function preProcessSpikeGLX(x)
% this funtion takes a single ephys folder and sorts and CARs raw data within it
    ksTempFolder = ['D:\ksTemp\' x.subject '\' x.expDate]; % local folder on ssd where I process the data ==rootH 
    ksOutFolder = ['D:\ksProcessed\' x.subject '\' x.expDate]; % rootZ
    servFolder = fileparts(x.kilosortOutput);
    % check which days from the mice's folder contain ephys data
    if ~exist(ksTempFolder, 'dir'); mkdir(ksTempFolder); end
    if ~exist(ksOutFolder, 'dir'); mkdir(ksOutFolder); end    
    
    apbinfiles = dir([servFolder '\**\*ap.bin*']); % get all the probe filenames   
    
    %% remove filenames that are in the kilosort folder 
    while contains(apbinfiles(end).folder,'kilosort')==1
        apbinfiles(end)=[];
    end

    %%
    for myProbe = size(apbinfiles,1):-1:1
        myAPbin = apbinfiles(myProbe).name;
        myAPfolder = apbinfiles(myProbe).folder;
        probeName = myAPbin((length(myAPbin)-11):(length(myAPbin)-7)); 
        probeSortedFolder = [ksOutFolder '\' probeName]; 
        myAPdata = [myAPfolder '\' apbinfiles(myProbe).name]; 
        
        
        if ~exist([x.kilosortOutput '\' probeName '\spike_templates.npy'], 'file') ...
                && ~exist([probeSortedFolder '\spike_templates.npy'], 'file')
            fprintf('need to process %s file\n',myAPbin); 

            % create channelmap
            [~] = kil.createMultiShankChannelMap(myAPdata,myAPfolder,1);
            channelMapFile = dir([myAPfolder '\\**\\*_channelmap.mat*']);
            channelMapDir = [channelMapFile(1).folder '\' channelMapFile(1).name]; % channelmap for the probe - should be in the same folder
            
            if ~exist(probeSortedFolder, 'dir')
               mkdir(probeSortedFolder)
            end 
            
            if ~exist([probeSortedFolder '\' myAPbin], 'file')
                disp('copying data to local SSD');          
                copyfile(myAPdata , probeSortedFolder);
                disp('copied data') 
            else 
                fprintf('Data already copied\n');
            end
            
            
            kil.runKilosortSpikeGLX(probeSortedFolder,ksTempFolder,channelMapDir)
            %delete([probeSortedFolder '\' myAPbin]);  delete if you also
            %CAR the data
            % copy all other output to znas
            fprintf('Copying results to server ... \n'); 
            if ~exist([x.kilosortOutput '\' probeName], 'dir')
                mkdir([x.kilosortOutput '\' probeName])
            end
            prc.syncfolder(probeSortedFolder,[x.kilosortOutput '\' probeName])
        else 
            fprintf('already processed %s file\n',myAPbin);
        end
    end
%% extract sync pulse if I don't CAR the data 
    for myProbe = 1:size(apbinfiles,1)    
        myAPbin = apbinfiles(myProbe).name;
        myAPdata = [apbinfiles(myProbe).folder '\' apbinfiles(myProbe).name]; 
        probeName = myAPbin((length(myAPbin)-11):(length(myAPbin)-7)); 
        probeSortedFolder = [x.kilosortOutput '\' probeName];
        d = dir([probeSortedFolder '\**\sync.mat']);
        if numel(d)<1            
            kil.syncFTSpikeGLX(myAPdata, 385, probeSortedFolder);
        else 
            disp('sync extracted already.');
        end
    end
    disp('Done processing phase data.');
end



