function powerSpectra = loadAndSpectrogramLFP(animal,day,sites2Process)
%% Load LFP
dataPaths = {prc.pathFinder('serverprobedata', animal, day, 1)};
savePaths = {[prc.pathFinder('serverprobedata', animal, day, 1) '\kilosort']};

% Check for multiple sites (based on the data path containing only folders)
dataPathDir = dir(dataPaths{1});
dataPathDir = dataPathDir(~contains({dataPathDir.name}', 'kilosort'));
if all([dataPathDir(3:end).isdir])
    dataPaths = cellfun(@(x) [dataPaths{1} filesep x],{dataPathDir(3:end).name},'uni',false);
    savePaths = cellfun(@(x) [savePaths{1} filesep x],{dataPathDir(3:end).name},'uni',false);
end
%%
for currSite = find(sites2Process')
    %%
    currDataPath = dataPaths{currSite};
    currSavePath = savePaths{currSite};
    
    headerFID = fopen([currSavePath '\dat_params.txt']);
    headerInfo = textscan(headerFID,'%s %s', 'delimiter',{' = '});
    fclose(headerFID);
    headerInfo = [headerInfo{1}'; headerInfo{2}'];
    header = struct(headerInfo{:});
    %%
    ephysSampleRate = str2double(header.apSampleRate);
    lfpSampleRate = str2double(header.lfpSampleRate);
    numChannels = str2double(header.numChannels);
    lfpFilename = [currDataPath '\experiment1\recording1\continuous\Neuropix-3a-100.1\continuous.dat'];

    spikeTimes = double(readNPY([currSavePath '\spike_times.npy']))./ephysSampleRate;
    totalLength = floor(max(spikeTimes));
    
    %%
    windowSize = 60;
    numOfSamplePoints = 10;
    buffer = 300;    
    samplePoints = round(buffer:(totalLength-buffer)/(numOfSamplePoints+1):(totalLength-buffer));
    fid = fopen(lfpFilename);
    numberOfBytes = 2;
    freqPoints = 1:200;
    
    %%
    powerSpectra = zeros(length(freqPoints),numChannels,numOfSamplePoints);
    fprintf('%s %s creating spectrogram of LFP channels ... \n', day,animal);
    for i = 1:length(samplePoints)
        if i == round(numOfSamplePoints*0.25); fprintf('25 percent done... \n'); end
        if i == round(numOfSamplePoints*0.50); fprintf('50 percent done... \n'); end
        if i == round(numOfSamplePoints*0.75); fprintf('75 percent done... \n'); end
        lfpLoadStart = (lfpSampleRate*samplePoints(i)*numChannels*numberOfBytes);
        fseek(fid,lfpLoadStart,'bof');
        loadedLFP = fread(fid,[numChannels,windowSize*lfpSampleRate],'int16');
        powerSpectrum = pwelch(loadedLFP', [], [], freqPoints, lfpSampleRate);
        powerSpectra(:,:,i) = powerSpectrum;
    end
    
    save([currSavePath '\lfpPowerSpectra.mat'], 'powerSpectra', 'freqPoints');
    fclose(fid);
end