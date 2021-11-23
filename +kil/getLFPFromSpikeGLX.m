function getLFPFromSpikeGLX(x)
% check which days from the mice's folder contain ephys data
apbinfiles = dir([fileparts(fileparts(x.serverFolder)) '\lfpSamp\**\*ap.bin*']); % get all the probe filenames
if isempty(apbinfiles)
    warning('No LFPSamp files detected...!');
end
%%
disp('Calculating LFP from lfp samples...');
nChan = 385;
lfpSpecta = struct;
probeNames = arrayfun(@(x) x.name((length(x.name)-11):(length(x.name)-7)), apbinfiles, 'uni', 0);
for myProbe = 1:size(apbinfiles,1)
    disp(['Sample ' num2str(myProbe) ' of ' num2str(size(apbinfiles,1)) '...']);
    myAPbin = apbinfiles(myProbe).name;
    myAPfolder = apbinfiles(myProbe).folder;    
    probeName = myAPbin((length(myAPbin)-11):(length(myAPbin)-7));
    probeSortedFolder = [x.kilosortOutput '\' probeName];
    d = dir([probeSortedFolder '\**\lfpPowerSpectra.mat']);
    if ~isempty(d); continue; end
    
    if ~isfield(lfpSpecta, probeName)
        lfpSpecta.(probeName).powerSpectra = []; 
        lfpSpecta.(probeName).xCoord = []; 
        lfpSpecta.(probeName).yCoord = []; 
    end
    myAPdata = [myAPfolder '\' apbinfiles(myProbe).name];
    
    myMeta = strrep(myAPdata, '.ap.bin', '.ap.meta');
    metaText = fileread(myMeta);
    expr = '[^\n]*imSampRate[^\n]*';
    matches = cell2mat(regexp(metaText,expr,'match'));    
    apSampleRate = str2double(matches(strfind(matches, '=')+1:end-1));

    [~, channelMapLoc] = kil.createMultiShankChannelMap(myAPdata,myAPfolder,1);
    load(channelMapLoc, 'xcoords', 'ycoords');
    
    nSamps = apbinfiles(myProbe).bytes/2/nChan;
    fsubsamp = 2500;
    lowPassCutoff = 300; % Hz
    [b1, a1] = butter(5, lowPassCutoff/apSampleRate, 'low');

    mmf = memmapfile(myAPdata,'Format',{'int16', [nChan nSamps],'x'});
    freqPoints = 1:200;
    windowSize = 10*fsubsamp;
    subData = zeros(10, windowSize);
    powerSpectrum = zeros(nChan-1, length(freqPoints));
    for i = 1:nChan-1
        chanData = double(mmf.Data.x(i,:));
        filtData = filtfilt(b1,a1, chanData);
        resampData = interp1(1:nSamps, filtData, 1:apSampleRate/fsubsamp:nSamps);
        
        for j = 1:10
            subData(j,:) = resampData(floor((j-1)*(length(resampData)-windowSize)/10)+(1:windowSize));
        end
         
        powerSpectrum(i,:) = mean(pwelch(subData', [], [], freqPoints, fsubsamp),2);
    end
    lfpSpecta.(probeName).powerSpectra = [lfpSpecta.(probeName).powerSpectra; powerSpectrum];
    lfpSpecta.(probeName).xCoord = [lfpSpecta.(probeName).xCoord; xcoords];
    lfpSpecta.(probeName).yCoord = [lfpSpecta.(probeName).yCoord; ycoords];
    
    if myProbe == find(contains(probeNames, probeName), 1, 'last')
        lfpPowerSpectra = lfpSpecta.(probeName);
        save([probeSortedFolder filesep 'lfpPowerSpectra.mat'], 'lfpPowerSpectra');
    end
end
