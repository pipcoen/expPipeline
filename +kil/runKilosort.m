function runKilosort(dataFilename, sampleRate, tRange)
% AP_run_kilosort(dataFilename,sampleRate)
%
% dataFilename = .dat flat binary file of all channels together
% 
% New version (old version in _old)
% Runs kilosort (modified from master_file_example_MOVEME)

%% Assume the data is already copied locally for now
[dataPath,~,~] = fileparts(dataFilename);

%% Run Kilosort

% Run config script to get options
kil.kilosortConfigIMECP3O2
tic;
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

% Run kilosort
rez = preprocessDataSub(ops);
rez = clusterSingleBatches(rez);
rez = learnAndSolve8b(rez);
rez = find_merges(rez, 1);
rez = splitAllClusters(rez, 1);
rez = splitAllClusters(rez, 0);
rez = set_cutoff(rez);

% Convert results to phy, save
resultsDir = [dataPath '\results'];
mkdir(resultsDir);

save([resultsDir '\rez.mat'], 'rez', '-v7.3');
rezToPhy(rez, resultsDir);


%% Delete temporary data on local SSD

disp('Deleting temporary files...')

delete(ops.fproc);

disp('Done.')




