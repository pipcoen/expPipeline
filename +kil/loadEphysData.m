function loadEphysData(x)
%% Load timeline and associated inputs
if ~exist(x.rawTimeline, 'file'); error('No timeline file exists for requested ephys session'); end
fprintf('Loading timeline... \n');
timeline = load(x.rawTimeline); timeline = timeline.Timeline;
%%
inputNames = {timeline.hw.inputs.name}';
% Get acqLive signal
acqLiveTrace = timeline.rawDAQData(:,strcmp(inputNames, 'acqLive'));
acqLiveTrace = acqLiveTrace>(max(acqLiveTrace)/2);
acqLiveSrtEnd = timeline.rawDAQTimestamps(find(diff(acqLiveTrace))+1);

% Get wheel position %%%%NEVER USED?
% wheelPosition = timeline.rawDAQData(:,strcmp(inputNames, 'rotaryEncoder'));
% wheelPosition(wheelPosition > 2^31) = wheelPosition(wheelPosition > 2^31) - 2^32;

%Get whether stim was flickering %%%%NEVER USED?
% if contains('stimScreen', inputNames); stimScreenFlicker = range(timeline.rawDAQData(:,strcmp('stimScreen', inputNames))) > 2; end

% Get flipper flips
flipperTrace = timeline.rawDAQData(:,strcmp(inputNames, 'flipper')) > 2;
flipperFlips = sort([strfind(flipperTrace', [0 1]), strfind(flipperTrace', [1 0])])'+1;
flipperFlipTimelineTimes = timeline.rawDAQTimestamps(flipperFlips)';

%%%%%% Load mpep protocol: Ignoring for now. See AP_load_experiment for original %%%%%%%%%

%% Load task/behavior
fprintf('Converting block file... \n');
block  = load(x.processedData, 'blk'); block  = block.blk;

%% Load face/eyecam processing (with eyeGUI)

% Don't load if no timeline
if exist('Timeline','var') && load_parts.cam
    
    % Get cam sync from timeline
    camSync_idx = strcmp(inputNames, 'camSync');
    camSync_thresh = max(timeline.rawDAQData(:,camSync_idx))/2;
    camSync = timeline.rawDAQData(:,camSync_idx) > camSync_thresh;
    camSync_up = find((~camSync(1:end-1) & camSync(2:end)))+1;
    
    % EYECAM
    [eyecam_dir,eyecam_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam');
    
    if eyecam_exists
        if verbose; disp('Loading eyecam...'); end
        
        % Load camera processed data
        [eyecam_processed_filename,eyecam_processed_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam_processed');
        if eyecam_processed_exists
            eyecam = load(eyecam_processed_filename);
        end
        
        % Get camera times
        eyecam_fn = AP_cortexlab_filename(animal,day,experiment,'eyecam');
        eyecam_dir = fileparts(eyecam_fn);
        eyecam_t_savefile = [eyecam_dir filesep 'eyecam_t.mat'];
        
        if exist(eyecam_fn,'file') && ~exist(eyecam_t_savefile,'file')
            % Get facecam strobes
            eyeCamStrobe_idx = strcmp(inputNames, 'eyeCameraStrobe');
            eyeCamStrobe_thresh = max(timeline.rawDAQData(:,eyeCamStrobe_idx))/2;
            eyeCamStrobe = timeline.rawDAQData(:,eyeCamStrobe_idx) > eyeCamStrobe_thresh;
            eyeCamStrobe_up = find((~eyeCamStrobe(1:end-1) & eyeCamStrobe(2:end)))+1;
            eyeCamStrobe_up_t = timeline.rawDAQTimestamps(eyeCamStrobe_up);
            
            % Get sync times for cameras (or load if already done)
            [eyecam_sync_frames,n_eyecam_frames] = AP_get_cam_sync_frames(eyecam_fn);
            
            if ~isempty(eyecam_sync_frames)
                % Get the closest facecam strobe to sync start, find offset and frame idx
                [~,eyecam_strobe_sync] = min(abs(camSync_up(1) - eyeCamStrobe_up));
                eyecam_frame_offset = eyecam_sync_frames(1) - eyecam_strobe_sync;
                eyecam_frame_idx = [1:length(eyeCamStrobe_up)] + eyecam_frame_offset;
                
                % Get times of facecam frames in timeline
                eyecam_t = nan(n_eyecam_frames,1);
                eyecam_t(eyecam_frame_idx(eyecam_frame_idx > 0)) = eyeCamStrobe_up_t(eyecam_frame_idx > 0);
                
                save(eyecam_t_savefile,'eyecam_t');
            end
        elseif exist(eyecam_fn,'file') && exist(eyecam_t_savefile,'file')
            load(eyecam_t_savefile);
        end
        
    end
    
    % FACECAM
    [facecam_dir,facecam_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam');
    
    if facecam_exists
        if verbose; disp('Loading facecam...'); end
        
        [facecam_processed_filename,facecam_processed_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam_processed');
        if facecam_processed_exists
            facecam = load(facecam_processed_filename);
        end
        
        % Get camera times
        facecam_fn = AP_cortexlab_filename(animal,day,experiment,'facecam');
        facecam_dir = fileparts(facecam_fn);
        facecam_t_savefile = [facecam_dir filesep 'facecam_t.mat'];
        
        if exist(facecam_fn,'file') && ~exist(facecam_t_savefile,'file')
            % Get facecam strobes
            faceCamStrobe_idx = strcmp(inputNames, 'faceCamStrobe');
            faceCamStrobe_thresh = max(timeline.rawDAQData(:,faceCamStrobe_idx))/2;
            faceCamStrobe = timeline.rawDAQData(:,faceCamStrobe_idx) > faceCamStrobe_thresh;
            faceCamStrobe_up = find((~faceCamStrobe(1:end-1) & faceCamStrobe(2:end)))+1;
            faceCamStrobe_up_t = timeline.rawDAQTimestamps(faceCamStrobe_up);
            
            % Get sync times for cameras (or load if already done)
            [facecam_sync_frames,n_facecam_frames] = AP_get_cam_sync_frames(facecam_fn);
            
            if ~isempty(facecam_sync_frames)
                % Get the closest facecam strobe to sync start, find offset and frame idx
                [~,facecam_strobe_sync] = min(abs(camSync_up(1) - faceCamStrobe_up));
                facecam_frame_offset = facecam_sync_frames(1) - facecam_strobe_sync;
                facecam_frame_idx = [1:length(faceCamStrobe_up)] + facecam_frame_offset;
                
                % Get times of facecam frames in timeline
                facecam_t = nan(n_facecam_frames,1);
                facecam_t(facecam_frame_idx) = faceCamStrobe_up_t;
                
                save(facecam_t_savefile,'facecam_t');
            end
        elseif exist(facecam_fn,'file') && exist(facecam_t_savefile,'file')
            load(facecam_t_savefile);
        end
        
    end
    
end

%% Load imaging data

[data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'imaging',site);
experiment_path = [data_path filesep num2str(experiment)];

% (check for specific imaging file since data path is just root)
spatialComponents_fns = dir([data_path filesep 'svdSpatialComponents*']);
imaging_exists = ~isempty(spatialComponents_fns);

if imaging_exists && load_parts.imaging
    if verbose; disp('Loading imaging data...'); end
    
    % Get the imaging file locations
    spatialComponents_dir = dir([data_path filesep 'svdSpatialComponents*']);
    temporalComponents_dir = dir([experiment_path filesep 'svdTemporalComponents*']);
    imaging_timestamps_idx = cellfun(@any,strfind({temporalComponents_dir.name},'timestamps'));
    
    meanImage_dir = dir([data_path filesep 'meanImage*']);
    
    cam_color_n = length(spatialComponents_dir);
    cam_color_signal = 'blue';
    cam_color_hemo = 'purple';
    
    if cam_color_n == 1
        
        U = readUfromNPY([data_path filesep spatialComponents_dir.name]);
        V = readVfromNPY([experiment_path filesep temporalComponents_dir(~imaging_timestamps_idx).name]);
        frame_t = readNPY([experiment_path filesep temporalComponents_dir(imaging_timestamps_idx).name]);
        
        framerate = 1./nanmean(diff(frame_t));
        
        % Detrend and high-pass filter
        highpassCutoff = 0.01; % Hz
        [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
        dV = detrend(V', 'linear')';
        fV = single(filter(b100s,a100s,double(dV)')');
        
        avg_im = readNPY([data_path filesep meanImage_dir.name]);
        
    elseif cam_color_n == 2
        
        % Load in all things as neural (n) or hemodynamic (h)
        
        tn = readNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_signal '.timestamps.npy']);
        Un = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_signal '.npy']);
        Vn = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_signal '.npy']);
        dataSummary_n = load([data_path filesep 'dataSummary_' cam_color_signal '.mat']);
        avg_im_n = readNPY([data_path filesep 'meanImage_' cam_color_signal '.npy']);
        
        th = readNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.timestamps.npy']);
        Uh = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_hemo '.npy']);
        Vh = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.npy']);
        dataSummary_h = load([data_path filesep 'dataSummary_' cam_color_signal '.mat']);
        avg_im_h = readNPY([data_path filesep 'meanImage_' cam_color_hemo '.npy']);
        
        framerate = 1./nanmean(diff(tn));
        
        % Correct hemodynamic signal in blue from green
        % First need to shift alternating signals to be temporally aligned
        % (shifts neural to hemo)
        if verbose; disp('Correcting hemodynamics...'); end
        
        % Check if number of timeline frames matches imaged frames
        cam_tl_imaged_diff = length(cam_time) - (size(Vn,2) + size(Vh,2));
        if cam_tl_imaged_diff ~= 0
            warning(sprintf( ...
                '\n %s %s: %d timeline-imaged frames: assuming dropped at end', ...
                animal,day,cam_tl_imaged_diff));
        end
        
        % Eliminate odd frames out (unpaired colors)
        min_frames = min(size(Vn,2),size(Vh,2));
        Vn = Vn(:,1:min_frames);
        tn = tn(1:min_frames);
        
        Vh = Vh(:,1:min_frames);
        th = th(1:min_frames);
        
        % This was to get rid of bad exposures: not sure I want this though
        %         cam_expose_time_reshape = ...
        %             reshape(cam_expose_times(1:end-mod(length(cam_expose_times),2)),2,[]);
        %         bad_cam_expose = any(cam_expose_time_reshape > ...
        %             median(cam_expose_time_reshape(:))*2,1);
        %
        %         if any(bad_cam_expose)
        %             warning(['Bad cam expose time: ' num2str(find(bad_cam_expose)) '/' num2str(min_frames)]);
        %             Vn = Vn(:,~bad_cam_expose);
        %             tn = tn(~bad_cam_expose);
        %
        %             Vh = Vh(:,~bad_cam_expose);
        %             th = th(~bad_cam_expose);
        %         end
        
        Vn_th = SubSampleShift(Vn,1,2);
        
        Vh_Un = ChangeU(Uh,Vh,Un);
        
        hemo_tform_fn = [experiment_path filesep 'hemo_tform.mat'];
        if exist(hemo_tform_fn,'file')
            % If the hemo tform matrix has been computed, load and fix
            if verbose; disp('Using old hemo tform...'); end
            load(hemo_tform_fn)
            zVh_Un = bsxfun(@minus, Vh_Un, nanmean(Vh_Un,2));
            Vn_hemo = transpose(Vn_th' - zVh_Un'*hemo_tform');
        else
            % If no p hemo tform matrix, compute and save
            if verbose; disp('Computing hemo tform...'); end
            %hemo_freq = [0.1,1];
            hemo_freq = [7,13];
            [Vn_hemo,hemo_tform] = HemoCorrectLocal(Un,Vn_th,Vh_Un,framerate,hemo_freq,3);
            save(hemo_tform_fn,'hemo_tform');
            % Close the figures (hacky - but function isn't mine)
            close(gcf)
            close(gcf)
        end
        
        if verbose; disp('Filtering...'); end
        % Don't bother filtering heartbeat, just detrend and highpass
        % fVn_hemo = detrendAndFilt(Vn_hemo, framerate);
        highpassCutoff = 0.01; % Hz
        [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
        
        dVn_hemo = detrend(Vn_hemo', 'linear')';
        
        % non-zero-lag filter, but causal (only moves forwards in time)
        fVn_hemo = filter(b100s,a100s,dVn_hemo,[],2);
        % non-causal but zero-lag filter: changed because can introduce
        % artifacts with single wonky points, also big changes propogate
        % backwards in time which potentially gives bad causality
        %fVn_hemo = single(filtfilt(b100s,a100s,double(dVn_hemo)')');
        
        % Do this for the colors individually, in case they're used
        dVn = detrend(Vn', 'linear')';
        fVn = single(filter(b100s,a100s,double(dVn)')');
        
        dVh = detrend(Vh', 'linear')';
        fVh = single(filter(b100s,a100s,double(dVh)')');
        
        % set final U/V to use
        fV = fVn_hemo;
        U = Un;
        avg_im = avg_im_n;
        frame_t = th; % shifted to use hemo color times
        
    end
    if verbose; disp('Done.'); end
    
    % Make dF/F
    [Udf,fVdf] = dffFromSVD(U,fV,avg_im);
    % zero out NaNs in the Udfs (from saturated pixels?)
    Udf(isnan(Udf)) = 0;
end

%% Load ephys data (single long recording)
[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys',site);

if ephys_exists && load_parts.ephys
    
    if verbose; disp('Loading ephys...'); end
    
    % These are the digital channels going into the FPGA
    photodiode_sync_idx = 1;
    acqLive_sync_idx = 2;
    pcoExposure_sync_idx = 3;
    flipper_sync_idx = 4;
    
    load_lfp = false;
    
    % Load clusters, if they exist
    cluster_filename = [ephys_path filesep 'cluster_groups.csv'];
    if exist(cluster_filename,'file')
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
    end
    
    % Apparently now sometimes it's a different filename/type, if that
    % exists overwrite the other one
    cluster_filename = [ephys_path filesep 'cluster_group.tsv'];
    if exist(cluster_filename,'file')
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
    end
    
    % Load sync/photodiode
    load(([ephys_path filesep 'sync.mat']));
    
    % Read header information
    header_path = [ephys_path filesep 'dat_params.txt'];
    header_fid = fopen(header_path);
    header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
    fclose(header_fid);
    
    header = struct;
    for i = 1:length(header_info{1})
        header.(header_info{1}{i}) = header_info{2}{i};
    end
    
    % Load spike data
    if isfield(header,'sample_rate')
        ephys_sample_rate = str2num(header.sample_rate);
    elseif isfield(header,'ap_sample_rate')
        ephys_sample_rate = str2num(header.ap_sample_rate);
    end
    spike_times = double(readNPY([ephys_path filesep 'spike_times.npy']))./ephys_sample_rate;
    spike_templates = readNPY([ephys_path filesep 'spike_templates.npy']);
    templates = readNPY([ephys_path filesep 'templates.npy']);
    channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
    channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
    winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
    template_amplitudes = readNPY([ephys_path filesep 'amplitudes.npy']);
    
    % Flip channel map and positions if banks are reversed
    % (this was only for phase 2, so setting false by default)
    flipped_banks = false;
    if flipped_banks
        channel_map = [channel_map(61:end);channel_map(1:60)];
        channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
    end
    
    % Default channel map/positions are from end: make from surface
    channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);
    
    % Load LFP
    n_channels = str2num(header.n_channels);
    %lfp_filename = [ephys_path filesep 'lfp.dat']; (this is old)
    [data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'ephysraw',site);
    lfp_dir = dir([data_path 'experiment*-1_0.dat']);
    lfp_filename = [data_path lfp_dir.name];
    if load_lfp && exist(lfp_filename,'file')
        lfp_sample_rate = str2num(header.lfp_sample_rate);
        lfp_cutoff = str2num(header.filter_cutoff);
        
        fid = fopen(lfp_filename);
        % define where/how much of LFP to load
        lfp_skip_minutes = 10; % move to N minutes after recording start
        lfp_load_start = (lfp_sample_rate*60*lfp_skip_minutes*n_channels);
        lfp_load_samples = 1e6;
        % load LFP
        fseek(fid,lfp_load_start,'bof');
        lfp_all = fread(fid,[n_channels,lfp_load_samples],'int16'); % pull snippet
        fclose(fid);
        % eliminate non-connected channels
        lfp = lfp_all(channel_map+1,:);
        clear lfp_all;
        
        lfp_t = [(lfp_load_start/n_channels):(lfp_load_start/n_channels)+lfp_load_samples-1]/lfp_sample_rate;
    end
    
    % Get sync points for alignment
    
    % Get experiment index by finding numbered folders
    protocols = AP_list_experiments(animal,day);
    experiment_idx = experiment == [protocols.experiment];
    
    if exist('flipperFlipTimelineTimes','var')
        % (if flipper, use that)
        % (at least one experiment the acqLive connection to ephys was bad
        % so it was delayed - ideally check consistency since it's
        % redundant)
        bad_flipper = false;
        
        % Get flipper experiment differences by long delays
        flip_diff_thresh = 1; % time between flips to define experiment gap (s)
        flipper_expt_idx = [1;find(diff(sync(flipper_sync_idx).timestamps) > ...
            flip_diff_thresh)+1;length(sync(flipper_sync_idx).timestamps)+1];
        
        flipper_flip_times_ephys = sync(flipper_sync_idx).timestamps( ...
            flipper_expt_idx(find(experiment_idx)):flipper_expt_idx(find(experiment_idx)+1)-1);
        
        % Check that number of flipper flips in timeline matches ephys
        if length(flipper_flip_times_ephys) ~= length(flipperFlipTimelineTimes)
            warning([animal ' ' day ':Flipper flip times different in timeline/ephys'])
            bad_flipper = true;
        end
        
        %         % Plot stim aligned by flipper times for sanity check
        %         figure; hold on;
        %         stim_ephys_timestamps = interp1(flipper_flip_times_ephys,flipperFlipTimelineTimes,sync(1).timestamps,'linear','extrap');
        %         plot(stim_ephys_timestamps,sync(1).values*5,'.g')
        %         plot(stimOn_times,5,'ob');
        
        sync_timeline = flipperFlipTimelineTimes;
        sync_ephys = flipper_flip_times_ephys;
    end
    
    if ~exist('flipperFlipTimelineTimes','var') || bad_flipper
        % (if no flipper or flipper problem, use acqLive)
        
        % Get acqLive times for current experiment
        experiment_ephys_starts = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 1);
        experiment_ephys_stops = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 0);
        acqlive_ephys_currexpt = [experiment_ephys_starts(experiment_idx), ...
            experiment_ephys_stops(experiment_idx)];
        
        sync_timeline = acqLiveSrtEnd;
        sync_ephys = acqlive_ephys_currexpt;
        
        % Check that the experiment time is the same within threshold
        % (it should be almost exactly the same)
        if abs(diff(acqLiveSrtEnd) - diff(acqlive_ephys_currexpt)) > 1
            error([animal ' ' day ': acqLive duration different in timeline and ephys']);
        end
    end
    
    % Get the spike/lfp times in timeline time (accounts for clock drifts)
    spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times,'linear','extrap');
    if load_lfp && exist(lfp_filename,'file')
        lfp_t_timeline = interp1(sync_ephys,sync_timeline,lfp_t,'linear','extrap');
    end
    
    % Get the depths of each template
    % (by COM - this used to not work but now looks ok)
    [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
        templatePositionsAmplitudes(templates,winv,channel_positions(:,2),spike_templates,template_amplitudes);
    %     % (by max waveform channel)
    %     template_abs = permute(max(abs(templates),[],2),[3,1,2]);
    %     [~,max_channel_idx] =  max(template_abs,[],1);
    %     templateDepths = channel_positions(max_channel_idx,2);
    %     % Get each spike's depth
    %     spikeDepths = templateDepths(spike_templates+1);
    
    % Get the waveform duration of all templates (channel with largest amp)
    [~,max_site] = max(max(abs(templates),[],2),[],3);
    templates_max = nan(size(templates,1),size(templates,2));
    for curr_template = 1:size(templates,1)
        templates_max(curr_template,:) = ...
            templates(curr_template,:,max_site(curr_template));
    end
    waveforms = templates_max;
    
    % Get trough-to-peak time for each template
    templates_max_signfix = bsxfun(@times,templates_max, ...
        sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));
    
    [~,waveform_trough] = min(templates_max,[],2);
    [~,waveform_peak_rel] = arrayfun(@(x) ...
        max(templates_max(x,waveform_trough(x):end),[],2), ...
        transpose(1:size(templates_max,1)));
    waveform_peak = waveform_peak_rel + waveform_trough;
    
    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;
    
    % Eliminate spikes that were classified as not "good"
    if exist('cluster_groups','var')
        
        if verbose; disp('Removing non-good/MUA templates'); end
        
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
        
        % Throw out all non-good template data
        templates = templates(good_templates,:,:);
        templateDepths = templateDepths(good_templates);
        waveforms = waveforms(good_templates,:);
        templateDuration = templateDuration(good_templates);
        templateDuration_us = templateDuration_us(good_templates);
        
        % Throw out all non-good spike data
        good_spike_idx = ismember(spike_templates,good_templates_idx);
        spike_times = spike_times(good_spike_idx);
        spike_templates = spike_templates(good_spike_idx);
        template_amplitudes = template_amplitudes(good_spike_idx);
        spikeDepths = spikeDepths(good_spike_idx);
        spike_times_timeline = spike_times_timeline(good_spike_idx);
        
        % Re-name the spike templates according to the remaining templates
        % (and make 1-indexed from 0-indexed)
        new_spike_idx = nan(max(spike_templates)+1,1);
        new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
        spike_templates = new_spike_idx(spike_templates+1);
        
    elseif ~exist('cluster_groups','var')
        if verbose; disp('Clusters not yet sorted'); end
    end
    
end

%% Regress light artifact out of LFP (unused at the moment)

% if ephys_exists && load_parts.ephys
%
%     if verbose; disp('Cleaning LFP...'); end;
%
%     % LFP correlation
%     % (fix the light artifact on each channel with regression)
%     light_timeline = interp1(acqlive_ephys_currexpt,acqLiveSrtEnd,sync(3).timestamps,'linear','extrap');
%     light_on = light_timeline(sync(3).values == 1);
%     light_off = light_timeline(sync(3).values == 0);
%
%     blue_on = light_on(1:2:end);
%     blue_off = light_off(1:2:end);
%     violet_on = light_on(2:2:end);
%     violet_off = light_off(2:2:end);
%
%     lfp_t_bins = [lfp_t_timeline-0.5/lfp_sample_rate,lfp_t_timeline(end)+0.5/lfp_sample_rate];
%     blue_on_vector = histcounts(blue_on,lfp_t_bins);
%     blue_off_vector = histcounts(blue_off,lfp_t_bins);
%     violet_on_vector = histcounts(violet_on,lfp_t_bins);
%     violet_off_vector = histcounts(violet_off,lfp_t_bins);
%
%     light_vectors = [blue_on_vector;blue_off_vector;violet_on_vector;violet_off_vector];
%
%     t_shift = round((1/35)*lfp_sample_rate*1.5);
%     t_shifts = [-t_shift:t_shift];
%     lambda = 0;
%     zs = [false,false];
%     cvfold = 1;
%
%     % (in chunks: necessary memory-wise, also allows changing light)
%     n_chunks = 10;
%     lfp_t_chunk = round(linspace(1,size(lfp,2),n_chunks+1));
%
%     lfp_lightfix = nan(size(lfp));
%     for curr_chunk = 1:n_chunks
%         curr_chunk_t = lfp_t_chunk(curr_chunk):lfp_t_chunk(curr_chunk+1);
%         [light_k,artifact_lfp] = AP_regresskernel(light_vectors(:,curr_chunk_t),lfp(:,curr_chunk_t),t_shifts,lambda,zs,cvfold);
%
%         lfp_lightfix(:,curr_chunk_t) = lfp(:,curr_chunk_t)-artifact_lfp;
%         %             AP_print_progress_fraction(curr_chunk,n_chunks);
%     end
%     lfp_lightfix(isnan(lfp_lightfix)) = 0;
%
%     % (group channels by depth)
%     channel_depth_grp = discretize(channel_positions(:,2),depth_group_edges);
%     lfp_depth_median = grpstats(lfp_lightfix,channel_depth_grp,'median');
%
%     % (low-pass filter: sometimes bunch of junk at high freq?)
%     freqCutoff = 300; % Hz
%     [b100s, a100s] = butter(2,freqCutoff/(lfp_sample_rate/2),'low');
%     lfp_depth_median_filt = single(filtfilt(b100s,a100s,double(lfp_depth_median)')');
%
% end

%% Estimate striatal boundaries on probe

if ephys_exists && load_parts.ephys
    if verbose; disp('Estimating striatum boundaries on probe...'); end
    
    % str_align = alignment method ('none', 'depth', or 'kernel')
    
    % requires n_aligned_depths for alignment, default 4
    if ~exist('n_aligned_depths','var')
        n_aligned_depths = 4;
    end
    
    % if no alignment specified, default kernel
    if ~exist('str_align','var')
        str_align = 'kernel';
    end
    
    AP_align_striatum_ephys
    
    if verbose
        figure;
        imagesc(depth_group_centers,depth_group_centers,mua_corr);
        axis tight equal;
        colormap(hot)
        line([str_depth(1),str_depth(1)],ylim,'color','b','linewidth',2);
        line([str_depth(2),str_depth(2)],ylim,'color','b','linewidth',2);
        line(xlim,[str_depth(1),str_depth(1)],'color','b','linewidth',2);
        line(xlim,[str_depth(2),str_depth(2)],'color','b','linewidth',2);
        xlabel('Probe depth (\mum)');
        ylabel('Probe depth (\mum)');
        title('MUA correlation: striatum location');
    end
    
end

%% Classify spikes

if ephys_exists && load_parts.ephys && exist('cluster_groups','var')
    if verbose; disp('Classifying spikes...'); end
    
    str_templates = templateDepths >= str_depth(1) & templateDepths <= str_depth(2);
    non_str_templates = ~str_templates;
    
    % Define the window to look for spiking statistics in (spikes go in and
    % out, so take the bin with the largest firing rate for each cell and work
    % with that one)
    % spiking_stat_window = 60*5; % seconds
    % spiking_stat_bins = min(spike_times_timeline):spiking_stat_window: ...
    %     max(spike_times_timeline);
    
    % % (for whole session)
    spiking_stat_window = max(spike_times_timeline)-min(spike_times_timeline);
    spiking_stat_bins = [min(spike_times_timeline),max(spike_times_timeline)];
    
    % Get firing rate across the session
    bin_spikes = nan(max(spike_templates), ...
        length(spiking_stat_bins)-1);
    for curr_template = unique(spike_templates)'
        bin_spikes(curr_template,:) = ...
            histcounts(spike_times_timeline(spike_templates == curr_template), ...
            spiking_stat_bins);
    end
    min_spikes = 10;
    use_spiking_stat_bins = bsxfun(@ge,bin_spikes,prctile(bin_spikes,80,2)) & bin_spikes > min_spikes;
    spike_rate = sum(bin_spikes.*use_spiking_stat_bins,2)./ ...
        (sum(use_spiking_stat_bins,2)*spiking_stat_window);
    
    % Get proportion of ISI > 2s (Yamin/Cohen 2013) and CV2 (Stalnaker/Schoenbaum 2016)
    prop_long_isi = nan(max(spike_templates),1);
    cv2 = nan(max(spike_templates),1);
    for curr_template = unique(spike_templates)'
        
        long_isi_total = 0;
        isi_ratios = [];
        for curr_bin = find(use_spiking_stat_bins(curr_template,:))
            curr_spike_times = spike_times_timeline( ...
                spike_times_timeline > spiking_stat_bins(curr_bin) & ...
                spike_times_timeline < spiking_stat_bins(curr_bin+1) & ...
                spike_templates == curr_template);
            curr_isi = diff(curr_spike_times);
            
            long_isi_total = long_isi_total + sum(curr_isi(curr_isi > 2));
            
            isi_ratios = [isi_ratios;(2*abs(curr_isi(2:end) - curr_isi(1:end-1)))./ ...
                (curr_isi(2:end) + curr_isi(1:end-1))];
        end
        
        prop_long_isi(curr_template) = long_isi_total/ ...
            (sum(use_spiking_stat_bins(curr_template,:))*spiking_stat_window);
        cv2(curr_template) = nanmean(isi_ratios);
        
    end
    
    % Cortical classification (like Bartho JNeurophys 2004)
    waveform_duration_cutoff = 400;
    narrow = non_str_templates & templateDuration_us <= waveform_duration_cutoff;
    wide = non_str_templates & templateDuration_us > waveform_duration_cutoff;
    
    % Striatum classification
    prop_long_isi_cutoff = 0.35;
    cv2_cutoff = 0.8;
    
    msn = str_templates & ...
        templateDuration_us > waveform_duration_cutoff & ...
        prop_long_isi >= prop_long_isi_cutoff;
    
    fsi = str_templates & ...
        templateDuration_us <= waveform_duration_cutoff & ...
        prop_long_isi < prop_long_isi_cutoff;
    
    tan = str_templates & ...
        templateDuration_us > waveform_duration_cutoff & ...
        prop_long_isi < prop_long_isi_cutoff;
    
    uin = str_templates & ~msn & ~fsi & ~tan;
    
    waveform_t = 1e3*((0:size(templates,2)-1)/ephys_sample_rate);
    
    if verbose
        
        % Plot the waveforms and spike statistics
        figure;
        
        if any(non_str_templates)
            subplot(2,2,1); hold on;
            p = plot(waveform_t,waveforms(non_str_templates,:)');
            set(p(wide(non_str_templates)),'color','k')
            set(p(narrow(non_str_templates)),'color','r')
            xlabel('Time (ms)')
            title('Not striatum');
            legend([p(find(wide(non_str_templates),1)),p(find(narrow(non_str_templates),1))],{'Wide','Narrow'})
        end
        
        subplot(2,2,2); hold on;
        p = plot(waveform_t,waveforms(str_templates,:)');
        set(p(msn(str_templates)),'color','m')
        set(p(fsi(str_templates)),'color','b')
        set(p(tan(str_templates)),'color','g')
        set(p(uin(str_templates)),'color','c')
        xlabel('Time (ms)')
        title('Striatum');
        legend([p(find(msn(str_templates),1)),p(find(fsi(str_templates),1)), ...
            p(find(tan(str_templates),1)),p(find(uin(str_templates),1))],{'MSN','FSI','TAN','UIN'});
        
        subplot(2,2,3); hold on;
        
        stem3( ...
            templateDuration_us(wide)/1000, ...
            prop_long_isi(wide), ...
            spike_rate(wide),'k');
        
        stem3( ...
            templateDuration_us(narrow)/1000, ...
            prop_long_isi(narrow), ...
            spike_rate(narrow),'r');
        
        xlabel('waveform duration (ms)')
        ylabel('frac long ISI')
        zlabel('spike rate')
        
        set(gca,'YDir','reverse')
        set(gca,'XDir','reverse')
        view(3);
        grid on;
        axis vis3d;
        
        subplot(2,2,4); hold on;
        stem3( ...
            templateDuration_us(msn)/1000, ...
            prop_long_isi(msn), ...
            spike_rate(msn),'m');
        
        stem3( ...
            templateDuration_us(fsi)/1000, ...
            prop_long_isi(fsi), ...
            spike_rate(fsi),'b');
        
        stem3( ...
            templateDuration_us(tan)/1000, ...
            prop_long_isi(tan), ...
            spike_rate(tan),'g');
        
        stem3( ...
            templateDuration_us(uin)/1000, ...
            prop_long_isi(uin), ...
            spike_rate(uin),'c');
        
        xlabel('waveform duration (ms)')
        ylabel('frac long ISI')
        zlabel('spike rate')
        
        set(gca,'YDir','reverse')
        set(gca,'XDir','reverse')
        view(3);
        grid on;
        axis vis3d;
        
        % Plot cell type by depth
        celltype_labels = {'Wide','Narrow','MSN','FSI','TAN','UIN'};
        celltypes = wide.*1 + narrow.*2 + msn.*3 + fsi.*4 + tan.*5 + uin.*6;
        use_colors = ...
            [0,0,0;
            1,0,0;
            1,0,1;
            0,0,1;
            0,1,0;
            0,1,1];
        
        plot_celltypes = any([wide,narrow,msn,fsi,tan,uin],1);
        
        figure('Position',[94,122,230,820]);
        scatter(rand(size(templateDepths))-0.5,templateDepths,10,use_colors(celltypes,:),'filled');
        xlim([-1,1])
        set(gca,'XTick',[]);
        set(gca,'YDir','reverse');
        ylabel('Depth (\mum)');
        legend(celltype_labels(plot_celltypes));
        ylim([0,max(channel_positions(:,2))])
        
        drawnow;
        
    end
    
end


%% Finished
if verbose; disp('Finished loading experiment.'); end








