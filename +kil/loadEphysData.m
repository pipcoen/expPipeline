function loadEphysData(x)
%% Load timeline and associated inputs
if ~exist(x.rawTimeline, 'file'); error('No timeline file exists for requested ephys session'); end
fprintf('Loading timeline... \n');
timeline = load(x.rawTimeline); timeline = timeline.Timeline;

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
flipperFlipTimesTimeline = timeline.rawDAQTimestamps(flipperFlips)';

%% Load task/behavior
fprintf('Loading block file... \n');
block  = load(x.processedData, 'blk'); block  = block.blk;

%% Load ephys data (single long recording)
fprintf('Loading ephys... \n');
% These are the digital channels going into the FPGA
acqLiveSyncIdx = 2;
flipperSyncIdx = 4;
%%
% Load clusters and sync/photodiode
clusterGroups = tdfread([x.kilosortOutput '\cluster_group.tsv']);
sync = load([x.kilosortOutput '\sync.mat']); sync = sync.sync;

%%
% Read header information
headerFID = fopen([x.kilosortOutput '\dat_params.txt']);
headerInfo = textscan(headerFID,'%s %s', 'delimiter',{' = '});
fclose(headerFID);
headerInfo = [headerInfo{1}'; headerInfo{2}'];
header = struct(headerInfo{:});
%%
% Load spike data
ephysSampleRate = header.apSampleRate;
spikeTimes = double(readNPY([x.kilosortOutput '\spike_times.npy']))./ephysSampleRate;
spikeTemplates = readNPY([x.kilosortOutput '\spike_templates.npy']);
templates = readNPY([x.kilosortOutput '\templates.npy']);
channelPositions = readNPY([x.kilosortOutput '\channel_positions.npy']);
channelMap = readNPY([x.kilosortOutput '\channel_map.npy']); %#ok<NASGU>
winv = readNPY([x.kilosortOutput '\whitening_mat_inv.npy']);
templateAmplitudes = readNPY([x.kilosortOutput '\amplitudes.npy']);
%%
% Default channel map/positions are from end: make from surface
channelPositions(:,2) = max(channelPositions(:,2)) - channelPositions(:,2);

% Get experiment index by finding numbered folders
sessionNum = str2double(x.sessionNum);

if exist('flipperFlipTimesTimeline','var')    
    % Get flipper experiment differences by long delays
    flipThresh = 1; % time between flips to define experiment gap (s)
    flipTimes = sync(flipperSyncIdx).timestamps;
    flipperExptIdx = [1;find(diff(flipTimes) > flipThresh)+1;length(flipTimes)+1];
    
    flipperFlipTimesEphys = flipTimes(flipperExptIdx(sessionNum):flipperExptIdx(sessionNum+1)-1);
    
    % Check that number of flipper flips in timeline matches ephys
    if length(flipperFlipTimesEphys) ~= length(flipperFlipTimesTimeline)
        warning([animal ' ' day ':Flipper flip times different in timeline/ephys'])
        bad_flipper = true;
    end
    
    %         % Plot stim aligned by flipper times for sanity check
    %         figure; hold on;
    %         stim_ephys_timestamps = interp1(flipperFlipTimesEphys,flipperFlipTimesTimeline,sync(1).timestamps,'linear','extrap');
    %         plot(stim_ephys_timestamps,sync(1).values*5,'.g')
    %         plot(stimOn_times,5,'ob');
    
    sync_timeline = flipperFlipTimesTimeline;
    sync_ephys = flipperFlipTimesEphys;
end

if ~exist('flipperFlipTimesTimeline','var') || bad_flipper
    % (if no flipper or flipper problem, use acqLive)
    
    % Get acqLive times for current experiment
    experiment_ephys_starts = sync(acqLiveSyncIdx).timestamps(sync(acqLiveSyncIdx).values == 1);
    experiment_ephys_stops = sync(acqLiveSyncIdx).timestamps(sync(acqLiveSyncIdx).values == 0);
    acqlive_ephys_currexpt = [experiment_ephys_starts(sessionNum), ...
        experiment_ephys_stops(sessionNum)];
    
    sync_timeline = acqLiveSrtEnd;
    sync_ephys = acqlive_ephys_currexpt;
    
    % Check that the experiment time is the same within threshold
    % (it should be almost exactly the same)
    if abs(diff(acqLiveSrtEnd) - diff(acqlive_ephys_currexpt)) > 1
        error([animal ' ' day ': acqLive duration different in timeline and ephys']);
    end
end

% Get the spike/lfp times in timeline time (accounts for clock drifts)
spikeTimes_timeline = interp1(sync_ephys,sync_timeline,spikeTimes,'linear','extrap');
if load_lfp && exist(lfp_filename,'file')
    lfp_t_timeline = interp1(sync_ephys,sync_timeline,lfp_t,'linear','extrap');
end

% Get the depths of each template
% (by COM - this used to not work but now looks ok)
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(templates,winv,channelPositions(:,2),spikeTemplates,templateAmplitudes);
%     % (by max waveform channel)
%     template_abs = permute(max(abs(templates),[],2),[3,1,2]);
%     [~,max_channel_idx] =  max(template_abs,[],1);
%     templateDepths = channelPositions(max_channel_idx,2);
%     % Get each spike's depth
%     spikeDepths = templateDepths(spikeTemplates+1);

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
templateDuration_us = (templateDuration/ephysSampleRate)*1e6;

% Eliminate spikes that were classified as not "good"
if exist('clusterGroups','var')
    
    if verbose; disp('Removing non-good/MUA templates'); end
    
    good_templates_idx = uint32(clusterGroups{1}( ...
        strcmp(clusterGroups{2},'good') | strcmp(clusterGroups{2},'mua')));
    good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
    
    % Throw out all non-good template data
    templates = templates(good_templates,:,:);
    templateDepths = templateDepths(good_templates);
    waveforms = waveforms(good_templates,:);
    templateDuration = templateDuration(good_templates);
    templateDuration_us = templateDuration_us(good_templates);
    
    % Throw out all non-good spike data
    good_spike_idx = ismember(spikeTemplates,good_templates_idx);
    spikeTimes = spikeTimes(good_spike_idx);
    spikeTemplates = spikeTemplates(good_spike_idx);
    templateAmplitudes = templateAmplitudes(good_spike_idx);
    spikeDepths = spikeDepths(good_spike_idx);
    spikeTimes_timeline = spikeTimes_timeline(good_spike_idx);
    
    % Re-name the spike templates according to the remaining templates
    % (and make 1-indexed from 0-indexed)
    new_spike_idx = nan(max(spikeTemplates)+1,1);
    new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
    spikeTemplates = new_spike_idx(spikeTemplates+1);
    
elseif ~exist('clusterGroups','var')
    if verbose; disp('Clusters not yet sorted'); end
end


%% Classify spikes

if ephys_exists && load_parts.ephys && exist('clusterGroups','var')
    if verbose; disp('Classifying spikes...'); end
    
    str_templates = templateDepths >= str_depth(1) & templateDepths <= str_depth(2);
    non_str_templates = ~str_templates;
    
    % Define the window to look for spiking statistics in (spikes go in and
    % out, so take the bin with the largest firing rate for each cell and work
    % with that one)
    % spiking_stat_window = 60*5; % seconds
    % spiking_stat_bins = min(spikeTimes_timeline):spiking_stat_window: ...
    %     max(spikeTimes_timeline);
    
    % % (for whole session)
    spiking_stat_window = max(spikeTimes_timeline)-min(spikeTimes_timeline);
    spiking_stat_bins = [min(spikeTimes_timeline),max(spikeTimes_timeline)];
    
    % Get firing rate across the session
    bin_spikes = nan(max(spikeTemplates), ...
        length(spiking_stat_bins)-1);
    for curr_template = unique(spikeTemplates)'
        bin_spikes(curr_template,:) = ...
            histcounts(spikeTimes_timeline(spikeTemplates == curr_template), ...
            spiking_stat_bins);
    end
    min_spikes = 10;
    use_spiking_stat_bins = bsxfun(@ge,bin_spikes,prctile(bin_spikes,80,2)) & bin_spikes > min_spikes;
    spike_rate = sum(bin_spikes.*use_spiking_stat_bins,2)./ ...
        (sum(use_spiking_stat_bins,2)*spiking_stat_window);
    
    % Get proportion of ISI > 2s (Yamin/Cohen 2013) and CV2 (Stalnaker/Schoenbaum 2016)
    prop_long_isi = nan(max(spikeTemplates),1);
    cv2 = nan(max(spikeTemplates),1);
    for curr_template = unique(spikeTemplates)'
        
        long_isi_total = 0;
        isi_ratios = [];
        for curr_bin = find(use_spiking_stat_bins(curr_template,:))
            curr_spikeTimes = spikeTimes_timeline( ...
                spikeTimes_timeline > spiking_stat_bins(curr_bin) & ...
                spikeTimes_timeline < spiking_stat_bins(curr_bin+1) & ...
                spikeTemplates == curr_template);
            curr_isi = diff(curr_spikeTimes);
            
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
    
    waveform_t = 1e3*((0:size(templates,2)-1)/ephysSampleRate);
    
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
        ylim([0,max(channelPositions(:,2))])
        
        drawnow;
        
    end
    
end


%% Finished
if verbose; disp('Finished loading experiment.'); end








