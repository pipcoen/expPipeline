function [sig] = sigsort(obj, spikeTimes, clu, eventTimes, timeWinPre, timeWinPost)
        %Assigning default values
        if strcmp(obj.blocks{1,1}.expDef, 'multiSpaceWorld') 
            fltblk = prc.filtStruct(obj.blocks{1,1}, obj.blocks{1,1}.responseMade~=0);
        else
            fltblk = obj.blocks;
        end
        if ~exist('spikeTimes', 'var'); spikeTimes = fltblk.ephSpikeTimes; end
        if ~exist('clu', 'var'); clu = fltblk.ephSpikeTemplates; end
        if ~exist('eventTimes', 'var'); eventTimes = fltblk.stimPeriodStart; end
        if ~exist('timeWinPre', 'var'); timePreLow = 0; timePreUp = 1; 
        elseif length(timeWinPre) ~= 1; timePreLow = timeWinPre(1,1); timePreUp = timeWinPre(end);
        else
        timePreLow = 0; timePreUp = timeWinPre;
        end

        if ~exist('timeWinPost', 'var'); timePostLow = 0; timePostUp = 1; 
        elseif length(timeWinPost) ~= 1; timePostLow = timeWinPost(1,1); timePostUp = timeWinPost(end);
        else
        timePostLow = 0; timePostUp = timeWinPost;
        end

        eventTimes = eventTimes(~isnan(eventTimes));
        %Selecting unique clusters
        cluList = unique(clu);
        sig = zeros(length(cluList),1);
        %Testing for significance
        for i = 1:length(cluList)
            sigSpikes = spikeTimes(clu == cluList(i));
            eventIndex = 1:length(eventTimes);
            allSpikes = arrayfun(@(x) sigSpikes(sigSpikes >= (eventTimes(x)-timePreUp) & sigSpikes <= (eventTimes(x)+timePostUp))-eventTimes(x), eventIndex, 'uni',0);
            preSpikes = cellfun(@(x) x(x <= 0 - timePreLow), allSpikes, 'uni',0);
            postSpikes = cellfun(@(x) x(x >= 0 - timePostLow), allSpikes, 'uni',0);
            preSpikeIdx = cellfun(@length, preSpikes)/(timePreUp-timePreLow); 
            postSpikeIdx = cellfun(@length, postSpikes)/(timePostUp-timePostLow);
            z = ttest(postSpikeIdx, preSpikeIdx, 'alpha', 0.05/sqrt(length(cluList)));
            if z == 1; sig(i) = cluList(i); end    
        end
        %Storing and sorting significant clusters
        sig = sort(sig(sig~=0));
        end

   
    
        
