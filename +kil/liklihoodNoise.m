function liklihoodNoise(ksDir)
    s = loadKSdir(ksDir,struct('excludeNoise',false));
    
    %Immediate set all clusters where numSpikes<=100 to noise
    tab=tabulate(s.clu);
    nIdx = tab(:,2)<=100;
    writePhyTSV(ksDir, 'group', tab(nIdx), repmat({'noise'},sum(nIdx),1));

    templateIDs = unique(s.spikeTemplates);
    for clu = 1:length(templateIDs)
        cluID = templateIDs(clu);
        
        template = double(squeeze(s.temps(cluID + 1,:,:))); %+1 because cluID is 0-indexed 
        marginal = mean(abs(template),1);
        marginal = marginal/max(marginal);
        
        maxMarginal = max(abs(template),[],1);
        maxMarginal = maxMarginal/max(maxMarginal);

        probabilityNoise(clu,1) = trapz(marginal)/length(marginal); %Proportion of AUC      
%         tempScalingAmp(clu,1) = median(s.tempScalingAmps(s.spikeTemplates==cluID));
    end

    writePhyTSV(ksDir, 'pNoise', templateIDs, probabilityNoise);
%     writePhyTSV(ksDir, 'Amps', templateIDs, tempScalingAmp);

end