classdef behaviorAnalysis
    %% behaviorAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    % that can be used to plot various aspects of animal behavior.
    % NOTE: This function is designed to operate on the output of convertExpFiles and most methods are specific to multisensoty spatial integration plots.
    
    % Inputs(default values)
    % subjects({'PC011';'PC012';'PC013';'PC016'})------A cell array of subjects collect parameter and block files for
    % expDate('last')----------------------------------A cell array of dates, one for all subjects or one for each subject
    % processingTag('none')----------------------------A tag to indicate whether you want to combine data across mice or retain individuals.
    
    properties (Access=public)
        subjects;                %Cell array of subject names
        expDate;                 %Recording dates to use for each subject
        blocks;                  %Block files loaded for each subject (one cell per subject)
        params;                  %Parameter files loaded for each subject (one cell per subject)
        uniqueConditions;        %The unique conditions presented to each subject (one cell per subject)
        audValues;               %List of auditory values (amplitude or azimuth changes) pesented (one cell per subject)
        visValues;               %List of visual values pesented (one cell per subject)
        audType;                 %Azimuthal changes or auditory changes (one cell per subject)
        currentAxes;             %Handle to current axis being used for plotting
    end
    
    %%
    methods
        function obj = behaviorAnalysis(subjects, expDate, processingTag)
            % Initialize fields with default values if no vaules are provided.
            if ~exist('subjects', 'var') || isempty(subjects); subjects = {'PC011';'PC012';'PC013';'PC015';'PC010';'PC017'}; end
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~exist('processingTag', 'var'); processingTag = 'none'; end
            if ~iscell(subjects); subjects = {subjects}; end
            if ~iscell(expDate); expDate = {expDate}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            obj = changeMouse(obj, subjects, expDate, processingTag);
        end
        
        function viewBoxPlots(obj, plotType)
            figure;
            if ~exist('plotType', 'var'); plotType = 'response'; end
            maxGrid = max(cellfun(@length, [obj.audValues obj.visValues]), [], 1);
            if  strcmpi(plotType(1:3), 'las')
                obj.subjects = {'Left Visual';'Both Visual';'Right Visual';'Left M2';'Both M2';'Right M2';'Left Control';'Both Control';'Right Control'};
                obj.blocks = repmat(obj.blocks, length(obj.subjects),1);
                obj.audValues = repmat(obj.audValues, length(obj.subjects),1);
                obj.visValues = repmat(obj.visValues, length(obj.subjects),1);
                obj.audType = repmat(obj.audType, length(obj.subjects),1);
            end
            
            for i  = 1:length(obj.subjects)
                colorYTick = [0 1];
                colorMap = redblue(64);
                colorDirection = 'normal';
                axisLimits = [0 1];
                reactionTimes = vertcat(obj.blocks{i}.reactionTime) < diff(vertcat(obj.blocks{i}.laserOnOff), [], 2);
                blkOff = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserSession)~=0 & vertcat(obj.blocks{i}.laserType)==0 & reactionTimes);
                blkNorm = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserSession)==0);
                
                if strcmpi(plotType(1:3), 'abs') || (~strcmpi(plotType(1:3), 'las') && length(blkNorm.response) > length(blkOff.response))
                    blkOff = blkNorm;
                end
                colorLabel = 'Fraction of right turns';
                switch lower(plotType(1:3))
                    case 'res'
                        blk = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==0);
                        plotData = makeGrid(blk, blk.response==2, @mean, 1);
                    case 'rea'
                        blk = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==0);
                        plotData = makeGrid(blk, blk.reactionTime, @median, 1)*1000;
                        axisLimits = [round(min(plotData(:))-5,-1) round(max(plotData(:))+5,-1)];
                        colorLabel = 'Median rection time';
                        colorYTick = {'Fast'; 'Slow'};
                        colorMap = flipud(colorMap);
                        colorDirection = 'reverse';
                    case 'ina'
                        blk = concatinateBlocks(obj.blocks{i},  vertcat(obj.blocks{i}.laserSession)==1);
                        plotData =  makeGrid(blk, blk.laserType>0, @mean, 1);
                        colorLabel = 'Fraction laser trials';
                    case {'las'}
                        axisLimits = [-1 1];
                        corticalLocations = vertcat(obj.blocks{i}.galvoPosition);
                        if contains(obj.subjects{i}, 'Left'); hemiMod = -1; else, hemiMod = 1; end
                        if contains(obj.subjects{i}, 'Both'); laserType = 2; corticalLocations(:,1) = abs(corticalLocations(:,1)); else, laserType = 1; end
                        if contains(obj.subjects{i}, 'Visual')
                            laserTrials = corticalLocations(:,1)==(2.5*hemiMod) & corticalLocations(:,2)==-3;
                        elseif contains(obj.subjects{i}, 'M2')
                            laserTrials = corticalLocations(:,1)==(1*hemiMod) & corticalLocations(:,2)==2;
                        elseif contains(obj.subjects{i}, 'Control')
                            laserTrials = (corticalLocations(:,1)==(1.5*hemiMod) | corticalLocations(:,1)==(2.5*hemiMod)) & corticalLocations(:,2)==0;
                        end
                        blkLas= concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==laserType & laserTrials & reactionTimes);
                        
                        if length(blkLas.trialEnd) < 50; continue; end
                        fracRightTurnsLaser = makeGrid(blkLas, blkLas.response==2, @mean, 1);
                        fracRightTurnsNorm = makeGrid(blkOff, blkOff.response==2, @mean, 1);
                        plotData = fracRightTurnsLaser-fracRightTurnsNorm;
                        trialEst = floor(length(blkLas.response)/8);
                        maxGrid = [1 1];
                end
                plotIdx = ~all(isnan(plotData));
                obj.getAxes(obj.subjects{i}, [60 30], maxGrid(2)/(1.3*maxGrid(1)), [50 50], [50 100]);
                imsc(plotData(:,plotIdx), axisLimits, colorMap, 'k');
                daspect([1 1 1]); axis xy;
                title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, trialEst, length(obj.blocks{i})));
                set(gca, 'xTick', 1:size(plotData(:,plotIdx),2), 'xTickLabel', obj.visValues{i}(plotIdx)*100, 'fontsize', 14)
                set(gca, 'yTick', 1:size(plotData(:,plotIdx),1), 'yTickLabel', obj.audValues{i}*10, 'fontsize', 14)
                ylabel(obj.audType{i});
                xlabel('Visual Contast');
                box off;
            end
            currentAxisPotision = get(gca, 'position');
            colorBar = colorbar;
            set(colorBar,'Ticks', get(colorBar, 'Limits'), 'TickLabels', colorYTick, 'YDir', colorDirection);
            set(gca, 'position', currentAxisPotision);
            figureSize = get(gcf, 'position');
            colorLabel = ylabel(colorBar, colorLabel);
            set(colorLabel, 'position', [1 mean(get(colorBar, 'Limits')) 0], 'FontSize', 14)
            set(colorBar, 'position', [1-75/figureSize(3), 0.2, 30/figureSize(3), 0.6])
        end
        
        function viewPsychoCurves(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            PF = @PAL_Logistic;
            searchGrid.alpha = -20:1:20;
            searchGrid.beta = 0:0.01:2;
            searchGrid.gamma = 0:.01:.1;
            searchGrid.lambda = 0:.01:.1;
            
            if  strcmpi(plotType(1:3), 'las')
                obj.subjects = {'Left Visual';'Both Visual';'Right Visual';'Left M2';'Both M2';'Right M2';'Left Control';'Both Control';'Right Control'};
                obj.blocks = repmat(obj.blocks, length(obj.subjects),1);
                obj.audValues = repmat(obj.audValues, length(obj.subjects),1);
            end
            
            for i  = 1:length(obj.subjects)
                reactionTimes = vertcat(obj.blocks{i}.reactionTime) < diff(vertcat(obj.blocks{i}.laserOnOff), [], 2);
                blkOff = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserSession)~=0 & vertcat(obj.blocks{i}.laserType)==0 & reactionTimes);
                blkNorm = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserSession)==0);
                
                if strcmpi(plotType(1:3), 'abs') || (~strcmpi(plotType(1:3), 'las') && length(blkNorm.response) > length(blkOff.response))
                    blkOff = blkNorm;
                end
                corticalLocations = vertcat(obj.blocks{i}.galvoPosition);
                obj.getAxes(obj.subjects{i}, [200 40], [], [60 80], [60 20]); hold on; box off;
                visGrid = blkOff.visGrid*100;
                colorOpt = 'bkr';
                xlim([min(visGrid(:))-5, max(visGrid(:))+5]);
                for j = 1:length(obj.audValues{i})
                    if j == 1 && ~strcmp(plotType(1:3), 'rea')
                        yL = ylim; hold on; plot([0 0],yL, '--k', 'linewidth', 1.5);
                        xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                    end
                    numTrials = makeGrid(blkOff, blkOff.feedback*0+1, @sum, 1);
                    trialType = makeGrid(blkOff, blkOff.trialType, @mean, 1);
                    tkIdx = numTrials>0 & blkOff.audGrid == obj.audValues{i}(j);
                    if sum(tkIdx(:)) < 3; continue; end
                    conflicts = tkIdx.*trialType==4;
                    nonConflict = (visGrid.*blkOff.audGrid)>=0 & (tkIdx);
                    useColor =colorOpt(sign(obj.audValues{i}(j))+2);
                    
                    switch lower(plotType(1:3))
                        case {'res'; 'las'; 'abs'}
                            StimLevelsFineGrain = min(visGrid(tkIdx)):(max(visGrid(tkIdx))/1000):max(visGrid(tkIdx));
                            fracRightTurns = makeGrid(blkOff, blkOff.response==2, @mean, 1);
                            numRightTurns = makeGrid(blkOff, blkOff.response==2, @sum, 1);
                            numTrials = makeGrid(blkOff, blkOff.response==2, @length, 1);
                            
                            paramsValues = PAL_PFML_Fit(visGrid(tkIdx)', numRightTurns(tkIdx)', numTrials(tkIdx)', ...
                                searchGrid, [1,1,1,1], PF,'lapseLimits',[0 1],'guessLimits', [0  1]);
                            
                            plot(StimLevelsFineGrain, PF(paramsValues, StimLevelsFineGrain),'LineWidth',1.1, 'color', useColor);
                            %%
                            if strcmpi(plotType(1:3), 'las')
                                if contains(obj.subjects{i}, 'Left'); hemiMod = -1; else, hemiMod = 1; end
                                if contains(obj.subjects{i}, 'Both'); laserType = 2; corticalLocations(:,1) = abs(corticalLocations(:,1)); else, laserType = 1; end
                                if contains(obj.subjects{i}, 'Visual')
                                    laserTrials = corticalLocations(:,1)==(2.5*hemiMod) & corticalLocations(:,2)==-3;
                                elseif contains(obj.subjects{i}, 'M2')
                                    laserTrials = corticalLocations(:,1)==(1*hemiMod) & corticalLocations(:,2)==2;
                                elseif contains(obj.subjects{i}, 'Control')
                                    laserTrials = (corticalLocations(:,1)==(1.5*hemiMod) | corticalLocations(:,1)==(2.5*hemiMod)) & corticalLocations(:,2)==0;
                                end
                                blkLas = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==laserType & laserTrials & reactionTimes);
                                if length(blkLas.trialEnd) < 50; continue; end
                                fracRightTurns = makeGrid(blkLas, blkLas.response==2, @mean, 1);
                                noNanIdx = (tkIdx.*~isnan(fracRightTurns))>0;
                                plot(visGrid((conflicts.*noNanIdx)>0),fracRightTurns((conflicts.*noNanIdx)>0),'^','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                plot(visGrid((nonConflict.*noNanIdx)>0),fracRightTurns((nonConflict.*noNanIdx)>0),'o','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                trialEst = floor(length(blkLas.response)/8);
                            else
                                plot(visGrid(conflicts),fracRightTurns(conflicts),'^','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                plot(visGrid(nonConflict),fracRightTurns(nonConflict),'o','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                trialEst = length(blkOff.response);
                            end
                            yAxLabel = 'Fraction of right turns';
                        case 'rea'
                            yAxLabel = 'Rection time (ms)';
                            reactionTimesGrid = makeGrid(blkOff, blkOff.reactionTime, @median, 1)*1000-500;
                            plot(visGrid(tkIdx), reactionTimesGrid(tkIdx),'LineWidth',1.1, 'color', useColor);
                            plot(visGrid(conflicts), reactionTimesGrid(conflicts),'^','markerfacecolor', useColor, 'markeredgecolor', useColor);
                            plot(visGrid(nonConflict), reactionTimesGrid(nonConflict),'o','markerfacecolor', useColor, 'markeredgecolor', useColor);
                            ylim([0 500])
                    end
                end
                title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, trialEst, blkOff.nSessions));
                set(gca, 'XTick', -45:15:45);
            end
            figureSize = get(gcf, 'position');
            mainAxes = [60./figureSize(3:4) 1-2*(50./figureSize(3:4))];
            suplabel('\fontsize{16} Auditory on {\color{red}right, \color{black}centre, or \color{blue}left}', 't', mainAxes);
            suplabel(['\fontsize{16} ' yAxLabel], 'y', mainAxes);
            suplabel('\fontsize{16} Visual Contrast', 'x', mainAxes);
        end
        
        function viewJitterPlot(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'per'; end
            chosenSubjects = 1:length(obj.subjects);
            switch lower(plotType(1:3))
                case 'amp'
                    maxPlotLength = max(ceil(cellfun(@length, obj.audValues)/2));
                    mainTitle = '\fontsize{16} Effect of auditory amplitude on performance';
                    mainYLabel = '\fontsize{16} Fraction correct';
                    mainXLabel = '\fontsize{16} Audio amplitude';
                    sigPlot = 0;
                case 'rig'
                    if length(plotType)==3; plotType = [plotType 'aud']; end
                    eval(['sufficientConditions = cellfun(@(x) sum(unique(abs(x))>0), obj.' lower(plotType(4:end)) 'Values);']);
                    chosenSubjects = find((cellfun(@(x) length(unique({x(:).rigName}')), obj.blocks)>1).*sufficientConditions>0)';
                    maxPlotLength = max(cellfun(@(x) length(unique({x(:).rigName}')), obj.blocks));
                    mainTitle = ['\fontsize{16} Effect of rig on ' plotType(4:end) '  performance'];
                    mainYLabel = '\fontsize{16} Fraction correct';
                    mainXLabel = '\fontsize{16} Rig Name';
                case 'mul'
                    mainTitle = '\fontsize{16} Percentage correct by condition';
                    mainYLabel = '\fontsize{16} Fraction correct';
                    mainXLabel = '\fontsize{16} Condition';
                    sigPlot = 0;
            end
            
            figure;
            for i  = chosenSubjects
                switch lower(plotType(1:3))
                    case 'amp'
                        blk = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==0);
                        percentCorrect = makeGrid(blk, blk.feedback>0, @mean, 2);
                        xData = blk.audGrid(blk.audGrid>=0 & blk.visGrid==0);
                        yData = percentCorrect(blk.audGrid>=0 & blk.visGrid==0);
                        [~, pVal] = cellfun(@(x) ttest(x-0.5), yData);
                        plotColors = repmat('k', length(pVal),1);
                        plotColors(pVal<0.05) = 'k';
                        
                    case 'rig'
                        blk = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==0);
                        [rigsUsed,~,rigIdx] = unique(blk.rigName);
                        for j = 1:length(rigsUsed)
                            percentCorrect = makeGrid(blk, blk.feedback>0, @mean, 2, find(rigIdx==j));
                            if strcmpi(plotType(4:end), 'vis')
                                visTest = obj.visValues{i}(find(obj.visValues{i}>0, 3));
                                yData{j,1} = cell2mat(percentCorrect(blk.audGrid==0 & blk.visGrid==visTest(end)));
                            else
                                audTest = obj.audValues{i}(find(obj.audValues{i}>0, 1));
                                yData{j,1} = cell2mat(percentCorrect(blk.audGrid==audTest & blk.visGrid==0));
                            end
                        end
                        xData = rigsUsed;
                        ttestCombos = num2cell(nchoosek(1:length(rigsUsed),2),2);
                        [~, significance] = cellfun(@(x) ttest2(yData{x(1)}, yData{x(2)}), ttestCombos);
                        if any(significance<0.05); sigPlot = 1; else, sigPlot = 0; end
                        testedPairs = ttestCombos(significance<0.05);
                        significance = significance(significance<0.05);
                        plotColors = 'kr'; plotColors = plotColors(:);
                        
                    case 'mul'
                        blk = concatinateBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==0);
                        trialType = makeGrid(blk, blk.trialType, @mean, 3);
                        multiTrials = trialType==3;
                        percentCorrect = makeGrid(blk, blk.feedback>0, @mean, 2);
                        numTrials = makeGrid(blk, blk.feedback*0+1, @sum, 3);
                        if ~any(multiTrials); continue; end
                        
                        idx = 0;
                        for j = find(multiTrials)'
                            yData(idx+1,1) = percentCorrect(blk.audGrid==blk.audGrid(j) & trialType==1);
                            yData(idx+2,1) = percentCorrect(blk.visGrid==blk.visGrid(j) & trialType==2);
                            yData(idx+3,1) = percentCorrect(j);
                            xData(idx+1:idx+3,1) = {'Aud'; 'Vis'; 'Mul'};
                            xDataPosition(idx+1:idx+3,1) = floor((idx+1)/3)+(idx+1:idx+3);
                            idx = idx+3;
                        end
                        [~, pVal] = cellfun(@(x) ttest(x-0.5), yData);
                        plotColors = repmat('k', length(pVal),1);
                        plotColors(pVal<0.05) = 'k';
                        yLimits = [0.4 1];
                end
                if ~exist('xDataPosition', 'var'); xDataPosition = 1:length(xData); end
                if ~exist('maxPlotLength', 'var'); maxPlotLength = length(xData)/2; end
                if ~exist('yLimits', 'var'); yLimits = [0 1]; end
                obj.getAxes(obj.subjects{i}, [80 70], maxPlotLength/5, [70 70], [60 20], obj.subjects(chosenSubjects)); hold on; box off;
                jitterPlot(yData, 'edgC', plotColors, 'xPos', xDataPosition, 'curA', 1); grid('on');
                
                set(gca, 'xTickLabels', xData);
                title(sprintf('%s: n = %d', obj.subjects{i}, length(percentCorrect{1})));
                ylim(yLimits);
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                
                if sigPlot
                    sigstar(testedPairs, significance);
                end
                
            end
            figureSize = get(gcf, 'position');
            mainAxes = [60./figureSize(3:4) 1-2*(50./figureSize(3:4))];
            suplabel(mainTitle, 't', mainAxes);
            suplabel(mainYLabel, 'y', mainAxes);
            suplabel(mainXLabel, 'x', mainAxes);
        end
        
        function viewInactivationEffects(obj, UniOrBilateral, individualFlies)
            imgBWOuline=imread('BrainOutlineBW.png');
            if ~exist('UniOrBilateral', 'var'); UniOrBilateral = 2; end
            if ~exist('individualFlies', 'var'); individualFlies = 1; end
            
            if UniOrBilateral==1; galvoType = 1; else, galvoType = 2; end
            
            switch individualFlies
                case 0
                    laserBlocks = cell(length(obj.subjects),1);
                    for i  = 1:length(obj.subjects)
                        laserBlocks{i,1} = obj.blocks{i}(arrayfun(@(x) any(all([x.laserType>0 x.galvoType==galvoType],2)), obj.blocks{i}));
                    end
                    blkOff = concatinateBlocks(cat(1,laserBlocks{:}), {vertcat(obj.blocks{i}.laserType)==0});
                    blkOn = concatinateBlocks(cat(1,laserBlocks{:}), 'onlyLaser');
                    
                    if UniOrBilateral==2
                        blkOff.galvoPosition(:,1) = abs(blkOff.galvoPosition(:,1));
                        blkOn.galvoPosition(:,1) = abs(blkOn.galvoPosition(:,1));
                    end
                    galvoGrid = unique(blkOn.galvoPosition, 'rows');
                    blkOn.response(blkOn.reactionTime>diff(blkOn.laserOnOff, [], 2)) = 0;
                    blkOff.response(blkOff.reactionTime>diff(blkOff.laserOnOff, [], 2)) = 0;
                    
                    idx = 0;
                    yOpt = {'VisLeft'; 'VisRight'; 'ConflictVisLeft'; 'ConflictVisRight'};
                    tOpt = {'Timeout'; 'ChooseLeft'; 'ChooseRight'};
                    for conditions = [unique(blkOn.conditions(blkOn.trialType == 2))' (unique(blkOn.conditions(blkOn.trialType == 4))*-1)']
                        idx = idx+1;
                        idx2 = 0;
                        for response = 0:2
                            idx2 = idx2+1;
                            baseline = mean(blkOff.response(blkOff.conditions==conditions)==response);
                            siteEffects = arrayfun(@(x,y) mean(blkOn.response(all(blkOn.galvoPosition==[x, y],2) & blkOn.conditions==conditions)==response), galvoGrid(:,1), galvoGrid(:,2));
                            subplot(4,3,(idx-1)*3 + idx2);
                            siteEffects = siteEffects-baseline;
                            imagesc(linspace(-4.5,4.5,1000),linspace(3.75,-5.2,1000),imgBWOuline); axis xy;
                            hold on;
                            scatter(galvoGrid(:,1),galvoGrid(:,2),150,siteEffects,'o','filled'); axis equal;  drawnow;
                            scatter(galvoGrid(:,1)*-1,galvoGrid(:,2),150,siteEffects,'o','filled'); axis equal;  drawnow;
                            colormap(redblue); caxis([-0.8 0.8]);
                            box off;
                            if idx ==1; title(tOpt{idx2}); end
                            if idx2 ==1; ylabel(yOpt{idx}); end
                        end
                    end
                    
                    
                    
                case 1
                    for i  = 1:length(obj.subjects)
                        laserBlocks = obj.blocks{i}(arrayfun(@(x) any(all([x.laserType>0 x.galvoType==galvoType],2)), obj.blocks{i}));
                        blk = concatinateBlocks(laserBlocks, 'none');
                        laserDuration = diff(blk.laserOnOff,[],2);
                        visonCorrect = sign(blk.visInitialAzimuth)>0;
                        noLaserConditions = arrayfun(@(x) ~any(blk.laserType(blk.conditions==x)), blk.uniqueConditionsIdx);
                        tkIdx = (blk.reactionTime<laserDuration);
                        tkIdx = tkIdx.*(ismember(blk.conditions, blk.uniqueConditionsIdx(~noLaserConditions)));
                        
                        if strcmp(UniOrBilateral(1:3), 'bil'); blk.galvoPosition(:,1) = abs(blk.galvoPosition(:,1)); end
                        gGrid = unique(blk.galvoPosition, 'rows');
                        trialLabel = {'auditory performance'; 'visual performance'; 'multisensory performance'; 'conflict visual choice'};
                        %%
                        for j = 1:4
                            figure(figHand{j}); obj.getAxes(obj.subjects{i}, [80 0], [], [70 70], [60 80]); hold on;
                            imagesc(linspace(-4.5,4.5,1000),linspace(3.75,-5.2,1000),imgBWOuline); axis xy;
                            box off;
                            
                            gPos = blk.galvoPosition(tkIdx>0 & blk.trialType == j,:);
                            dIdx = blk.sessionIdx(tkIdx>0 & blk.trialType == j);
                            lTyp = blk.laserType(tkIdx>0 & blk.trialType == j);
                            resp = blk.response(tkIdx>0 & blk.trialType == j);
                            visR = visonCorrect(tkIdx>0 & blk.trialType == j);
                            expPerform = dIdx*0;
                            
                            if j < 4; fBck = blk.feedback(tkIdx>0 & blk.trialType == j);
                            else, fBck = resp==2 == visR;
                            end
                            
                            for k = unique(blk.sessionIdx)'
                                expPerform(dIdx==k) = mean(fBck(dIdx==k & lTyp==0));
                            end
                            
                            
                            numberOfLaserTrials = arrayfun(@(x,y) sum(all(gPos==[x, y],2) & lTyp>0), gGrid(:,1), gGrid(:,2));
                            diffLaserOn = 100*arrayfun(@(x,y) mean(fBck(all(gPos==[x, y],2) & lTyp>0) - expPerform(all(gPos==[x, y],2) & lTyp>0)), gGrid(:,1), gGrid(:,2));
                            scatter(gGrid(:,1),gGrid(:,2),150,diffLaserOn,'o','filled'); axis equal;  drawnow;
                            colormap(redblue); caxis([-35 35]);
                            axis off;
                            title(sprintf('%s: %.0f Trials/Site', obj.subjects{i}, mean(numberOfLaserTrials))); axis off;
                            figureSize = get(gcf, 'position');
                            mainAxes = [60./figureSize(3:4) 1-2*(50./figureSize(3:4))];
                            suplabel(['\fontsize{16} Effects on ' trialLabel{j}], 't', mainAxes);
                        end
                    end
                    for j = 1:4
                        %%
                        figure(figHand{j});
                        currentAxisPotision = get(gca, 'position');
                        colorBar = colorbar;
                        set(colorBar,'Limits', [-35 35], 'Ticks', [-35 35]);
                        set(gca, 'position', currentAxisPotision);
                        figureSize = get(gcf, 'position');
                        colorLabel = ylabel(colorBar, 'Change with inactivation');
                        set(colorLabel, 'position', [1 mean(get(colorBar, 'Limits')) 0], 'FontSize', 14)
                        set(colorBar, 'position', [1-100/figureSize(3), 0.2, 30/figureSize(3), 0.6])
                    end
            end
        end
        
        function obj = getAxes(obj, subject, axisGap, figureRatio, botTopEdge, leftRightEdge, allSubjects)
            if ~exist('subject', 'var') || isempty(subject); subject = obj.subject{1}; end
            if ~iscell(subject); subject = {subject}; end
            if ~exist('figureRatio', 'var') || isempty(figureRatio); figureRatio = 1; end
            if ~exist('axisGap', 'var') || isempty(axisGap); axisGap = 25; end
            if ~exist('botTopEdge', 'var') || isempty(botTopEdge); botTopEdge = 50; end
            if ~exist('leftRightEdge', 'var') || isempty(leftRightEdge); leftRightEdge = 50; end
            if ~exist('allSubjects', 'var') || isempty(allSubjects); allSubjects = obj.subjects; end
            
            screenSize = get(0,'MonitorPositions');
            screenSize = screenSize(screenSize(:,1)==min(screenSize(:,1)),:);
            screenRatio = round(screenSize(3)/screenSize(4));
            numOfSubjects = length(allSubjects);
            
            if numOfSubjects < 4; numOfRows = 1;
            else, numOfRows = find(((1:5)*screenRatio.*(1:5))>=numOfSubjects,1);
            end
            numOfCols = ceil(numOfSubjects/numOfRows);
            figureSize = min([400*numOfCols*figureRatio, 400*numOfRows], screenSize(3:4));
            
            botTopEdge = botTopEdge/figureSize(2);
            leftRightEdge = leftRightEdge/figureSize(1);
            axisGap = axisGap./figureSize;
            
            obj.currentAxes = tightSubplot(numOfRows, numOfCols, find(contains(allSubjects,subject)), axisGap, botTopEdge, leftRightEdge);
            set(gcf, 'position', [screenSize(1:2)+screenSize(3:4)-figureSize-[0 75], figureSize])
        end
        
        function obj = changeMouse(obj, subjects, expDate, processingTag)
            %Get block and parameter files for the requested dates.
            [obj.blocks, obj.params]  = cellfun(@(x) getFilesFromDates(x, expDate, 'bloprm'), subjects, 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
            
            retainIdx = ones(length(obj.params),1)>0;
            switch lower(processingTag)
                case 'none'
                    %If there are multiple parameter sets in the requested date range for a mouse, use only the most common parameter set.
                    for i = 1:length(subjects)
                        mouseIdx = strcmp({obj.blocks.subject}', subjects{i});
                        [conditionSets, ~, setIdx] = unique(cellfun(@(x) num2str(x(:)'),{obj.blocks(mouseIdx).uniqueConditions}','uni',0));
                        if length(conditionSets)>1
                            fprintf('WARNING: Several parameter sets in date range for %s. Using mode\n', subjects{i});
                            retainIdx(mouseIdx) = retainIdx(mouseIdx).*(setIdx == mode(setIdx));
                        end
                    end
                    obj.blocks = obj.blocks(retainIdx);
                    obj.params = obj.params(retainIdx);
                    obj.subjects = unique({obj.blocks.subject}');
                    
                    [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
                    obj.blocks = arrayfun(@(x) obj.blocks(subjectIdx==x), 1:length(subjects), 'uni', 0)';
                    obj.params = arrayfun(@(x) obj.params(subjectIdx==x), 1:length(subjects), 'uni', 0)';
                case 'combine'
                    %Remove files without enough trials.
                    %Remove files with performance lower than 60%
            end
            
            for i = 1:length(obj.subjects)
                blocks = concatinateBlocks(obj.blocks{i});
                uniqueConditions = double(blocks.uniqueConditions);
                if size(unique(uniqueConditions, 'rows'),1) < size(blocks.uniqueConditions,1)
                    error('Unprepared to multisensory combinations of this nature');
                end
                if length(unique(uniqueConditions(:,1)))>1 && length(unique(abs(uniqueConditions(~isinf(uniqueConditions(:,3)),3))))>1
                    error('Detected changes in both audAmplitude and audInitialAzimuth so cannot plot');
                end
                
                if length(unique(uniqueConditions(:,1)))==1
                    obj.audType{i,1} = 'Auditory azimuth';
                    uniqueConditions = [uniqueConditions(:,3) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];
                else
                    obj.audType{i,1} = 'Auditory amplitude';
                    uniqueConditions = [uniqueConditions(:,1).*sign(uniqueConditions(:,3)) uniqueConditions(:,2).*sign(uniqueConditions(:,4))];
                end
                obj.audValues{i,1} = unique(uniqueConditions(:,1));
                obj.visValues{i,1} = unique(uniqueConditions(:,2));
                [visGridConditions, audGridConditions] = meshgrid(obj.visValues{i}, obj.audValues{i});
                [~, gridIdx] = ismember(uniqueConditions, [audGridConditions(:) visGridConditions(:)], 'rows');
                conditionsInGrid = nan*ones(length(obj.audValues{i}), length(obj.visValues{i}));
                conditionsInGrid(gridIdx) = obj.blocks{i}(1).uniqueConditionsIdx;
                obj.uniqueConditions{i,1} = uniqueConditions;
                
                [obj.blocks{i}.conditionsInGrid] = deal(conditionsInGrid);
                [obj.blocks{i}.visGrid] = deal(visGridConditions);
                [obj.blocks{i}.audGrid] = deal(audGridConditions);
            end
            obj.currentAxes = [];
        end
    end
end