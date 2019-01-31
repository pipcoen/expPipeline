classdef spatialAnalysis
    %% spatialAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    % that can be used to plot various aspects of animal behavior. NOTE: This function is designed to operate on the output of convertExpFiles and
    % most methods are specific to multisensoty spatial integration plots.
    %
    % Inputs(default values) subjects({'PC011';'PC012';'PC013';'PC016'})------A cell array of subjects collect parameter and block files for
    % expDate('last')----------------------------------A cell array of dates, one for all subjects or one for each subject
    % processingTag('none')----------------------------A tag to indicate whether you want to combine data across mice or retain individuals.
    
    properties (Access=public)
        subjects;                %Cell array of subject names
        expDate;                 %Recording dates to use for each subject
        blocks;                  %Block files loaded for each subject (one cell per subject)
        params;                  %Parameter files loaded for each subject (one cell per subject)
        splits;                  %Block splits blocks that are commonly used
        currentAxes;             %Handle to current axis being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate, processingTag)
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
            maxGrid = max(cell2mat(cellfun(@(x) [length(x.audValues) length(x.visValues)], obj.blocks, 'uni', 0)), [], 1);
            colorYTick = [0 1]; colorMap = redblue(64); colorDirection = 'normal'; axisLimits = [0 1]; colorLabel = 'Fraction of right turns';
            
            for i  = 1:length(obj.subjects)
                allReactionTimes = obj.blocks{i}.reactionTime < diff(obj.blocks{i}.laserOnOff, [], 2);
                if contains(plotType(1:3), {'las'; 'ina'}) || mode(obj.blocks{i}.laserSession>0)==1
                    laserOff = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserSession~=0 & obj.blocks{i}.laserType==0);
                else, laserOff = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserSession==0);
                end
                
                switch lower(plotType(1:3))
                    case 'res'
                        plotData = prc.makeGrid(laserOff, laserOff.response==2, @mean, 1);
                        trialNumber = length(laserOff.response);
                    case 'rea'
                        plotData = prc.makeGrid(laserOff, laserOff.reactionTime, @median, 1)*1000;
                        axisLimits = [round(min(plotData(:))-5,-1) round(max(plotData(:))+5,-1)];
                        colorLabel = 'Median rection time';
                        colorYTick = {'Fast'; 'Slow'};
                        colorMap = flipud(colorMap);
                        colorDirection = 'reverse';
                    case 'las'
                        blk = prc.combineBlocks(obj.blocks{i},  obj.blocks{i}.laserSession==1);
                        plotData =  prc.makeGrid(blk, blk.laserType>0, @mean, 1);
                        colorLabel = 'Fraction laser trials';
                    case 'ina'
                        if ~isempty(obj.expDate)
                            close;
                            for j  = 1:length(obj.subjects)
                                areaNames = spatialAnalysis.getSpecificGalvoParadigm(plotType(1:4));
                                newObj = createMouseReplicate(obj, j, areaNames);
                                newObj.expDate = [];
                                newObj.viewBoxPlots(plotType);
                            end
                            return;
                        end
                        [~, galvoSites, laserType] = spatialAnalysis.getSpecificGalvoParadigm(plotType(1:4));
                        selectedGalvoPositions = ismember(obj.blocks{i}.galvoPosition, galvoSites{i}, 'rows');
                        if length(obj.blocks{i}.trialEnd) < 50; continue; end
                        if ~contains(plotType, 'tim')
                            laserOn = prc.combineBlocks(obj.blocks{i}, selectedGalvoPositions & allReactionTimes & obj.blocks{i}.laserType == laserType(i));
                            plotData = prc.makeGrid(laserOn, laserOn.response==2, @mean, 1) - prc.makeGrid(laserOff, laserOff.response==2, @mean, 1);
                            trialNumber = length(laserOn.response);
                            colorLabel = 'Change in fraction of right turns';
                            colorYTick = [-0.75 0.75]; axisLimits = [-0.75 0.75];
                        else
                            laserOn = prc.combineBlocks(obj.blocks{i}, selectedGalvoPositions & obj.blocks{i}.laserType == laserType(i));
                            laserOnTimeouts = laserOn.reactionTime > diff(laserOn.laserOnOff, [], 2);
                            laserOff = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserSession~=0 & obj.blocks{i}.laserType==0);
                            laserOffTimeouts = laserOff.reactionTime > diff(laserOff.laserOnOff, [], 2);
                            plotData = prc.makeGrid(laserOn, laserOnTimeouts, @mean, 1) - prc.makeGrid(laserOff, laserOffTimeouts, @mean, 1);
                            trialNumber = length(laserOn.response);
                            colorLabel = 'Change in fraction of timeouts';
                            colorYTick = [-0.75 0.75]; axisLimits = [-0.75 0.75];
                        end
                end
                plotIdx = ~all(isnan(plotData));
                obj.plt.getAxes(obj.subjects{i}, [60 60], maxGrid(2)/(1.3*maxGrid(1)), [50 50], [75 100]);
                imsc(plotData(:,plotIdx), axisLimits, colorMap, 'k');
                daspect([1 1 1]); axis xy;
                title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, trialNumber, obj.blocks{i}.nSessions));
                set(gca, 'xTick', 1:size(plotData(:,plotIdx),2), 'xTickLabel', obj.blocks{i}.visValues(plotIdx)*100, 'fontsize', 14)
                set(gca, 'yTick', 1:size(plotData(:,plotIdx),1), 'yTickLabel', obj.blocks{i}.audValues*10, 'fontsize', 14)
                ylabel(obj.blocks{i}.audType);
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
            for i  = 1:length(obj.subjects)
                audBlocks = obj.splits.audLvls{i};
                colorOpt = 'bkr';
                switch plotType(1:3)
                    
                    case 'res'
                        obj.plt.getAxes(obj.subjects{i}, [200 40], [], [60 80], [60 20]); hold on; box off;
                        
                        
                    otherwise
                    
                        
                        xlim([min(visGrid(:))-5, max(visGrid(:))+5]);
                        for j = 1:length(obj.blocks{i}.audValues)
                            if j == 1 && ~strcmp(plotType(1:3), 'rea')
                                yL = ylim; hold on; plot([0 0],yL, '--k', 'linewidth', 1.5);
                                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                            end
                            numTrials = prc.makeGrid(laserOff, laserOff.feedback*0+1, @sum, 1);
                            trialType = prc.makeGrid(laserOff, laserOff.trialType, @mean, 1);
                            tkIdx = numTrials>0 & laserOff.grids.audValues == obj.blocks{i}.audValues(j);
                            if sum(tkIdx(:)) < 3; continue; end
                            conflicts = tkIdx.*trialType==4;
                            nonConflict = (visGrid.*laserOff.grids.audValues)>=0 & (tkIdx);
                            useColor =colorOpt(sign(obj.blocks{i}.audValues(j))+2);
                            
                            switch lower(plotType(1:3))
                                case {'res'; 'las'; 'abs'}
                                    StimLevelsFineGrain = min(visGrid(tkIdx)):(max(visGrid(tkIdx))/1000):max(visGrid(tkIdx));
                                    fracRightTurns = prc.makeGrid(laserOff, laserOff.response==2, @mean, 1);
                                    numRightTurns = prc.makeGrid(laserOff, laserOff.response==2, @sum, 1);
                                    numTrials = prc.makeGrid(laserOff, laserOff.response==2, @length, 1);
                                    
                                    [paramsValues, fittingFunction] = fitPsychoCurve(visGrid(tkIdx)', numRightTurns(tkIdx)', numTrials(tkIdx)');
                                    
                                    plot(StimLevelsFineGrain, fittingFunction(paramsValues, StimLevelsFineGrain),'LineWidth',1.1, 'color', useColor);
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
                                        blkLas = prc.combineBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==laserType & laserTrials & reactionTimes);
                                        if length(blkLas.trialEnd) < 50; continue; end
                                        fracRightTurns = prc.makeGrid(blkLas, blkLas.response==2, @mean, 1);
                                        noNanIdx = (tkIdx.*~isnan(fracRightTurns))>0;
                                        plot(visGrid((conflicts.*noNanIdx)>0),fracRightTurns((conflicts.*noNanIdx)>0),'^','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                        plot(visGrid((nonConflict.*noNanIdx)>0),fracRightTurns((nonConflict.*noNanIdx)>0),'o','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                        trialEst = floor(length(blkLas.response)/8);
                                    else
                                        plot(visGrid(conflicts),fracRightTurns(conflicts),'^','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                        plot(visGrid(nonConflict),fracRightTurns(nonConflict),'o','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                        trialEst = length(laserOff.response);
                                    end
                                    yAxLabel = 'Fraction of right turns';
                                case 'rea'
                                    yAxLabel = 'Rection time (ms)';
                                    reactionTimesGrid = prc.makeGrid(laserOff, laserOff.reactionTime, @median, 1)*1000-500;
                                    plot(visGrid(tkIdx), reactionTimesGrid(tkIdx),'LineWidth',1.1, 'color', useColor);
                                    plot(visGrid(conflicts), reactionTimesGrid(conflicts),'^','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                    plot(visGrid(nonConflict), reactionTimesGrid(nonConflict),'o','markerfacecolor', useColor, 'markeredgecolor', useColor);
                                    ylim([0 500])
                            end
                        end
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, trialEst, laserOff.nSessions));
                        set(gca, 'XTick', -45:15:45);
                end
            end
            figureSize = get(gcf, 'position');
            mainAxes = [60./figureSize(3:4) 1-2*(50./figureSize(3:4))];
            plt.suplabel('\fontsize{16} Auditory on {\color{red}right, \color{black}centre, or \color{blue}left}', 't', mainAxes);
            plt.suplabel(['\fontsize{16} ' yAxLabel], 'y', mainAxes);
            plt.suplabel('\fontsize{16} Visual Contrast', 'x', mainAxes);
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
                switch lower(plotType)
                    case 'amp'
                        blk = prc.combineBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==0);
                        percentCorrect = prc.makeGrid(blk, blk.feedback>0, @mean, 2);
                        xData = blk.audGrid(blk.audGrid>=0 & blk.visGrid==0);
                        yData = percentCorrect(blk.audGrid>=0 & blk.visGrid==0);
                        [~, pVal] = cellfun(@(x) ttest(x-0.5), yData);
                        plotColors = repmat('k', length(pVal),1);
                        plotColors(pVal<0.05) = 'k';
                        
                    case 'rig'
                        blk = prc.combineBlocks(obj.blocks{i}, vertcat(obj.blocks{i}.laserType)==0);
                        [rigsUsed,~,rigIdx] = unique(blk.rigName);
                        for j = 1:length(rigsUsed)
                            percentCorrect = prc.makeGrid(blk, blk.feedback>0, @mean, 2, find(rigIdx==j));
                            if strcmpi(plotType(4:end), 'vis')
                                visTest = obj.visValues{i}(find(obj.visValues{i}>0, 3));
                                yData{j,1} = cell2mat(percentCorrect(blk.audGrid==0 & blk.visGrid==visTest(end)));
                            else
                                audTest = obj.blocks{i}.audValues(find(obj.blocks{i}.audValues>0, 1));
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
                        blk = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserType==0);
                        trialType = prc.makeGrid(blk, blk.trialType, @mean, 3);
                        multiTrials = trialType==3;
                        percentCorrect = prc.makeGrid(blk, blk.feedback>0, @mean, 2);
                        if ~any(multiTrials); continue; end
                        numTrials = prc.makeGrid(blk, blk.feedback*0+1, @sum, 2);
                        testedPairs = {};
                        significance = [];
                        idx = 0;
                        for j = find(multiTrials)'
                            yData(idx+1,1) = percentCorrect(blk.grids.audValues==blk.grids.audValues(j) & trialType==1);
                            yData(idx+2,1) = percentCorrect(blk.grids.visValues==blk.grids.visValues(j) & trialType==2);
                            yData(idx+3,1) = percentCorrect(j);
                            xData(idx+1:idx+3,1) = {'Aud'; 'Vis'; 'Mul'};
                            xDataPosition(idx+1:idx+3,1) = floor((idx+1)/3)+(idx+1:idx+3);
                            [~, pairedTestSig(1,1)] = ttest(yData{idx+1}, yData{idx+3});
                            [~, pairedTestSig(2,1)] = ttest(yData{idx+2}, yData{idx+3});
                            significance = [significance; pairedTestSig];
                            testedPairs = [testedPairs; {xDataPosition([idx+1 idx+3]); xDataPosition([idx+2 idx+3])}];
                            idx = idx+3;
                        end
                        [~, pVal] = cellfun(@(x) ttest(x-0.5), yData);
                        plotColors = repmat('k', length(pVal),1);
                        plotColors(pVal<0.05) = 'k';
                        yLimits = [0.4 1.05];
                        sigPlot = 1;
                        figureSize = 450;
                end
                if ~exist('xDataPosition', 'var'); xDataPosition = 1:length(xData); end
                if ~exist('maxPlotLength', 'var'); maxPlotLength = length(xData)/2; end
                if ~exist('yLimits', 'var'); yLimits = [0 1]; end
                if ~exist('figureSize', 'var'); figureSize = 400; end
                obj.plt.getAxes(obj.subjects{i}, [80 70], maxPlotLength/4, [70 70], [60 20], obj.subjects(chosenSubjects), figureSize); hold on; box off;
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
            plt.suplabel(mainTitle, 't', mainAxes);
            plt.suplabel(mainYLabel, 'y', mainAxes);
            plt.suplabel(mainXLabel, 'x', mainAxes);
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
                    laserOff = prc.combineBlocks(cat(1,laserBlocks{:}), {vertcat(obj.blocks{i}.laserType)==0});
                    blkOn = prc.combineBlocks(cat(1,laserBlocks{:}), 'onlyLaser');
                    
                    if UniOrBilateral==2
                        laserOff.galvoPosition(:,1) = abs(laserOff.galvoPosition(:,1));
                        blkOn.galvoPosition(:,1) = abs(blkOn.galvoPosition(:,1));
                    end
                    galvoGrid = unique(blkOn.galvoPosition, 'rows');
                    blkOn.response(blkOn.reactionTime>diff(blkOn.laserOnOff, [], 2)) = 0;
                    laserOff.response(laserOff.reactionTime>diff(laserOff.laserOnOff, [], 2)) = 0;
                    
                    idx = 0;
                    yOpt = {'VisLeft'; 'VisRight'; 'ConflictVisLeft'; 'ConflictVisRight'};
                    tOpt = {'Timeout'; 'ChooseLeft'; 'ChooseRight'};
                    for conditions = [unique(blkOn.conditions(blkOn.trialType == 2))' (unique(blkOn.conditions(blkOn.trialType == 4))*-1)']
                        idx = idx+1;
                        idx2 = 0;
                        for response = 0:2
                            idx2 = idx2+1;
                            baseline = mean(laserOff.response(laserOff.conditions==conditions)==response);
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
                        blk = prc.combineBlocks(laserBlocks, 'none');
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
                            figure(figHand{j}); obj.plt.getAxes(obj.subjects{i}, [80 0], [], [70 70], [60 80]); hold on;
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
                            plt.suplabel(['\fontsize{16} Effects on ' trialLabel{j}], 't', mainAxes);
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
        

        
        function obj = changeMouse(obj, subjects, expDate, processingTag)
            %Get block and parameter files for the requested dates.
            [obj.blocks, obj.params]  = cellfun(@(x) getFilesFromDates(x, expDate, 'bloprm'), subjects, 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
            obj.expDate = expDate;
            
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
                    idx = 1:length(subjects);
                    obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x)), idx, 'uni', 0)';
                    obj.params = arrayfun(@(x) prc.combineParams(obj.params(subjectIdx==x)), idx, 'uni', 0)';
                    
                    b = obj.blocks;
                    reactionTimes = arrayfun(@(x) b{x}.reactionTime < diff(b{x}.laserOnOff, [], 2), idx, 'uni', 0);
                    obj.splits.offLaser = arrayfun(@(x) prc.combineBlocks(b{x}, b{x}.laserSession~=0 & b{x}.laserType==0 & reactionTimes{x}), idx, 'uni', 0);
                    obj.splits.noLaser = arrayfun(@(x) prc.combineBlocks(b{x}, b{x}.laserSession==0), idx, 'uni', 0);            
                    obj.splits.audLvls = arrayfun(@(y) arrayfun(@(x) prc.combineBlocks(b{y}, b{y}.audDiff==x), b{y}.audValues, 'uni', 0), idx, 'uni', 0);
                    obj.splits.visLvls = arrayfun(@(y) arrayfun(@(x) prc.combineBlocks(b{y}, b{y}.visDiff==x), b{y}.visValues, 'uni', 0), idx, 'uni', 0);                    
            end
            obj.currentAxes = [];
        end
        
        function replicatedObj = createMouseReplicate(obj, idx, subjectNames)
            replicatedObj = obj;
            replicatedObj.subjects = subjectNames;
            replicatedObj.blocks = repmat(replicatedObj.blocks(idx), length(replicatedObj.subjects),1);
            replicatedObj.params = repmat(replicatedObj.params(idx), length(replicatedObj.subjects),1);
        end
    end
    
    methods (Static)
        function [areaNames, galvoSites, laserType] = getSpecificGalvoParadigm(siteReference)
            switch lower(siteReference)
                case 'ina3'
                    areaNames = {'Left Visual';'Both Visual';'Right Visual';'Left M2';'Both M2';'Right M2';'Left Control';'Both Control';'Right Control'};
                    galvoSites = {[-2.5,-3]; [2.5,-3]; [2.5,-3]; [-1,2]; [1,2]; [1,2]; [-1.5,0]; [1.5,0]; [1.5,0]};
                    laserType = repmat([1;2;1],3,1);
            end
        end
    end
%                 reactionTimes = arrayfun(@(x) b{x}.reactionTime < diff(b{x}.laserOnOff, [], 2), idx, 'uni', 0);
%             obj.splits.offLaser = arrayfun(@(x) prc.combineBlocks(b{x}, b{x}.laserSession~=0 & b{x}.laserType==0 & reactionTimes{x}), idx, 'uni', 0);
%             obj.splits.noLaser = arrayfun(@(x) prc.combineBlocks(b{x}, b{x}.laserSession==0), idx, 'uni', 0);
%             obj.splits.audLvls = arrayfun(@(y) arrayfun(@(x) prc.combineBlocks(b{y}, b{y}.audDiff==x), b{y}.audValues, 'uni', 0), idx, 'uni', 0);
%             obj.splits.visLvls = arrayfun(@(y) arrayfun(@(x) prc.combineBlocks(b{y}, b{y}.visDiff==x), b{y}.visValues, 'uni', 0), idx, 'uni', 0);
%             
%             for i = 1:length(obj.splits.offLaser) 
%                 if isempty(obj.splits.offLaser{i}); obj.splits.maxNoLaser = obj.splits.noLaser;
%                 elseif isempty(obj.splits.noLaser{i}); obj.splits.maxNoLaser = obj.splits.offLaser;
%                 elseif length(obj.splits.offLaser{i}.trialStart) > length(obj.splits.noLaser{i}.trialStart);  obj.splits.maxNoLaser = obj.splits.offLaser;
%                 else, obj.splits.maxNoLaser = obj.splits.noLaser;
%                 end
%             end
end