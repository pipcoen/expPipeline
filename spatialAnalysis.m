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
        currentAxes;             %Handle to current axis being used for plotting
        currentFigure;           %Handle to current figure being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate)
            % Initialize fields with default values if no vaules are provided.
            if ~exist('subjects', 'var') || isempty(subjects); subjects = {'PC011';'PC012';'PC013';'PC015';'PC010';'PC017'}; end
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~iscell(subjects); subjects = {subjects}; end
            if ~iscell(expDate); expDate = {expDate}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            obj = changeMouse(obj, subjects, expDate);
        end
        
        function viewBoxPlots(obj, plotType)
            figure;
            if ~exist('plotType', 'var'); plotType = 'res'; end
            maxGrid = max(cell2mat(cellfun(@(x) [length(x.audValues) length(x.visValues)], obj.blocks, 'uni', 0)), [], 1);
            boxPlot.colorMap = plt.redblue(64);
            boxPlot.axisLimits = [0 1];
            
            colorBar.colorLabel = 'Fraction of right turns';
            colorBar.colorDirection = 'normal';
            
            for i  = 1:length(obj.subjects)
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                boxPlot.subject = obj.subjects{i};
                boxPlot.trialNumber = length(normBlock.response);
                boxPlot.nSessions = obj.blocks{i}.nSessions;
                boxPlot.xyValues = {normBlock.visValues*100; normBlock.audValues};
                boxPlot.xyLabel = {normBlock.audType; 'VisualContrast'};
                switch lower(plotType(1:3))
                    case 'res'
                        boxPlot.plotData = prc.makeGrid(normBlock, normBlock.response==2, @mean, 1);
                    case {'svd'; 'mul'}
                        if ~isempty(obj.expDate)
                            obj.runMouseReplicate({'Original'; 'Model'; 'Difference'}, ['viewBoxPlots(''' plotType ''')']);
                            return;
                        end
                        results = fit.outerProduct(normBlock);
                        if contains(plotType, 'svd'); results = results.svd; else; results = results.mul; end
                        boxPlot.xyValues{1} = results.visValues;
                        boxPlot.axisLimits = [-0.25 1.25];
                        boxPlot.plotData = results.model{i};
                        if i == 3; boxPlot.axisLimits = [-0.25 0.25]; end
                    case 'rea'
                        boxPlot.plotData = prc.makeGrid(normBlock, round(normBlock.reactionTime*1000), @median, 1);
                        boxPlot.axisLimits = [min(boxPlot.plotData(:)) max(boxPlot.plotData(:))];
                end
                plt.getAxes(i, length(obj.subjects), [], maxGrid(2)/(1.3*maxGrid(1)), [100 80 60 100], [100 40]);
                plt.boxPlot(boxPlot);
                colorBar.colorYTick = {'Min'; 'Max'};
            end
            currentAxisPotision = get(gca, 'position');
            figureSize = get(gcf, 'position');
            
            colorBar.handle = colorbar;
            set(colorBar.handle, 'Ticks', get(colorBar.handle, 'Limits'), 'TickLabels', colorBar.colorYTick, 'YDir', colorBar.colorDirection);
            set(gca, 'position', currentAxisPotision);
            colorBar.textHandle = ylabel(colorBar.handle, colorBar.colorLabel);
            set(colorBar.textHandle, 'position', [1 mean(get(colorBar.handle, 'Limits')) 0], 'FontSize', 14)
            set(colorBar.handle, 'position', [1-75/figureSize(3), 0.2, 30/figureSize(3), 0.6])
        end
        
        function viewGLMFit(obj, modelString)
            if ~exist('modelString', 'var'); modelString = 'C-subset-multiAud2'; end
            figure;
            for i  = 1:length(obj.subjects)
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                obj.currentAxes = plt.getAxes(i, length(obj.subjects), [], [], [80 100 80 40], [100 60]);
                glm = fit.GLMmulti(normBlock);
                glm = glm.setModel(modelString);
                glm = glm.fit;
                glm.plotFit;
                hold on; box off;
                plt.dataWithErrorBars(normBlock);
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
            end
        end
        
        function viewAltCurves(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'mul'; end
            figure;
            for i  = 1:length(obj.subjects)
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                switch plotType(1:3)
                    case {'svd'; 'mul'}
                        results = fit.outerProduct(normBlock);
                        if contains(plotType, 'svd'); results = results.svd; else; results = results.mul; end
                        results.model{1}(isnan(results.model{2})) = nan;
                        results.model{1}(isnan(results.model{1})) = nan;
                        obj.currentAxes = plt.getAxes(i, length(obj.subjects), [], [], [80 100 80 40], [100 60]);
                        plotOpt.Marker = '.'; plotOpt.MarkerSize = 20; plotOpt.lineStyle = 'none';
                        plt.gridSplitByRows(results.model{1}, results.visValues/100, normBlock.audValues, plotOpt);
                        plotOpt.Marker = '^'; plotOpt.MarkerSize = 5; plotOpt.lineWidth = 1.5; plotOpt.lineStyle = '-';
                        plt.gridSplitByRows(results.model{2}, results.visValues/100, normBlock.audValues, plotOpt);
                        title(sprintf('%s: %d Sess', obj.subjects{i}, normBlock.nSessions));
                end
                figureSize = get(gcf, 'position');
                mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
                plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
                plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
            end
        end
        
        function viewWheelMovement(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'mov'; end
            figure;
            for i  = 1:length(obj.subjects)
                tmp = obj.changeMouse({obj.subjects{i}}, {obj.expDate{i}}, 1);
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(tmp.blocks{1});
                switch plotType(1:3)
                    case {'mov'}
                        timeSeg = -0.1:0.01:0.5;
                        colorTake = 'rrkr';
                        tmp.currentAxes = plt.getAxes(i, length(obj.subjects), [], [], [80 100 80 40], [100 60]);
                        for tTyp = 1:4
                            for response = 1:2
                            blk = prc.combineBlocks(normBlock, normBlock.trialType==tTyp & normBlock.response == response & (abs(normBlock.visDiff) == abs(normBlock.visValues(3)) | tTyp == 1));
                            wheelMove = cellfun(@double, blk.wheelTimeValue(~cellfun(@isempty, blk.wheelTimeValue)), 'uni', 0);
                            wheelMove = cell2mat(cellfun(@(x) interp1(x(:,1), x(:,2)/350, timeSeg, 'nearest'), wheelMove, 'uni', 0));
                            plot(timeSeg,nanmean(wheelMove), colorTake(tTyp), 'linewidth', 2);
                            hold on;
                            end
                        end
                        xlim([timeSeg(1) timeSeg(end)]);
                        ylim([-1 1]);
                        set(gca, 'yTick', -1:0.5:1, 'yTickLabel', -360:180:360);
                        plot(xlim, [0 0], '--k', 'linewidth', 1.5);
                        box off;
                        title(sprintf('%s: %d Sess', obj.subjects{i}, normBlock.nSessions));

                end

            end
                            figureSize = get(gcf, 'position');

                mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
                plt.suplabel('\fontsize{20} Wheel Movement', 't', mainAxes);
                plt.suplabel('\fontsize{20} Degrees pf rotation', 'y', mainAxes);
                plt.suplabel('\fontsize{20} Time from stimulus on (s)', 'x', mainAxes);
        end
        
        function viewPsychoCurves(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            for i  = 1:length(obj.subjects)
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                switch plotType(1:3)
                    case 'res'
                        obj.currentAxes = plt.getAxes(i, length(obj.subjects), [], [], [80 100 80 40], [100 60]);
                        plt.psychoFits(normBlock);
                        plt.dataWithErrorBars(normBlock)
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(normBlock.response), normBlock.nSessions));
                        xlim([min(normBlock.visValues), max(normBlock.visValues)]);
                    case 'vis'
                        if ~isempty(obj.expDate)
                            obj.runMouseReplicate({'UniLeft'; 'Bilateral'; 'UniRight'}, ['viewPsychoCurves(''' plotType ''')']);
                            return;
                        end
                        allReactionTimes = obj.blocks{i}.reactionTime < diff(obj.blocks{i}.laserOnOff, [], 2);
                        normReactionTimes = normBlock.reactionTime < diff(normBlock.laserOnOff, [], 2);
                        galvoPos = obj.blocks{i}.galvoPosition;
                        laserType = [1 2 1]; hemiMod = [-1 1 1];
                        laserBlock = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserType == laserType(i) & galvoPos(:,1) == 2.5*hemiMod(i)&allReactionTimes);
                        obj.currentAxes = plt.getAxes(i, length(obj.subjects), 500, 0.8, [80 100 80 40], [100 30]);
                        plt.psychoFits(prc.combineBlocks(normBlock, normReactionTimes));
                        plt.dataWithErrorBars(laserBlock)
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(laserBlock.response), normBlock.nSessions));
                        xlim([min(normBlock.visValues), max(normBlock.visValues)]);
                        %                         set(gca, 'XTick', min(normBlock.visValues):1:max(normBlock.visValues));
                end
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel('\fontsize{20} Auditory on {\color{red}right, \color{gray}centre, or \color{blue}left}', 't', mainAxes);
            plt.suplabel('\fontsize{20} Fraction of right choices', 'y', mainAxes);
            plt.suplabel('\fontsize{20} Visual Contrast', 'x', mainAxes);
        end
        
        function viewJitterPlot(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'coh'; end
            figure;
            for i  = 1:length(obj.subjects)
                switch lower(plotType)
                    case {'coh'; 'con'}
                        normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                        jPlot = prc.getMultiTriplets(normBlock, strcmpi(plotType, 'coh'));
                        
                        jPlot.figureSize = 400;
                        if strcmpi(plotType, 'coh')
                            jPlot.mainTitle = '\fontsize{20} Fraction correct by condition';
                            jPlot.mainYLabel = '\fontsize{20} Fraction correct';
                            jPlot.yLimits = [0.2 1.05]; 
                        else
                            jPlot.mainTitle = '\fontsize{20} Fraction of unisensory choices in conflict';
                            jPlot.mainYLabel = '\fontsize{20} Fraction choices in unisensory direction';
                            jPlot.yLimits = [0 1.05];
                        end
                        jPlot.mainXLabel = '\fontsize{20} Condition';
                    case {'mod'}
                        normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                        jPlot = prc.compareModels(normBlock);
                        
                        jPlot.figureSize = 400;
                        if strcmpi(plotType, 'coh')
                            jPlot.mainTitle = '\fontsize{20} Fraction correct by condition';
                            jPlot.mainYLabel = '\fontsize{20} Fraction correct';
                            jPlot.yLimits = [0.2 1.05];
                        else
                            jPlot.mainTitle = '\fontsize{20} Fraction of unisensory choices in conflict';
                            jPlot.mainYLabel = '\fontsize{20} Fraction choices in unisensory direction';
                            jPlot.yLimits = [0 1.05];
                        end
                        jPlot.mainXLabel = '\fontsize{20} Condition';
                end
                obj.currentAxes = plt.getAxes(i, length(obj.subjects), jPlot.figureSize, [], [100 100 100 20], [100 80]);
                jitterPlot(jPlot.yData, 'xPos', jPlot.xPosition); grid('on');
                set(gca, 'xTickLabels', jPlot.xLabels);
                title(sprintf('%s: n = %d', obj.subjects{i}, obj.blocks{i}.nSessions));
                ylim(jPlot.yLimits);
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                
                if ~isempty(jPlot.testedPairs); plt.sigstar(jPlot.testedPairs, jPlot.significance); end
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel(jPlot.mainTitle, 't', mainAxes);
            plt.suplabel(jPlot.mainYLabel, 'y', mainAxes);
            plt.suplabel(jPlot.mainXLabel, 'x', mainAxes);
        end
        
        function obj = changeMouse(obj, subjects, expDate, includeRaw)
            if ~exist('includeRaw', 'var'); includeRaw = 0; end
            %Get block and parameter files for the requested dates.
            if ~includeRaw; [obj.blocks, obj.params]  = arrayfun(@(x,y) prc.getFilesFromDates(x, y{1}, 'bloprm'), subjects, expDate, 'uni', 0); rawData = {};
            else, [obj.blocks, obj.params, rawData]  = arrayfun(@(x,y) prc.getFilesFromDates(x, y{1}, 'bloprmraw'), subjects, expDate, 'uni', 0);
            end
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
            rawData = vertcat(rawData{:});
            if ~includeRaw; rawData = cell(length(obj.params), 1); end
            obj.expDate = expDate;
            
            retainIdx = ones(length(obj.params),1)>0;
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
            rawData = rawData(retainIdx);
            obj.subjects = unique({obj.blocks.subject}');
            
            [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
            idx = 1:length(subjects);
            if includeRaw; obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x), [], rawData(subjectIdx==x)), idx, 'uni', 0)';
            else, obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x)), idx, 'uni', 0)';
            end
            obj.params = arrayfun(@(x) prc.combineParams(obj.params(subjectIdx==x)), idx, 'uni', 0)';
            obj.currentAxes = [];
            obj.currentFigure = [];
        end
        
        function runMouseReplicate(obj, subjectNames, funString)
            close;
            for i = 1:length(obj.subjects)
                replicatedObj = obj;
                replicatedObj.subjects = subjectNames;
                replicatedObj.blocks = repmat(replicatedObj.blocks(i), length(replicatedObj.subjects),1);
                replicatedObj.params = repmat(replicatedObj.params(i), length(replicatedObj.subjects),1);
                replicatedObj.expDate = [];
                eval(['replicatedObj.' funString ';']);
                suptitle(obj.subjects(i));
            end
            
        end
    end
    
    methods (Static)
        function normBlock = getMaxNumberOfNormalTrials(block)
            if mode(block.laserSession>0)==1
                normBlock = prc.combineBlocks(block, block.laserSession~=0 & block.laserType==0);
            else, normBlock = prc.combineBlocks(block, block.laserSession==0);
            end
        end
    end
end