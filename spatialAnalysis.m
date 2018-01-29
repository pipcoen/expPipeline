classdef spatialAnalysis
    %% spatialAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    % that can be used to plot various aspects of animal behavior. NOTE: This function is designed to operate on the output of convertExpFiles and
    % most methods are specific to multisensoty spatial integration plots.
    %
    % Inputs(default values) subjects({'PC011';'PC012';'PC013';'PC016'})------A cell array of subjects collect parameter and block files for
    % expDate('last')----------------------------------A cell array of dates, one for all subjects or one for each subject
    % combineMice(0)----------------------------A tag to indicate whether you want to combine data across mice or retain individuals.
    
    properties (Access=public)
        subjects;                %Cell array of subject names
        expDate;                 %Recording dates to use for each subject
        blocks;                  %Block files loaded for each subject (one cell per subject)
        params;                  %Parameter files loaded for each subject (one cell per subject)
        axesHandles;             %Handle to current axis being used for plotting
        figureHandles;           %Handle to current figure being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate, combineMice)
            % Initialize fields with default values if no vaules are provided. Then called changeMouse function to get requested data.
            if ~exist('subjects', 'var') || isempty(subjects); subjects = {'PC011';'PC012';'PC013';'PC015';'PC010';'PC017'}; end
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~exist('combineMice', 'var'); combineMice = 0; end
            if ~iscell(subjects); subjects = {subjects}; end
            if ~iscell(expDate); expDate = {expDate}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            obj = changeMouse(obj, subjects, expDate, 'bloprm', combineMice);
        end
        
        function viewBoxPlots(obj, plotType, alter)
            if ~exist('plotType', 'var'); plotType = 'mul'; end
            if ~exist('alter', 'var'); alter = 0; figure; else; axesOpt.reposition = 0; end
            if isgraphics(alter,'figure'); clf; end
            
            maxGrid = max(cell2mat(cellfun(@(x) [length(x.audValues) length(x.visValues)], obj.blocks, 'uni', 0)), [], 1);
            axesOpt.figureHWRatio = maxGrid(2)/(1.3*maxGrid(1));
            axesOpt.btlrMargins = [100 80 60 100];
            axesOpt.gapBetweenAxes = [100 40];
            
            boxPlot.colorMap = plt.redblue(64);
            boxPlot.axisLimits = [0 1];
            colorBar.colorLabel = 'Fraction of right turns';
            colorBar.colorDirection = 'normal';
            
            if isgraphics(alter,'figure') || ~alter; subjects2Run = 1:length(obj.subjects);
            elseif isaxes(alter); disp('NotFunctionalYet');
            end
            
            for i  = subjects2Run
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                boxPlot.subject = obj.subjects{i};
                boxPlot.trialNumber = length(normBlock.responseMade);
                boxPlot.nSessions = obj.blocks{i}.nSessions;
                boxPlot.xyValues = {normBlock.visValues*100; normBlock.audValues};
                boxPlot.xyLabel = {normBlock.audType; 'VisualContrast'};
                switch lower(plotType(1:3))
                    case 'res'
                        boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade==2, @mean, 1);
                        if isempty(obj.figureHandles) || ~ismember(obj.figureHandles, gcf); obj.figureHandles(end+1) = gcf; end
                        set(gcf, 'Tag', 'boxRes', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
                    case 'num'
                        normBlock = prc.combineBlocks(normBlock, normBlock.responseMade~=0);
                        boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade==2, @length, 1);
                        set(gcf, 'Tag', 'boxnum', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
                        colorBar.colorLabel = 'Relative Num of Trials';
                        boxPlot.axisLimits = [0 max(boxPlot.plotData(:))];
                    case 'rea'
                        boxPlot.plotData = prc.makeGrid(normBlock, round(normBlock.responseTime*1e3), @median, 1);
                        boxPlot.axisLimits = [min(boxPlot.plotData(:)) max(boxPlot.plotData(:))];
                end
                axesOpt.idx = i;
                axesOpt.totalNumOfAxes = length(obj.subjects);
                plt.getAxes(axesOpt);
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
            if ~exist('modelString', 'var'); modelString = 'SimpleLogistic'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [80 100 80 40];
            axesOpt.gapBetweenAxes = [100 60];
            for i  = 1:length(obj.subjects)
                axesOpt.idx = i;
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                obj.axesHandles = plt.getAxes(axesOpt);
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
        
        function viewDataWithoutFits(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [120 100 120 40];
            axesOpt.gapBetweenAxes = [100 60];
            axesOpt.figureSize = 500;
            for i  = 1:length(obj.subjects)
                axesOpt.idx = i;
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                plotOpt.Marker = '.'; plotOpt.MarkerSize = 20; plotOpt.lineStyle = '-';
                obj.axesHandles = plt.getAxes(axesOpt);
                switch plotType(1:3)
                    case 'rea'
                        gridData = prc.makeGrid(normBlock, round(normBlock.responseTime*1e3), @median, 1);
                        plt.gridSplitByRows(gridData, normBlock.visValues*100, normBlock.audValues, plotOpt);
                    case 'res'
                        gridData = prc.makeGrid(normBlock, normBlock.responseMade==2, @mean, 1);
                        plt.gridSplitByRows(gridData, normBlock.visValues*100, normBlock.audValues, plotOpt);
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
                        colorTake = 'rbck';
                        tmp.axesHandles = plt.getAxes(i, length(obj.subjects), [], [], [80 100 80 40], [100 60]);
                        for tTyp = 1:4
                            for response = 1:2
                                blk = prc.combineBlocks(normBlock, normBlock.trialType==tTyp & normBlock.responseMade == response & (abs(normBlock.visDiff) == abs(normBlock.visValues(3)) | tTyp == 1));
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
                        obj.axesHandles = plt.getAxes(i, length(obj.subjects), [], [], [80 100 80 40], [100 60]);
                        plt.psychoFits(normBlock);
                        plt.dataWithErrorBars(normBlock)
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(normBlock.responseMade), normBlock.nSessions));
                        xlim([min(normBlock.visValues), max(normBlock.visValues)]);
                    case 'flp'
                        obj.axesHandles = plt.getAxes(i, length(obj.subjects), [], [], [80 100 80 40], [100 60]);
                        plt.psychoFits(normBlock, normBlock.visValues);
                        plt.dataWithErrorBars(normBlock, 0, normBlock.visValues)
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(normBlock.responseMade), normBlock.nSessions));
                        xlim([min(normBlock.audValues), max(normBlock.audValues)]);
                    case 'vis'
                        if ~isempty(obj.expDate)
                            obj.runMouseReplicate({'UniLeft'; 'Bilateral'; 'UniRight'}, ['viewPsychoCurves(''' plotType ''')']);
                            return;
                        end
                        allresponseTimes = obj.blocks{i}.responseTime < diff(obj.blocks{i}.laserOnOff, [], 2);
                        normresponseTimes = normBlock.responseTime < diff(normBlock.laserOnOff, [], 2);
                        galvoPos = obj.blocks{i}.galvoPosition;
                        laserType = [1 2 1]; hemiMod = [-1 1 1];
                        laserBlock = prc.combineBlocks(obj.blocks{i}, obj.blocks{i}.laserType == laserType(i) & galvoPos(:,1) == 2.5*hemiMod(i)&allresponseTimes);
                        obj.axesHandles = plt.getAxes(i, length(obj.subjects), 500, 0.8, [80 100 80 40], [100 30]);
                        
                        plotOpt.lineStyle = '--';
                        plt.psychoFits(prc.combineBlocks(normBlock, normresponseTimes), [], plotOpt);
                        plt.psychoFits(laserBlock);
                        plt.dataWithErrorBars(laserBlock, 0)
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(laserBlock.responseMade), normBlock.nSessions));
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
                        plotOpt = prc.getMultiTriplets(normBlock, strcmpi(plotType, 'coh'));
                        plotOpt.figureSize = 400;
                        if strcmpi(plotType, 'coh')
                            plotOpt.mainTitle = '\fontsize{20} Fraction correct by condition';
                            plotOpt.mainYLabel = '\fontsize{20} Fraction correct';
                            plotOpt.yLimits = [0 1.2];
                        else
                            plotOpt.mainTitle = '\fontsize{20} Fraction of unisensory choices in conflict';
                            plotOpt.mainYLabel = '\fontsize{20} Fraction choices in unisensory direction';
                            plotOpt.yLimits = [0 1.2];
                        end
                        plotOpt.mainXLabel = '\fontsize{20} Condition';
                end
                obj.axesHandles = plt.getAxes(i, length(obj.subjects), plotOpt.figureSize, [], [100 100 100 20], [100 80]);
                plt.jitter(plotOpt.yData, plotOpt); grid('on');
                title(sprintf('%s: n = %d', obj.subjects{i}, obj.blocks{i}.nSessions));
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
            end
            figureSize = get(gcf, 'position');
            mainAxes = [80./figureSize(3:4) 1-2*(70./figureSize(3:4))];
            plt.suplabel(plotOpt.mainTitle, 't', mainAxes);
            plt.suplabel(plotOpt.mainYLabel, 'y', mainAxes);
            plt.suplabel(plotOpt.mainXLabel, 'x', mainAxes);
        end
        
        function viewInactivationResults(obj, plotType)
            imgBWOuline=imread('BrainOutlineBW.png');
            for i  = 1:length(obj.subjects)
                blkOff = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                blkOff = concatinateBlocks(blocks{i}, {vertcat(obj.blocks{i}.laserType)==0});
            end
            
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
            end
        end
        
        function obj = changeMouse(obj, subjects, expDate, dataType, combineMice)
            if ~exist('dataType', 'var'); dataType = 'bloprm'; end
            if ~exist('combineMice', 'var'); combineMice = 0; end
            %Get block and parameter files for the requested dates.
            [obj.blocks, obj.params, rawData]  = arrayfun(@(x,y) prc.getFilesFromDates(x, y{1}, dataType), subjects, expDate, 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
            rawData = vertcat(rawData{:});
            excessData = [isempty(obj.params), isempty(obj.blocks), isempty(rawData)];
            
            if combineMice
                combinedName = cell2mat(unique({obj.blocks.subject}')');
                [obj.blocks.subject] = deal(combinedName);
                [obj.params.subject] = deal(combinedName);
                obj.subjects = {combinedName};
                subjects = {combinedName};
            end
            
            maxLength = max([length(obj.blocks) length(obj.params) length(rawData)]);
            if excessData(1); obj.params = cell(maxLength, 1); end
            if excessData(2); obj.blocks = cell(maxLength, 1); end
            if excessData(3); rawData = cell(maxLength, 1); end
            obj.expDate = expDate;
            
            retainIdx = ones(length(obj.params),1)>0;
            %If there are multiple parameter sets in the requested date range for a mouse, use only the most common parameter set.
            for i = 1:length(subjects)
                if excessData(2); continue; end
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
            
            if ~excessData(1)
                obj.subjects = unique({obj.params.subject}');
                [~, subjectIdx] = ismember({obj.params.subject}', obj.subjects);
            elseif ~excessData(2)
                obj.subjects = unique({obj.blocks.subject}');
                [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
            end
            
            idx = 1:length(obj.subjects);
            if ~any(excessData(2:3)); obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x), [], rawData(subjectIdx==x)), idx, 'uni', 0)';
            elseif ~excessData(2), obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x)), idx, 'uni', 0)';
            end
            if ~excessData(1); obj.params = arrayfun(@(x) prc.combineParams(obj.params(subjectIdx==x)), idx, 'uni', 0)'; end
            obj.axesHandles = [];
            obj.figureHandles = [];
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
        function [normBlock, laserBlock] = getMaxNumberOfNormalTrials(block, inactivation)
            if ~exist('inactivation', 'var'); inactivation = 0; end
            if mode(block.laserSession>0)==1 || inactivation
                normBlock = prc.combineBlocks(block, block.laserSession~=0 & block.laserType==0 & block.responseMade~=0);
                laserBlock = prc.combineBlocks(block, block.laserSession~=0 & block.laserType~=0 & block.responseMade~=0);
            else, normBlock = prc.combineBlocks(block, block.laserSession==0 & block.responseMade~=0);
            end
        end
        
        function alterFigure(currentFigureHandle, ~)
            obj = get(currentFigureHandle.Number, 'UserData');
            switch lower(get(currentFigureHandle, 'Tag'))
                case 'boxres'
                    obj.viewBoxPlots('num', currentFigureHandle);
                case 'boxnum'
                    obj.viewBoxPlots('res', currentFigureHandle);
            end
        end
        
        function alterAxes(currentAxesHandle, ~)
            obj = get(currentAxesHandle.Number, 'UserData');
            switch lower(get(currentAxesHandle.Number, 'Tag'))
                case 'boxres'
                    obj.viewBoxPlots('num', currentAxesHandle);
            end
        end
    end
end