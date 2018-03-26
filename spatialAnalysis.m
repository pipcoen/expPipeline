 classdef spatialAnalysis < matlab.mixin.Copyable
    %% spatialAnalysis object that extracts beahvioral information for a specified animal or set of animals. The resulting object has a set of methods
    % that can be used to plot various aspects of animal behavior. NOTE: This function is designed to operate on the output of convertExpFiles and
    % most methods are specific to multisensoty spatial integration plots.
    %
    % Inputs(default values) subjects({'PC011';'PC012';'PC013';'PC016'})------A c   ell array of subjects collect parameter and block files for
    % expDate('last')----------------------------------A cell array of dates, one for all subjects or one for each subject
    % combineMice(0)----------------------------A tag to indicate whether you want to combine data across mice or retain individuals.
    
    properties (Access=public) 
        subjects;                %Cell array of subject names
        expDate;                 %Recording dates to use for each subject
        blocks;                  %Block files loaded for each subject (one cell per subject)
        params;                  %Parameter files loaded for each subject (one cell per subject)
        glmFit;
        axesHandles;             %Handle to current axis being used for plotting
        figureHandles;           %Handle to current figure being used for plotting
    end
    
    %%
    methods
        function obj = spatialAnalysis(subjects, expDate, combineMice, dataType)
            % Initialize fields with default values if no vaules are provided. Then called changeMouse function to get requested data.
            if ~exist('subjects', 'var') || isempty(subjects); subjects = {'PC011';'PC012';'PC013';'PC015';'PC010';'PC017'}; end
            if ~exist('expDate', 'var'); expDate = 'last'; end
            if ~exist('combineMice', 'var'); combineMice = 0; end
            if ~exist('dataType', 'var'); dataType = 'bloprm'; end
            if ~iscell(subjects); subjects = {subjects}; end
            if ~iscell(expDate); expDate = {expDate}; end
            if length(expDate) < length(subjects); expDate = repmat(expDate, length(subjects),1); end
            obj = changeMouse(obj, subjects, expDate, dataType, combineMice);
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
                [normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                boxPlot.subject = obj.subjects{i};
                boxPlot.trialNumber = length(normBlock.responseMade);
                boxPlot.nSessions = obj.blocks{i}.nSessions;
                boxPlot.xyValues = {normBlock.visValues*100; normBlock.audValues};
                boxPlot.xyLabel = {normBlock.audType; 'VisualContrast'};
                switch lower(plotType(1:3))
                    case 'res'
                        boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade==2, @mean);
                        if isempty(obj.figureHandles) || ~ismember(obj.figureHandles, gcf); obj.figureHandles(end+1) = gcf; end
                        set(gcf, 'Tag', 'boxRes', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
                    case 'gng'
                        [~,normBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i}, 1, -1);
                        normBlock = prc.combineBlocks(normBlock, normBlock.timeOutsBeforeResponse==0);
                        set(gcf, 'Tag', 'boxGNG', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
                        boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade~=0, @mean);
                    case 'num'
                        boxPlot.plotData = prc.makeGrid(normBlock, normBlock.responseMade==2, @length);
                        set(gcf, 'Tag', 'boxNum', 'userData', obj, 'ButtonDownFcn', @spatialAnalysis.alterFigure);
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
                normBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                obj.axesHandles = plt.getAxes(axesOpt);
                obj.glmFit{i} = fit.GLMmulti(normBlock);    
                obj.glmFit{i}.setModel(modelString);
                obj.glmFit{i}.fit;
                obj.glmFit{i}.plotFit;
                hold on; box off;
                plt.dataWithErrorBars(normBlock);
                xL = xlim; hold on; plot(xL,[0.5 0.5], '--k', 'linewidth', 1.5);
                yL = ylim; hold on; plot([0 0], yL, '--k', 'linewidth', 1.5);
            end
        end
        
        function getGLMPerturbations(obj, modelString)
            if ~exist('modelString', 'var'); modelString = 'SqrtLogisticSplit'; end
            if ~isempty(obj.expDate)
                runMouseReplicate(copy(obj), {'Bias', 'VR-VL', 'AR-AL'}, 'getGLMPerturbations')
                return;
            end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins =  [50 100 10 10];
            axesOpt.gapBetweenAxes = [40 0];
            axesOpt.numOfRows = 3;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.figureHWRatio = 0.8;
            axesOpt.figureSize = 400;
            for i  = 1:length(obj.subjects)
                axesOpt.idx = i;
                obj.axesHandles = plt.getAxes(axesOpt);
                [normBlock, laserBlock] = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                obj.glmFit{i} = fit.GLMmulti(normBlock);
                obj.glmFit{i}.setModel(modelString);
                obj.glmFit{i}.fit;
                [galvoBlocks, gridXY] = prc.makeGrid(laserBlock, laserBlock, [], 'galvouni', 2);
                plotData = nan([size(galvoBlocks), length(obj.glmFit{i}.prmInit)]);
                for j = find(~cellfun(@isempty, galvoBlocks))'
                    if length(galvoBlocks{j}.responseMade)<100; continue; end
                    galvoBlocks{j}.prmInit = obj.glmFit{i}.prmFits;
                    subGLM = fit.GLMmulti(galvoBlocks{j});
                    subGLM.setModel(modelString);
                    hold on; box off;
                    subGLM.fit;
                    [xIdx, yIdx] = ind2sub(size(galvoBlocks), j);
                    plotData(xIdx, yIdx, :) = subGLM.prmFits-(obj.glmFit{i}.prmFits);
                    numTrials(xIdx, yIdx, 1) = length(galvoBlocks{j}.responseMade);
                end
                plt.inactivatedBrainPerturbations(plotData, gridXY, obj.subjects{i}, numTrials);
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
                normBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
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
        
        function viewJitterPlot(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'coh'; end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins = [100 100 100 20];
            axesOpt.gapBetweenAxes = [100 80];
            axesOpt.figureSize = 500;
            for i  = 1:length(obj.subjects)
                axesOpt.idx = i;
                switch lower(plotType)
                    case {'coh'; 'con'}
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
                    case 'whe'
                        normBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{i});
                        plotOpt = prc.getCoherentConflictPairs(normBlock);
                        plotOpt.mainXLabel = '\fontsize{20} Condition';
                        plotOpt.mainTitle = '\fontsize{20} Difference in reation times for coherent/conflict conditions (ms)';
                        plotOpt.mainYLabel = '\fontsize{20} Time to wheel movement (ms)';
                end
                obj.axesHandles = plt.getAxes(axesOpt);
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
        
        function viewInactivationResults(obj)
            if ~isempty(obj.expDate)
                runMouseReplicate(copy(obj), {'VisUni(L-R)', 'AudUni(L-R)', 'CohUni(VL-VR)', 'ConUni(AL-AR)'}, 'viewInactivationResults')
                return;
            end
            figure;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.btlrMargins =  [50 100 10 10];
            axesOpt.gapBetweenAxes = [40 0];
            axesOpt.numOfRows = 2;
            axesOpt.totalNumOfAxes = length(obj.subjects);
            axesOpt.figureHWRatio = 0.8;
            axesOpt.figureSize = 400;
            respBlock = spatialAnalysis.getMaxNumberOfTrials(obj.blocks{1}, 2);
            for i  = 1:length(obj.subjects)
                [brainPlot.plotData, brainPlot.pVals, brainPlot.gridXY] = fit.testLaserEffectAtEachSite(respBlock, obj.subjects{i});
                axesOpt.idx = i;
                obj.axesHandles = plt.getAxes(axesOpt);
                brainPlot.title = obj.subjects{i};
                plt.scanningBrainEffects(brainPlot);
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
            for i = 1:length(obj.subjects)
                replicatedObj = copy(obj);
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
        function [normBlock, laserBlock] = getMaxNumberOfTrials(block, inactivation, excludedResponses)
            if ~exist('inactivation', 'var'); inactivation = 0; end
            if ~exist('excludedResponses', 'var'); excludedResponses = 0; end
            if inactivation == 2
                normBlock = prc.combineBlocks(block, block.laserSession~=0 & block.responseMade~=excludedResponses);
            elseif mode(block.laserSession>0)==1 || inactivation
                normBlock = prc.combineBlocks(block, block.laserSession~=0 & block.laserType==0 & block.responseMade~=excludedResponses);
                laserBlock = prc.combineBlocks(block, block.laserSession~=0 & block.laserType~=0 & block.responseMade~=excludedResponses);
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