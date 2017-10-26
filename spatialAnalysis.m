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
            colorBar.colorYTick = [0 1];
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
                [obj.currentAxes, obj.currentFigure] = plt.getAxes(i, length(obj.subjects), [], [], [60 80 60 20], [200 40]);
                glm = fit.GLMmulti(obj.blocks{i});
                glm = glm.setModel(modelString);
                glm = glm.fit;
                glm.plotFit;              
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
                        plt.dataWithErrorBars(normBlock, 0)
                        boxPlot.xyValues{1} = results.visValues;
                        boxPlot.axisLimits = [-0.25 1.25];
                        boxPlot.plotData = results.model{i};
                        if i == 3; boxPlot.axisLimits = [-0.25 0.25]; end
                end
            end
        end
        
        function viewPsychoCurves(obj, plotType)
            if ~exist('plotType', 'var'); plotType = 'res'; end
            figure;
            for i  = 1:length(obj.subjects)
                normBlock = spatialAnalysis.getMaxNumberOfNormalTrials(obj.blocks{i});
                switch plotType(1:3)
                    case 'res'
                        obj.currentAxes = plt.getAxes(i, length(obj.subjects), [], [], [60 80 60 20], [200 40]);
                        plt.psychoFits(normBlock);
                        plt.dataWithErrorBars(normBlock, 0)
                        title(sprintf('%s: %d Tri, %d Sess', obj.subjects{i}, length(normBlock.response), normBlock.nSessions));
                        set(gca, 'XTick', -45:15:45);
                end
            end
            figureSize = get(gcf, 'position');
            mainAxes = [60./figureSize(3:4) 1-2*(50./figureSize(3:4))];
            plt.suplabel('\fontsize{16} Auditory on {\color{red}right, \color{black}centre, or \color{blue}left}', 't', mainAxes);
            plt.suplabel(['\fontsize{16} ' yAxLabel], 'y', mainAxes);
            plt.suplabel('\fontsize{16} Visual Contrast', 'x', mainAxes);
        end    
        
        function obj = changeMouse(obj, subjects, expDate)
            %Get block and parameter files for the requested dates.
            [obj.blocks, obj.params]  = cellfun(@(x) prc.getFilesFromDates(x, expDate, 'bloprm'), subjects, 'uni', 0);
            obj.blocks = vertcat(obj.blocks{:});
            obj.params = vertcat(obj.params{:});
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
            obj.subjects = unique({obj.blocks.subject}');
            
            [~, subjectIdx] = ismember({obj.blocks.subject}', obj.subjects);
            idx = 1:length(subjects);
            obj.blocks = arrayfun(@(x) prc.combineBlocks(obj.blocks(subjectIdx==x)), idx, 'uni', 0)';
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