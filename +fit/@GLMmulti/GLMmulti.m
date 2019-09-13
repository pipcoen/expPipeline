classdef GLMmulti < matlab.mixin.Copyable
    properties (Access=public)
        modelString;
        prmLabels;
        prmFits;
        prmBounds;
        prmInit;
        blockData;
        pHat;
        logLik;
        evalPoints;
        initGuess;
    end
    
    methods
        function obj = GLMmulti(inputBlockData, modelString)
            %% Input blockData must be a struct with fields: conditions and responseMade
            inputBlockData.origMax = [max(abs(inputBlockData.tri.stim.visDiff)) max(abs(inputBlockData.tri.stim.audDiff))];
            inputBlockData.visDiff = inputBlockData.tri.stim.visDiff;%./inputBlockData.origMax(1);
            inputBlockData.audDiff = inputBlockData.tri.stim.audDiff./inputBlockData.origMax(2);
            obj.blockData = inputBlockData;
            obj.blockData.selectedTrials = ones(size(inputBlockData.audDiff,1),1);
            tab = tabulate(obj.blockData.tri.outcome.responseMade)/100;
            obj.initGuess = sum(tab(:,3).*log2(tab(:,3)));
            if exist('modelString', 'var'); obj.GLMMultiModels(modelString); end
        end
        
        function fit(obj)
            %Non crossvalidated fitting
            if isempty(obj.modelString); error('Set model first'); end
            options = optimset('algorithm','interior-point','MaxFunEvals',100000,'MaxIter',10000);
            fittingObjective = @(b) (obj.calculateLogLik(b));
            [obj.prmFits,~,exitflag] = fmincon(fittingObjective, obj.prmInit, [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
            if ~any(exitflag == [1,2])
                obj.prmFits = nan(1,length(obj.prmLabels));
            end
            obj.pHat = obj.calculatepHat(obj.prmFits);
            obj.logLik = obj.calculateLogLik(obj.prmFits);%./length(obj.blockData.tri.outcome.responseMade);
        end
        
        function fitCV(obj,nFolds)
            %Crossvalidated fitting
            if isempty(obj.modelString); error('Set model first'); end
            if ~exist('nFolds', 'var') || isempty(nFolds); nFolds = 10; end
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',2000);
            cvObj = cvpartition(obj.blockData.tri.outcome.responseMade,'KFold',nFolds);
            obj.prmFits = nan(cvObj.NumTestSets,length(obj.prmLabels));
            obj.pHat = [];
            obj.logLik = nan(cvObj.NumTestSets,1);
            for i = 1:cvObj.NumTestSets
                cvTrainObj = copy(obj); cvTrainObj.blockData = prc.filtStruct(cvTrainObj.blockData, cvObj.training(i));
                disp(['Model: ' obj.modelString '. Fold: ' num2str(i) '/' num2str(cvObj.NumTestSets)]);
                
                fittingObjective = @(b) (cvTrainObj.calculateLogLik(b));
                [obj.prmFits(i,:),~,exitflag] = fmincon(fittingObjective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
                if ~any(exitflag == [1,2]); obj.prmFits(i,:) = nan(1,length(obj.prmLabels)); end
                
                cvTestObj = copy(obj); cvTestObj.blockData = prc.filtStruct(cvTestObj.blockData, cvObj.test(i));
                pHatTested = cvTestObj.calculatepHat(obj.prmFits(i,:));
                if min(cvTestObj.blockData.tri.outcome.responseMade) == 0; idxMod = 1; else, idxMod = 0; end
                obj.pHat(cvObj.test(i)) = pHatTested(sub2ind(size(pHatTested),(1:size(pHatTested,1))', cvTestObj.blockData.tri.outcome.responseMade+idxMod));
                obj.logLik(i) = -mean(log2(obj.pHat(cvObj.test(i))));
            end
        end
        
        function h = plotblockData(obj)
            %%
            numTrials = prc.makeGrid(obj.blockData, obj.blockData.tri.outcome.responseMade~=0, @length, 1);
            numRightTurns = prc.makeGrid(obj.blockData, obj.blockData.tri.outcome.responseMade==2, @sum, 1);
            
            audValues = [obj.blockData.audValues]./abs(max(obj.blockData.audValues));
            colorChoices = plt.selectRedBlueColors(audValues);
            
            [prob,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), numRightTurns, numTrials, 'uni', 0);
            prob = cell2mat(cellfun(@(x) permute(x, [3,1,2]), prob, 'uni', 0));
            lowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
            highBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));
            
            for audVal = audValues(:)'
                idx = find(sign(obj.blockData.grids.audValues)==audVal & numTrials>0);
                err = [prob(idx)-lowBound(idx), highBound(idx) - prob(idx)];
                errorbar(obj.blockData.grids.visValues(idx),prob(idx),err(:,1),err(:,2),'.','MarkerSize',20, 'Color', colorChoices(audValues==audVal,:));
                hold on;
            end
            xlabel('Contrast');
            ylabel('P( choice | contrast)');
            set(gca,'box','off');
            h=gca;
            set(gcf,'color','w');
        end
        
        function figureHand = plotFit(obj)
            if isempty(obj.prmFits); error('Model not fitted (non-crossvalidated) yet'); end
            params2use = mean(obj.prmFits,1);
            hold on;
            colorChoices = plt.selectRedBlueColors(obj.blockData.audValues);
            pHatCalculated = obj.calculatepHat(params2use,'eval');
            for audVal = obj.blockData.audValues(:)'
                plotIdx = obj.evalPoints(:,2)==audVal;
%                 plot(gca, obj.evalPoints(plotIdx,1)*obj.blockData.origMax(1), pHatCalculated(plotIdx,2), ...
                plot(gca, obj.evalPoints(plotIdx,1), pHatCalculated(plotIdx,2), ...
                    'Color', colorChoices(obj.blockData.audValues==audVal,:), 'linewidth', 2);
            end
            maxContrast =obj.blockData.origMax(1);
            xlim([-maxContrast maxContrast])
            set(gca, 'xTick', (-maxContrast):(maxContrast/4):maxContrast, 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)*100));
            title({cell2mat(unique(obj.blockData.exp.subject)'); obj.modelString});
            hold off;
            figureHand = gca;
        end
        
        
        function h = plotParams(obj)
            if size(obj.prmFits,1)~=1; return; end
            bar(obj.prmFits);
            set(gca,'XTickLabel',obj.prmLabels,'XTick',1:numel(obj.prmLabels));
            title(obj.modelString);
            h=gca;
        end
        
        function pHatCalculated = calculatepHat(obj, P, tag)
            if isempty(obj.modelString); error('Set model first'); end
            if ~exist('tag', 'var'); tag = 'runModel'; end
            logOdds = obj.GLMMultiModels(tag, P);
            probRight = exp(logOdds)./(1+exp(logOdds));
            pHatCalculated = [1-probRight probRight];
        end
        
        function logLik = calculateLogLik(obj,testParams)
            pHatCalculated = obj.calculatepHat(testParams);
            responseMade = obj.blockData.tri.outcome.responseMade+(min(obj.blockData.tri.outcome.responseMade(:)) == 0);
            logLik = -mean(log2(pHatCalculated(sub2ind(size(pHatCalculated),(1:size(pHatCalculated,1))', responseMade))));
        end
        
    end
end