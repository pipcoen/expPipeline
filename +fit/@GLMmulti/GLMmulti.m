classdef GLMmulti < matlab.mixin.Copyable
    properties (Access=public)
        modelString;
        prmLabels;
        prmFits;
        prmBounds;
        prmInit;
        blockData;
        pHat;
        evalPoints;
        initGuess;
    end
    
    methods
        function obj = GLMmulti(inputblockData)
            %% Input blockData must be a struct with fields: conditions and responseMade
            inputblockData.origMax = [max(abs(inputblockData.visDiff)) max(abs(inputblockData.audDiff))];
            inputblockData.visDiff = inputblockData.visDiff./inputblockData.origMax(1);
            inputblockData.audDiff = inputblockData.audDiff./inputblockData.origMax(2);
            obj.blockData = inputblockData;
            tab = tabulate(obj.blockData.responseMade)/100;
            obj.initGuess = sum(tab(:,3).*log2(tab(:,3)));
        end
        
        function obj = fit(obj)
            %Non crossvalidated fitting
            if isempty(obj.modelString); error('Set model first'); end
            options = optimset('algorithm','interior-point','MaxFunEvals',100000,'MaxIter',10000);
            fittingObjective = @(b) (obj.calculateLogLik(b));
            [obj.prmFits,~,exitflag] = fmincon(fittingObjective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
            if ~any(exitflag == [1,2])
                obj.prmFits = nan(1,length(obj.prmLabels));
            end
        end
        
        function [obj,logLik] = fitCV(obj,nFolds)
            %Crossvalidated fitting
            if isempty(obj.modelString); error('Set model first'); end
            if ~exist('nFolds', 'var') || isempty(nFolds); nFolds = length(obj.blockData.responseMade); end
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',2000);
            cvObj = cvpartition(obj.blockData.responseMade,'KFold',nFolds);
            
            obj.prmFits = nan(cvObj.NumTestSets,length(obj.prmLabels));
            obj.pHat = [];
            allInputs = obj.inputFun(obj.blockData);
            logLik = nan(cvObj.NumTestSets,1);
            for i = 1:cvObj.NumTestSets
                disp(['Model: ' obj.modelString '. Fold: ' num2str(i) '/' num2str(cvObj.NumTestSets)]);
                
                fittingObjective = @(b) (obj.calculateLogLik(b,  allInputs(cvObj.training(i),:), obj.blockData.responseMade(cvObj.training(i))));
                [obj.prmFits(i,:),~,exitflag] = fmincon(fittingObjective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
                if ~any(exitflag == [1,2]); obj.prmFits(i,:) = nan(1,length(obj.prmLabels)); end
                
                pHatTrained = obj.calculatepHat(obj.prmFits(i,:), allInputs(cvObj.test(i),:));
                
                testIdx = sub2ind(size(pHatTrained), 1:size(pHatTrained, 1), obj.blockData.responseMade(cvObj.test(i))');
                obj.pHat(cvObj.test(i)) = pHatTrained(testIdx);
                logLik(i)=mean(-log2(obj.pHat(cvObj.test(i))));
            end
        end
        
        function h = plotblockData(obj)
            %%
            numTrials = prc.makeGrid(obj.blockData, obj.blockData.responseMade~=0, @length, 1);
            numRightTurns = prc.makeGrid(obj.blockData, obj.blockData.responseMade==2, @sum, 1);
            
            audValues = obj.blockData.audValues;
            colorChoices = plt.selectRedBlueColors(audValues);
            
            [prob,confInterval] = arrayfun(@(x,z) binofit(x, z, 0.05), numRightTurns, numTrials, 'uni', 0);
            prob = cell2mat(cellfun(@(x) permute(x, [3,1,2]), prob, 'uni', 0));
            lowBound = cell2mat(cellfun(@(x) permute(x(:,1), [3,2,1]), confInterval, 'uni', 0));
            highBound = cell2mat(cellfun(@(x) permute(x(:,2), [3,2,1]), confInterval, 'uni', 0));
            
            for audVal = audValues(:)'
                idx = find(obj.blockData.grids.audValues==audVal & numTrials>0);
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
            if size(obj.prmFits,1)~=1; error('Model not fitted (non-crossvalidated) yet'); end
            hold on;
            colorChoices = plt.selectRedBlueColors(obj.blockData.audValues);
            pHatCalculated = obj.calculatepHat(obj.prmFits,'eval');
            for audVal = obj.blockData.audValues(:)'
                plotIdx = obj.evalPoints(:,2)==audVal;
                plot(gca, obj.evalPoints(plotIdx,1)*obj.blockData.origMax(1), pHatCalculated(plotIdx,2), ...
                    'Color', colorChoices(obj.blockData.audValues==audVal,:), 'linewidth', 2);
            end
            maxContrast =obj.blockData.origMax(1);
            xlim([-maxContrast maxContrast])
            set(gca, 'xTick', (-maxContrast):(maxContrast/4):maxContrast, 'xTickLabel', round(((-maxContrast):(maxContrast/4):maxContrast)*100));
            title(obj.modelString);
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
            logLik = -sum(log2(pHatCalculated(:,1).*(obj.blockData.responseMade==1)+pHatCalculated(:,2).*(obj.blockData.responseMade==2)));
        end
        
    end
end