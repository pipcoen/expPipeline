classdef GLMmultiNest < fit.GLMmulti
    
    methods
        function obj = GLMmultiNest(inputBlockData, modelString)
            %% Input blockData must be a struct with fields: conditions and responseMade
            obj@fit.GLMmulti(inputBlockData);
            if exist('modelString', 'var'); obj.GLMMultiModels(modelString); end            
        end
        
        function figureHand = plotFit(obj)
            if size(obj.prmFits,1)~=1; error('Model not fitted (non-crossvalidated) yet'); end
            
            plotFit@fit.GLMmulti(obj);
            hold on;
            colorChoices = plt.selectRedBlueColors(obj.blockData.audValues);
            pHatCalculated = obj.calculatepHat(obj.prmFits,'eval');
            for audVal = obj.blockData.audValues(:)'
                plotIdx = obj.evalPoints(:,2)==audVal;
                plot(gca, obj.evalPoints(plotIdx,1)*obj.blockData.origMax(1), pHatCalculated(plotIdx,1), ...
                    'Color', colorChoices(obj.blockData.audValues==audVal,:), 'linewidth', 1.5, 'linestyle', '--');
            end
            figureHand = gca;
        end
        
        function fitCV(obj,nFolds)
            %Crossvalidated fitting
            if isempty(obj.modelString); error('Set model first'); end
            if ~exist('nFolds', 'var') || isempty(nFolds); nFolds = 10; end
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',2000);
            cvObj = cvpartition(obj.blockData.responseMade,'KFold',nFolds);
            obj.prmFits = nan(cvObj.NumTestSets,length(obj.prmLabels));
            obj.pHat = [];
            obj.logLik = nan(cvObj.NumTestSets,1);
            for i = 1:cvObj.NumTestSets
                cvTrainObj = copy(obj); cvTrainObj.blockData = prc.combineBlocks(cvTrainObj.blockData, cvObj.training(i));
                disp(['Model: ' obj.modelString '. Fold: ' num2str(i) '/' num2str(cvObj.NumTestSets)]);
                
                fittingObjective = @(b) (cvTrainObj.calculateLogLik(b));
                [obj.prmFits(i,:),~,exitflag] = fmincon(fittingObjective, obj.prmInit(), [], [], [], [], obj.prmBounds(1,:), obj.prmBounds(2,:), [], options);
                if ~any(exitflag == [1,2]); obj.prmFits(i,:) = nan(1,length(obj.prmLabels)); end
                
                cvTestObj = copy(obj); cvTestObj.blockData = prc.combineBlocks(cvTestObj.blockData, cvObj.test(i));
                pHatTested = cvTestObj.calculatepHat(obj.prmFits(i,:));
                obj.pHat(cvObj.test(i)) = pHatTested(sub2ind(size(pHatTested),(1:size(pHatTested,1))', cvTestObj.blockData.responseMade));
                obj.logLik(i) = mean(-log2(obj.pHat(cvObj.test(i))));
            end
        end

        
        function pHatCalculated = calculatepHat(obj, P, tag)
            if isempty(obj.modelString); error('Set model first'); end
            if ~exist('tag', 'var'); tag = 'runModel'; end
            [logOddsLR, logOddsTO] = obj.GLMMultiModels(tag, P);
            probTimeOut = exp(logOddsTO)./(1+exp(logOddsTO));
            probMakeResponse = 1-probTimeOut;
            probRespondRight = probMakeResponse.*(exp(logOddsLR)./(1+exp(logOddsLR)));
            probRespondLeft = probMakeResponse.*(1./(1+exp(logOddsLR)));
            pHatCalculated = [probTimeOut probRespondLeft probRespondRight];
        end
        
    end
end