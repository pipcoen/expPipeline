classdef GLMmultiNest < fit.GLMmulti
    
    methods
        function obj = GLMmultiNest(inputBlockData, modelString)
            %% Input blockData must be a struct with fields: conditions and responseMade
            obj@fit.GLMmulti(inputBlockData);
            if exist('modelString', 'var'); obj.GLMMultiModels(modelString); end            
        end
        
        function figureHand = plotFit(obj)
            if isempty(obj.prmFits); error('Model not fitted (non-crossvalidated) yet'); end
            params2use = mean(obj.prmFits,1);
            plotFit@fit.GLMmulti(obj);
%             hold on;
%             colorChoices = plt.selectRedBlueColors(obj.blockData.audValues);
%             pHatCalculated = obj.calculatepHat(params2use,'eval');
%             for audVal = obj.blockData.audValues(:)'
%                 plotIdx = obj.evalPoints(:,2)==audVal;
%                 plot(gca, obj.evalPoints(plotIdx,1)*obj.blockData.origMax(1), pHatCalculated(plotIdx,1), ...
%                     'Color', colorChoices(obj.blockData.audValues==audVal,:), 'linewidth', 1.5, 'linestyle', '--');
%             end
            figureHand = gca;
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