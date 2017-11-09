function [paramsValues, fittingFunction] = psychoCurve(xData, numChoices, numTrials, fittingFunction)
if ~exist('fittingFunction', 'var'); fittingFunction = @PAL_Logistic; end
PF = fittingFunction;
searchGrid.alpha = -5:0.1:15; 
searchGrid.beta = 0:0.05:5; 
searchGrid.gamma = 0:0.025:0.5; 
searchGrid.lambda = 0:0.025:0.5;
paramsValues = PAL_PFML_Fit(xData(:), numChoices(:), numTrials(:), searchGrid, [1,1,1,1], PF,'lapseLimits',[0 1],'guessLimits', [0  1]);