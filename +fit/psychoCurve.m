function [paramsValues, fittingFunction] = psychoCurve(xData, numChoices, numTrials, fittingFunction)
if ~exist('fittingFunction', 'var'); fittingFunction = @PAL_Logistic; end
PF = fittingFunction;
searchGrid.alpha = -20:1:20; 
searchGrid.beta = 0:0.01:2; 
searchGrid.gamma = 0:.01:.1; 
searchGrid.lambda = 0:.01:.1;
paramsValues = PAL_PFML_Fit(xData(:), numChoices(:), numTrials(:), searchGrid, [1,1,1,1], PF,'lapseLimits',[0 1],'guessLimits', [0  1]);