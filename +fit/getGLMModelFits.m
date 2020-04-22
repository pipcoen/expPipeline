function getGLMModelFits(models2fit, crossVal, maxAndMin)
modelSet = {'biasOnly';'visOnly';'audOnly';'simpLog';'simpLogSplitV';'simpLogSplitA';'simpLogSplitVSplitA';...
    'fullEmp';'simpEmp';'visOnlyEmp';'audOnlyEmp'};  

sStart = spatialAnalysis('all', 'behavior', 0, 1);
if ~exist('models2fit', 'var'); models2fit = modelSet; end
if ~exist('crossVal', 'var'); crossVal = 5; end
if ~exist('maxAndMin', 'var'); maxAndMin = 1; end

expListDir = prc.pathFinder('expList');
saveDir = [expListDir(1) ':\Dropbox (Personal)\XMatlabProg\GitHub\expPipeline\data4Plots\GLMFits2Behavior\'];

for i = 1:length(models2fit)
    fprintf('Currently fitting model %s \n', models2fit{i});
    s = sStart;
    s.viewGLMFit(models2fit{i}, crossVal);
    save([saveDir models2fit{i} '_Cross' num2str(crossVal) '.mat'], 's');
    close;
end

if maxAndMin
    s = sStart;
    s.viewGLMFit('biasOnly');
    save([saveDir 'BiasOnlyPerformance.mat']);
    close;
    
    s = sStart;
    s.viewGLMFit('fullEmp');
    save([saveDir 'FullEmpMaxPerformance.mat']);
    close;
end

end
