function allGLMModels(redo)
if ~exist('redo', 'var'); redo = 0; end
s = spatialAnalysis({'all'}, 'behavior', 0);
models2Fit = {'BiasOnly';'AudOnly';'VisOnly';'SimpEmp';'FullEmp';'SimpLog';'ReducedLog';'AudDom'};
for i = 1:length(models2Fit)
    savePath = ['D:\Dropbox (Personal)\XMatlabProg\Data4Plots\GLMFits2Behavior\CrossVal_' models2Fit{i} '.mat'];
    if exist(savePath, 'file') && ~redo; continue; end
    s.viewGLMFit(models2Fit{i}, 5);
    save(savePath, 's'); 
end

models2Fit = {'FullEmp'; 'BiasOnly'};
fileNames = {'FullEmpMaxPerformance'; 'BiasOnlyPerformance'};
for i = 1:length(models2Fit)
    savePath = ['D:\Dropbox (Personal)\XMatlabProg\Data4Plots\GLMFits2Behavior\' fileNames{i} '.mat'];
    if exist(savePath, 'file') && ~redo; continue; end
    s.viewGLMFit(models2Fit{i});
    save(savePath, 's'); 
end