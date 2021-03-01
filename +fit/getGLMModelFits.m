function getGLMModelFits(models2fit, crossVal, maxAndMin)
modelSet = {'biasOnly';'visOnly';'audOnly';'simpLog';'simpLogSplitV';'simpLogSplitA';'simpLogSplitVSplitA';...
    'fullEmp';'simpEmp';'visOnlyEmp';'audOnlyEmp'};  

allBlks = spatialAnalysis('all', 'behavior', 1, 0);
allBlks.blks = prc.filtBlock(allBlks.blks, allBlks.blks.tri.stim.visContrast ~= 0.06);

sStart = spatialAnalysis('all', 'behavior', 0, 1);
for i = 1:length(sStart.blks)
    sStart.blks(i) = prc.filtBlock(sStart.blks(i), sStart.blks(i).tri.stim.visContrast ~= 0.06);
end

sStart.blks(end+1) = allBlks.blks;
allBlks.blks.exp.subject{1} = 'Combined';
if ~exist('models2fit', 'var'); models2fit = modelSet; end
if ~exist('crossVal', 'var'); crossVal = 5; end
if ~exist('maxAndMin', 'var'); maxAndMin = 1; end

expListDir = prc.pathFinder('expList');
saveDir = [expListDir(1) ':\Dropbox (Personal)\XMatlabProg\GitHub\expPipeline\data4Plots\GLMFits2Behavior\'];

for i = 1:length(models2fit)
    fprintf('Currently fitting model %s \n', models2fit{i});
    s = sStart;
    s.viewGLMFits(models2fit{i}, crossVal);
    save([saveDir models2fit{i} '_Cross' num2str(crossVal) '.mat'], 's');
    close;
end

if maxAndMin
    s = sStart;
    s.viewGLMFits('biasOnly');
    save([saveDir 'BiasOnlyPerformance.mat']);
    close;
    
    s = sStart;
    s.viewGLMFits('fullEmp');
    save([saveDir 'FullEmpMaxPerformance.mat']);
    close;
end

end
