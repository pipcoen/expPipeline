function exampleInactModelFitsSignificantTestShuffle
%%
figure;
axHeight = 250;
axWidth = 250;
nCols = 3;
nRows = 2;
figHeight = nRows*axHeight;
figWidth = nCols*axWidth;

axesGap = [50/figHeight 50/figWidth];
botTopMarg = [40, 40]/figHeight;
lftRgtMarg = [40, 40]/figWidth;
set(gcf, 'position', get(gcf, 'position').*[1 0 0 0] + [0 200 figWidth, figHeight]);
%%
load inactCompResults210802.mat;
logLikR = zeros(5,5);
nReg = 1:5;
selfTestlogLikR = zeros(4,1);
for i = nReg
    for j = nReg
        if i == j
            compIdx = nReg(~ismember(nReg,j));
            for k = 1:4
                idx = ismember(trainTestGroups, [i,compIdx(k)], 'rows');
                selfTestlogLikR(k,1) = mean(testGrp1LogLik{idx}(1:normEstRepeats));
            end
            logLikR(i,j) = mean(selfTestlogLikR);
            continue; 
        end
        idx = ismember(trainTestGroups, [i,j], 'rows');
        sIdx = normEstRepeats+1:length(testGrp1LogLik{idx});
        
        LLRSelf = mean(testGrp1LogLik{idx}(1:normEstRepeats));
        LLRTest = mean(testGrp2LogLik{idx}(1:normEstRepeats));       
        shuffTest = sort([(LLRSelf-LLRTest); testGrp1LogLik{idx}(sIdx)-testGrp2LogLik{idx}(sIdx)]);
        
        pVal(i,j) = find(shuffTest == LLRSelf-LLRTest)/length(sIdx)*5;        
        logLikR(i,j) = LLRTest;
    end
end

%%
plt.tightSubplot(nRows,nCols,[1 4],axesGap,botTopMarg,lftRgtMarg); cla;
lRef = {'F'; 'V'; 'L'; 'S'; 'C'};
xlim([0.5 5.5])
ylim([-1 -0.55])
cCol = {'r'; 'b'; 'g'; 'm'; 'k'};
for i = nReg
    for j = nReg
%         text(i,logLikR(i,j)*-1,lRef{j}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        plot(i,logLikR(i,j)*-1, ['.' cCol{j}])
        hold on;
    end
end
set(gca, 'XTick', 1:5, 'XTickLabel', {'Frontal'; 'Visual'; 'Lateral'; 'Somat.'; 'Control'});
set(gca, 'position', get(gca, 'position').*[2 2.2 1 0.9]);
xlabel('Test Region')
ylabel('Logliklihood')
xtickangle(45)
export_fig('D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\sigTest4InactivationSitesShuffle', '-pdf', '-painters');
%%
% 
% imagesc(logLikR*-1)
% axis square
% set(gca, 'XTick', 1:5, 'XTickLabels', posNames)
% set(gca, 'YTick', 1:5, 'YTickLabels', posNames)
% ylabel('Trained site');
% xlabel('Test site');
% box off
% colormap('gray'); colorbar; caxis([-0.1 0])
% 
% 
% plt.tightSubplot(nRows,nCols,4:5,axesGap,botTopMarg,lftRgtMarg); cla;
% set(gca, 'position', get(gca, 'position').*[1 1.5 1 1])
% plotDat = mean(logLikR*-1,3);
% plotDat = mean(cat(3,plotDat, plotDat'), 3);
% plotDat = triu(plotDat);
% 
% imagesc(plotDat)
% axis square
% set(gca, 'XTick', 1:5, 'XTickLabels', posNames)
% set(gca, 'YTick', 1:5, 'YTickLabels', posNames)
% ylabel('Inactivation site 1');
% xlabel('Inactivation site 2');
% box off
% colormap('gray'); colorbar; caxis([-0.1 0])

%%
