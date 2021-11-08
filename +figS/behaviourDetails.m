function behaviourDetails(behBlks)
%% This function plots the data panels for figure one of the ms
if ~exist('behBlks', 'var') || isempty(behBlks); behBlks = spatialAnalysis('all', 'behavior', 0, 1, ''); end
for i = [0.1 0.2]
    fig1.scatterAltPlots(behBlks, i);
    export_fig(['D:\OneDrive\Papers\Coen_2020\FigureParts\SupX_behScattersAlt_AudVer_' num2str(i)], '-pdf', '-painters'); 
    close
end
%%
if ~exist('behBlks', 'var') || isempty(behBlks); behBlks = spatialAnalysis('all', 'behavior', 0, 1, ''); end
for i = [0.1 0.2 0.4]
    fig1.scatterAltPlots(behBlks, i, 0);
    export_fig(['D:\OneDrive\Papers\Coen_2021\Revision\NewFigureParts\SupX_behScattersAlt_AudVer_RawTimes' num2str(i)], '-pdf', '-painters'); 
    close
end
end
