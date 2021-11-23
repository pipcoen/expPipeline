function allRTPlotsNew
%%
opt.pltType = 'both';
opt.prmType = 'rea';
opt.revNorm = 1;
opt.siteLoc = 'contra';
opt.offset = 1;
fig3.scatterPlotsDifferenceWithLMETest(opt);

opt.prmType = 'tim';
fig3.scatterPlotsDifferenceWithLMETest(opt)
opt.offset = 0;

opt.prmType = 'rea';
opt.revNorm = 0;
opt.offset = 1;
fig3.scatterPlotsDifferenceWithLMETest(opt)

opt.prmType = 'tim';
opt.offset = 0;
fig3.scatterPlotsDifferenceWithLMETest(opt)


opt.prmType = 'rea';
opt.siteLoc = 'ipsi';
opt.offset = 1;
fig3.scatterPlotsDifferenceWithLMETest(opt)

opt.prmType = 'tim';
opt.offset = 0;
fig3.scatterPlotsDifferenceWithLMETest(opt)
end