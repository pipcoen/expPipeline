function allRTPlotsComplete
%%
opt.prmType = 'rea';
opt.pltType = 1;
fig3.scatterPlotsDifferenceWithLMETestComplete(opt);
%%
opt.prmType = 'tim';
opt.pltType = 1;
fig3.scatterPlotsDifferenceWithLMETestComplete(opt)
%%
opt.prmType = 'rea';
opt.pltType = 2;
fig3.scatterPlotsDifferenceWithLMETestComplete(opt);
%%
opt.prmType = 'tim';
opt.pltType = 2;
fig3.scatterPlotsDifferenceWithLMETestComplete(opt);
end