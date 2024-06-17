function error2 = ErrorValCalcFun(para)

aSosval = para(1);
bSosval = para(2);

Lissajous_Plots_GenerationFun([aSosval, bSosval ]);
error2 = CreateMultiPanelPlotViscousFun();

end