aSoso = 49.1698; 
aSosf = 1400.0;

bSoso = 1.00261;
bSosf = 1.90;

intervals = 5;


aSosVals = aSoso:(aSosf-aSoso)/intervals :aSosf; 
bSosVals = bSoso:(bSosf-bSoso)/intervals :bSosf;  



dim1 = size(aSosVals, 2);
dim2 = size(bSosVals, 2);

ErrorMat = zeros(dim1, dim2);

for i = 1: dim1
    for j = 1: dim2
        
        Lissajous_Plots_GenerationFun([aSosVals(i), bSosVals(j) ]);
        ErrorMat(i,j) = CreateMultiPanelPlotViscousFun();
        
    end
end