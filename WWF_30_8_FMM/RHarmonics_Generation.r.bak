omegaV = [1, 3, 10, 30, 100];
AmpV = [0.001, 0.005, 0.01,0.03];
dim1 = size(omegaV,2);
dim2 = size(AmpV, 2);

for i= 1:dim1
    for j= 1:dim2
omega = omegaV(i); Amp = AmpV(j);
AmpP = Amp*100;
[harmonicsI, harmonicsR] =  ...
    KBKZmodelSimpsonLaunHarmonics([omega, Amp, 0.6658, 10.35, 278.9]);
FileName = strcat('NonlinearHarmonicsParameters', num2str(AmpP), ...
    '%Om',num2str(omega) , '.txt');

e1 = harmonicsI(1); e3 = - harmonicsI(2);
GpM = e1 - 3*e3;
GpL = e1 + e3;
v1 = harmonicsR(1)/omegaV(i);
v3 = harmonicsR(2)/omegaV(i);
etapM = v1 - 3*v3;
etapL = v1 + v3;
S = (GpL - GpM)/GpL;
T = (etapL - etapM)/etapL;



save(FileName, 'omega', 'Amp', 'e1', 'e3', 'GpM', 'GpL', 'v1', 'v3', 'etapM', ...
'etapL', 'S', 'T' ,  '-ascii');

disp(i);
disp(j);

    end
end








