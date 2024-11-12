function Lissajous_Plots_GenerationFun(para)

aSosval = para(1);
bSosval = para(2);
omegaV = [1, 10,50, 100];
AmpV = [0.05, 0.25, 0.5, 1.0, 2.0];

dim = size(omegaV, 2)*size(AmpV, 2);

parfor ij=1:dim
   i = floor((ij-1)/size(AmpV, 2) )+1;
   j = ij-(i-1)*size(AmpV,2);
   disp(ij);
omega = omegaV(i); Amp = AmpV(j);
AmpP = Amp*100;
[time2, strain, strainRt , tau12t2, ElStr, VisStr] =  ...
    KBKZmodelSimpsonSoskeyLissajous([omega, Amp, aSosval, bSosval]);
disp(ij);
tau12t2 = ElStr+VisStr;
figure('Name',strcat('ElasticSt', num2str(AmpP), '%Om',num2str(omega) ) );
plot(strain/max(strain), tau12t2/(max(tau12t2)) ,'k', 'LineWidth',1);
hold on
plot(strain/max(strain), ElStr/(max(tau12t2)) ,'--r', 'LineWidth',0.8);
axis off
hold off
saveas(gcf,strcat('ElasticSt', num2str(AmpP), '%Om',num2str(omega) ),'fig')


figure('Name',strcat('ViscousSt', num2str(AmpP), '%Om',num2str(omega) ) );

data = [strain; tau12t2; strainRt; VisStr; ElStr]';

dlmwrite(strcat('Theo_WGF62%_', num2str(AmpP), '%Om',num2str(omega), ...
    '.txt' ), data, 'delimiter', '\t', 'newline', 'pc');

plot(strainRt/max(strainRt), tau12t2/max(tau12t2),'k', 'LineWidth',1 );
hold on
plot(strainRt/max(strainRt) , VisStr/(max(tau12t2)),'--r', 'LineWidth',0.8 );
axis off
hold off
saveas(gcf,strcat('ViscousSt', num2str(AmpP), '%Om',num2str(omega) ),'fig');

   
end

end