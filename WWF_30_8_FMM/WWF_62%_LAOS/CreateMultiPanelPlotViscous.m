FilePath=['C:\Users\wahab\Dropbox\My PC (Arsenal)\ErrorCalc\Round2\WWF_30_8_FMM\WWF_62%_LAOS\'];



omegaV = [1, 10,50, 100];
AmpV = [0.05, 0.25, 0.5, 1.0, 2.0];



% Create plots
t = tiledlayout(size(AmpV, 2),size(omegaV, 2));

for j=size(AmpV, 2):-1:1
for i=1:size(omegaV, 2)
    
        

        
        
omega = omegaV(i); Amp = AmpV(j);
AmpP = Amp*100;

FileName = strcat('Exp_WGF62%_', num2str(AmpP), '%Om',num2str(omega), ...
    '.txt' );



q = load(strcat(FilePath, FileName) );


nexttile;
plot(q(:,3)/max(q(:,3)), q(:,2)/max(q(:,2) ), 'LineWidth',0.8 );
hold on
plot(q(:,3)/max(q(:,3)), q(:,4)/max(q(:,2) ), '--r', 'LineWidth',0.8 ) ;

title(['($\gamma$ = ',num2str(AmpP),'\%',', $\omega$ = ',num2str(omega),' rad/s)'],'interpreter','latex')

axis off


end

end


