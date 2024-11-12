FilePath=[strcat(pwd,'\BF62% TEXT 24-5-2021 Batch2\')];

FilePath2 = [strcat(pwd, '\')];

omegaV = [1, 10,50, 100];
AmpV = [0.05, 0.25, 0.5, 1.0, 2.0];



% Create plots
t = tiledlayout(size(AmpV, 2),size(omegaV, 2));

for j=size(AmpV, 2):-1:1
for i=1:size(omegaV, 2)
    
        

        
        
omega = omegaV(i); Amp = AmpV(j);
AmpP = Amp*100;

FileName = strcat('Exp_BF62%_', num2str(AmpP), '%Om',num2str(omega), ...
    '.txt' );

FileName2 = strcat('Theo_BF62%_', num2str(AmpP), '%Om',num2str(omega), ...
    '.txt' );

q = load(strcat(FilePath, FileName) );
qtheo = load(strcat(FilePath2, FileName2) );

nexttile;
plot(q(:,1)/max(q(:,1)), q(:,2)/max(q(:,2) ), 'LineWidth',0.8 );
hold on
plot(q(:,1)/max(q(:,1)), q(:,5)/max(q(:,2) ), '--r', 'LineWidth',0.8 ) ;

plot(qtheo(:,1)/max(qtheo(:,1)), qtheo(:,2)/max(qtheo(:,2) ) ... 
    ,'g', 'LineWidth',0.8 );
hold on
plot(qtheo(:,1)/max(qtheo(:,1)), qtheo(:,5)/max(qtheo(:,2) ) ...
    , '--g', 'LineWidth',0.8 ) ;

axis off


    end
end
