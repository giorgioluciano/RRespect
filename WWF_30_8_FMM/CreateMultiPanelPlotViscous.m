FilePath=['C:\Users\wahab\Dropbox\Esam Thesis\' ... 
    'Rheology of dough\BF62%\BF62%-EXP-TXT\'];

FilePath2 = [strcat(pwd, '\')];


omegaV = [1, 10,50, 100];
AmpV = [0.05, 0.25, 0.5, 1.0, 2.0];


error = 0.0; 
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
plot(q(:,3)/max(q(:,3)), q(:,2)/max(q(:,2) ), 'LineWidth',0.8 );
hold on
plot(q(:,3)/max(q(:,3)), q(:,4)/max(q(:,2) ), '--r', 'LineWidth',0.8 ) ;

plot(qtheo(:,3)/max(qtheo(:,3)), qtheo(:,2)/max(qtheo(:,2) ) ... 
    ,'g', 'LineWidth',0.8 );
hold on
plot(qtheo(:,3)/max(qtheo(:,3)), qtheo(:,4)/max(qtheo(:,2) ) ...
    , '--g', 'LineWidth',0.8 ) ;

axis off


xinterp = [1:1:size(q, 1)];
theoElStrvals = (interp1([0:1:size(qtheo)-1], qtheo(:, 5) ./ ...
    max(qtheo(:, 2) ), (xinterp-1).*(size(qtheo,1)-1)./size(q,1) ) )';
error = error+ sum((theoElStrvals-q(:, 5)./ max(q(:,2) )).^2 )/size(q,1);

theoVisStrvals = (interp1([0:1:size(qtheo)-1], qtheo(:, 4) ./ ...
    max(qtheo(:, 2) ), (xinterp-1).*(size(qtheo,1)-1)./size(q,1) ) )';

error = error+ sum((theoVisStrvals-q(:, 4)./ max(q(:,2) ) ).^2  )/size(q,1);


    end
end
error = error/size(AmpV, 2)/size(omegaV, 2)/2;
disp(['error = ', num2str(sqrt(error) )] );
