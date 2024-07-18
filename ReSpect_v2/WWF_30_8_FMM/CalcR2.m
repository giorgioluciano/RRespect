function R2 = CalcR2()



FilePath=[strcat(pwd,'\WWF_62%_LAOS\')];

FilePath2 = [strcat(pwd, '\')];




omegaV = [1, 10,50, 100];
AmpV = [0.05, 0.25, 0.5, 1.0, 2.0];


error = 0.0; 
ErrorAroundMean = 0.0;
R2 = 0.0;






for j=size(AmpV, 2):-1:1
for i=1:size(omegaV, 2)
    
        

        
        
omega = omegaV(i); Amp = AmpV(j);
AmpP = Amp*100;

FileName = strcat('Exp_WGF62%_', num2str(AmpP), '%Om',num2str(omega), ...
    '.txt' );

FileName2 = strcat('Theo_WGF62%_', num2str(AmpP), '%Om',num2str(omega), ...
    '.txt' );

q = load(strcat(FilePath, FileName) );
qtheo = load(strcat(FilePath2, FileName2) );



xinterp = [1:1:size(q, 1)];
theoElStrvals = (interp1([0:1:size(qtheo)-1], qtheo(:, 5) ./ ...
    max(qtheo(:, 2) ), (xinterp-1).*(size(qtheo,1)-1)./size(q,1) ) )';
error = error+ sum((theoElStrvals-q(:, 5)./ max(q(:,2) )).^2 )/size(q,1);

MeanEl = mean(q(:, 5)./ max(q(:,2) ) );
ErrorAroundMean = ErrorAroundMean + sum((MeanEl ...
    -q(:, 5)./ max(q(:,2) )).^2 )/size(q,1);

R2 = R2 + 0.5*(1-error/ErrorAroundMean);


theoVisStrvals = (interp1([0:1:size(qtheo)-1], qtheo(:, 4) ./ ...
    max(qtheo(:, 2) ), (xinterp-1).*(size(qtheo,1)-1)./size(q,1) ) )';

error = error+ sum((theoVisStrvals-q(:, 4)./ max(q(:,2) ) ).^2  )/size(q,1);

MeanVis = mean(q(:, 4)./ max(q(:,2) ));

ErrorAroundMean = ErrorAroundMean + sum((MeanVis ...
    -q(:, 4)./ max(q(:,2) ) ).^2  )/size(q,1);

R2 = R2 + 0.5*(1-error/ErrorAroundMean);


    end
end


% create plots elastic



error = error/size(AmpV, 2)/size(omegaV, 2)/2;

R2 = R2/size(AmpV, 2)/size(omegaV, 2);

% disp(['error = ', num2str(sqrt(error) )] );


end
