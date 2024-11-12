function [taui, gi] =taugcal(p, Gp1, r, taumin, N)
taui = zeros(1, N);
H = zeros(1, N);
gi = zeros(1,N);
taui(1) = taumin;
for i=1:N-1
    taui(i+1) = taui(i)*r;
end
G1 = 2*Gp1*gamma(1+p)/p/pi*sin(p*pi/2);
for i = 1:N
    H(i) = G1*taui(i)^(-p)/gamma(p);
    gi(i) = H(i)*(r^(p/2) - r^(-p/2))/p;
end
end
    

function [Eabz] =Eabzcal(a, b, z, k)

Eabz = z.^k ./gamma(a .*k + b);

end
    

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
