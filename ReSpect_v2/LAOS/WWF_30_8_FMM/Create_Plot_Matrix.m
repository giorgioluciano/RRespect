Amps = [0.05, 0.1, 0.3, 1, 2];
omegas = [1, 10, 100];
dim1 = size(Amps, 2);
dim2 = size(omegas, 2);
figs = zeros(dim1, dim2);


for i = 1:dim1
    for j=1:dim2
        AmpP = Amps(i)*100; 
        openfig(strcat('ElasticSt', num2str(AmpP), ...
         '%Om',num2str(omegas(j) ),'.fig' ));
     
     % Load saved figures
        figs(dim1-i+1 , j) = hgload(strcat('ElasticSt', num2str(AmpP), ...
         '%Om',num2str(omegas(j) ),'.fig' ));
     end
end

h = zeros(1, dim1*dim2);
figure
for i = 1:dim1
    for j=1:dim2
% Prepare subplots

h((i-1)*dim2+j)=subplot(dim1,dim2,(i-1)*dim2+j );
% Paste figures on the subplots
copyobj(allchild(get(figs(i, j),'CurrentAxes')),h((i-1)*dim2+j ));
axis off
pbaspect([1 1 1]);
axis([-1 1 -1 1])
    end
end


