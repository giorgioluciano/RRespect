Amps = [0.001, 0.005, 0.01,0.03, 0.05, 0.1, 0.3,0.6, 1.0,1.5, 2.0];
omegas = [1,3, 10, 30, 100];
dim1 = size(Amps, 2);
dim2 = size(omegas, 2);


dim3 = 55;
dim4 = 12;

Alldata = zeros(dim3, dim4);





for i = 1:dim1
    for j=1:dim2
        
        AmpP = Amps(i)*100; 
        omega = omegas(j);
        Table = importdata(strcat('NonlinearHarmonicsParameters', ...
            num2str(AmpP), '%Om',num2str(omega) , '.txt'));
        
        Alldata((j-1)*dim1 + i, :) = Table;
     
     
     end
end

xdata = zeros(1, dim1);
ydata = zeros(1, dim1);
legends = cell(1, dim2);

figure('Name', 'e3');
hold on;
for j = 1:2:dim2
    for i=1:dim1
        xdata(i) = Alldata( (j-1)*dim1 + i, 2 )*100;
        ydata(i) = Alldata( (j-1)*dim1 + i, 4 );
        omega = Alldata((j-1)*dim1 + i, 1 );
    end
    xi = linspace(min(xdata), max(xdata), 150);  % Evenly-Spaced Interpolation Vector
    yi = interp1(xdata, ydata, xi, 'spline', 'extrap');
    plot(xi, yi, 'LineWidth',1.5);    
    % legends(j) = strcat('Omega = ',num2str(omega) ) ;
end


figure('Name', 'v3');
hold on;
for j = 1:2:dim2
    for i=1:dim1
        xdata(i) = Alldata( (j-1)*dim1 + i, 2 )*100;
        ydata(i) = Alldata( (j-1)*dim1 + i, 8 );
        omega = Alldata((j-1)*dim1 + i, 1 );
    end
    xi = linspace(min(xdata), max(xdata), 150);  % Evenly-Spaced Interpolation Vector
    yi = interp1(xdata, ydata, xi, 'spline', 'extrap');
    plot(xi, yi, 'LineWidth',1.5);    
    % legends(j) = strcat('Omega = ',num2str(omega) ) ;
end


figure('Name', 'S data');
hold on;
for j = 1:2:dim2
    for i=1:dim1
        xdata(i) = Alldata( (j-1)*dim1 + i, 2 )*100;
        ydata(i) = Alldata( (j-1)*dim1 + i, 11 );
        omega = Alldata((j-1)*dim1 + i, 1 );
    end
    xi = linspace(min(xdata), max(xdata), 150);  % Evenly-Spaced Interpolation Vector
    yi = interp1(xdata, ydata, xi, 'spline', 'extrap');
    plot(xi, yi, 'LineWidth',1.5);    
    % legends(j) = strcat('Omega = ',num2str(omega) ) ;
end


figure('Name', 'T data');
hold on;
for j = 1:2:dim2
    for i=1:dim1
        xdata(i) = Alldata( (j-1)*dim1 + i, 2 )*100;
        ydata(i) = Alldata( (j-1)*dim1 + i, 12 );
        omega = Alldata((j-1)*dim1 + i, 1 );
    end
    xi = linspace(min(xdata), max(xdata), 150);  % Evenly-Spaced Interpolation Vector
    yi = interp1(xdata, ydata, xi, 'spline', 'extrap');
    plot(xi, yi, 'LineWidth',1.5);    
    % legends(j) = strcat('Omega = ',num2str(omega) ) ;
end
        


    


