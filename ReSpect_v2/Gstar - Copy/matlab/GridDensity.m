%
%  PROGRAM: GridDensity(input)
%
%	Takes in a PDF or density function, and spits out a bunch of points in
%       accordance with the PDF
%
%  Input:
%       x  = vector of points. It need *not* be equispaced,
%       px = vector of same size as x: probability distribution or
%            density function. It need not be normalized but has to be positive.
%  	    N  = Number of points >= 3. The end points of "x" are included
%  	     necessarily,
%       Pt = Optional argument. If present, then some plotting.
% 
%  Output:
%       z  = Points distributed according to the density
%       hz = width of the "intervals" - useful to apportion domain to points
%            if you are doing quadrature with the results, for example.
%
%  (c) Sachin Shanbhag, March 5, 2012
%

function [z h] = GridDensity(x, px, N, Pt)

	npts = 100;                              % can potentially change
	xi   = linspace(min(x),max(x),npts)';    % reinterpolate on equi-spaced axis
	pint = interp1(x,px,xi,'spline');        % smoothen using splines
	ci   = cumtrapz(xi,pint);                
	pint = pint/ci(npts);
	ci   = ci/ci(npts);                      % normalize ci

	alfa = 1/(N-1);                          % alfa/2 + (N-1)*alfa + alfa/2
	zij  = zeros(N,1);                       % quadrature interval end marker
	z    = zeros(N,1);                       % quadrature point

	z(1)  = min(x);  
	z(N)  = max(x); 

	%
	% ci(Z_j,j+1) = (j - 0.5) * alfa
	%
	
	beta  = [0.5:1:N-1.5]'*alfa;
	zij   = [z(1); interp1(ci,xi,beta,'spline'); z(N)];
	h        = diff(zij);
	clear beta;
	
	%
	% Quadrature points are not the centroids, but rather the center of masses
	% of the quadrature intervals
	%
	beta     = [1:1:N-2]'*alfa;
	z(2:N-1) = interp1(ci,xi,beta,'spline');

%
% Some plotting if required
%
	if(nargin>3)

		subplot(2,1,1)

		plot(xi,ci,'b-','LineWidth',2,xi,pint,'r-','LineWidth',2);
		axis([min(xi),max(xi)]);
		legend('CDF','PDF')
		title('Visualization')
		xlabel('x');
		ylabel('CDF(x)/PDF(x)')

		subplot(2,1,2)
		
		plot(z,0.5*ones(size(z)),'ro');
		axis([min(xi),max(xi)]);
		xlabel('z');
		hold on;
		for i = 2:N
			X = [zij(i); zij(i)];
			Y = [0; 1];
			plot(X,Y,'b')
		end
		hold off;

	end

end


  
  
