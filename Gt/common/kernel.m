%
% Function: kernel(input)
%
% outputs the n*1 dimensional vector K(H)(t) which is comparable to Gexp = Gt
% modifying kernel for unevenly spaced s_i
%
% Input: H = substituted CRS,
%        t = n*1 vector contains times
%        s = relaxation modes
%

function K = kernel(H, t, s)

	ns = length(s);
	hs = zeros(ns,1);

	%
	% integration uses trapezoidal rule
	%
	
	hs(1)      = 0.5*log(s(2)/s(1));
	hs(ns)     = 0.5*log(s(ns)/s(ns-1));
	hs(2:ns-1) = 0.5 * (log(s(3:ns))-log(s(1:ns-2)));

	[S, T]  = meshgrid(s,t);
	kern    = exp(-T./S);
	clear S T;

	K       = kern * (hs .* exp(H));

end
