%
% Function: kernel(input)
%
% outputs the 2n*1 dimensional vector K(H)(w) which is comparable to Gexp
% modifying kernel for unevenly spaced s_i
%
% Input: H = substituted CRS,
%        w = n*1 vector contains frequencies,
%        s = relaxation modes
%

function K = kernel(H, w, s)

	ns = length(s);
	hs = zeros(ns,1);

	%
	% uses trapezoidal rule for integration
	%
	hs(1)  = 0.5*log(s(2)/s(1));
	hs(ns) = 0.5*log(s(ns)/s(ns-1));

	hs(2:ns-1) = 0.5 * (log(s(3:ns))-log(s(1:ns-2)));

	[S, W]  = meshgrid(s,w);
	ws      = S.*W;
	ws2     = ws.^2;
	clear S W;

	K       = [ (ws2./(1+ws2)) ;(ws./(1+ws2)) ] * (hs .* exp(H));

end
