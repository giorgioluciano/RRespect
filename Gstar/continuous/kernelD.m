%
% Function: kernelD(input)
%
% outputs the 2n*ns dimensional vector DK(H)(w)
% approximates dK/dHj
%
% Input: H = substituted CRS,
%        w = n*1 vector contains frequencies,
%        s = relaxation modes,
%
% Output: DK = Jacobian of H
%

function DK = kernelD(H, w, s)

	ns         = length(s);
	hs         = zeros(ns,1);
	hs(1)      = 0.5*log(s(2)/s(1));
	hs(ns)     = 0.5*log(s(ns)/s(ns-1));
	hs(2:ns-1) = 0.5 * (log(s(3:ns))-log(s(1:ns-2)));

	n  = length(w);

	[S, W]  = meshgrid(s,w);
	ws      = S.*W;
	clear S W;

	ws2     = ws.^2;
	Hsuper  = ones(2*n,1)*(hs' .* exp(H'));       % A 2n*ns matrix with all the rows = H'

	DK      = [ (ws2./(1+ws2)) ;(ws./(1+ws2)) ] .* Hsuper;

end
