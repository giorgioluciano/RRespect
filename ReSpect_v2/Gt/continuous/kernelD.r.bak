%
% Function: kernelD(input)
%
% outputs the n*ns dimensional vector DK(H)(t)
% approximates dK/dHj
%
% Input: H = substituted CRS,
%        t = n*1 vector of times,
%        s = relaxation modes,
%
% Output: DK = Jacobian of H
%

function DK = kernelD(H, t, s)

	ns         = length(s);
	hs         = zeros(ns,1);
	hs(1)      = 0.5*log(s(2)/s(1));
	hs(ns)     = 0.5*log(s(ns)/s(ns-1));
	hs(2:ns-1) = 0.5 * (log(s(3:ns))-log(s(1:ns-2)));

	n  = length(t);

	[S, T]  = meshgrid(s,t);
	kern    = exp(-T./S);
	clear S T;

	Hsuper  = ones(n,1)*(hs' .* exp(H'));       % A n*ns matrix with all the rows = H'
	DK      = kern  .* Hsuper;

end
