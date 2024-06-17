%
% Function: GetWeights(input)
%
% Finds the weight of "each" mode by taking a weighted average of its contribution
% to G(t)
%
% Input: H = CRS
%        t = n*1 vector contains times
%        s = relaxation modes
%
% Output: wt = weight of each mode
%

function wt = GetWeights(H, t, s)
  
	ns         = length(s);
	n          = length(t);

	hs         = zeros(ns,1);
	wt         = hs;
	hs(1)      = 0.5*log(s(2)/s(1));
	hs(ns)     = 0.5*log(s(ns)/s(ns-1));
	hs(2:ns-1) = 0.5 * (log(s(3:ns))-log(s(1:ns-2)));


	[S, T]  = meshgrid(s,t);
	kern    = exp(-T./S);
	clear S T;

	wij  =  kern * diag(hs .* exp(H));
	K    =  kern * (hs .* exp(H));

	for i = 1:n
		wij(i,:) = wij(i,:) ./ K(i);
	end

	for j = 1:ns
		wt(j) = sum(wij(:,j));
	end

end
