%
% Function: GetWeights(input)
%
% Finds the weight of "each" mode by taking a weighted average of its contribution
% to Gp and Gpp
%
% Input: H = substituted CRS
%        w = n*1 vector contains frequencies
%        s = relaxation modes
%
% Output: wt = weight of each mode
%

function wt = GetWeights(H,w,s)
  
	ns         = length(s);
	n          = length(w);

	hs         = zeros(ns,1);
	wt         = hs;
	hs(1)      = 0.5*log(s(2)/s(1));
	hs(ns)     = 0.5*log(s(ns)/s(ns-1));
	hs(2:ns-1) = 0.5 * (log(s(3:ns))-log(s(1:ns-2)));


	[S, W]  = meshgrid(s,w);
	ws      = S.*W;
	ws2     = ws.^2;

	clear X Y;

	wij  =  [ (ws2./(1+ws2)) ;(ws./(1+ws2)) ] * diag(hs .* exp(H));
	K    =  [ (ws2./(1+ws2)) ;(ws./(1+ws2)) ] * (hs .* exp(H));

	for i = 1:2*n
		wij(i,:) = wij(i,:) ./ K(i);
	end

	for j = 1:ns
		wt(j) = sum(wij(:,j));
	end

end
