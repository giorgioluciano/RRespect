% 
% Function: lcurve(input)
%
% Input: Gexp = n*1 vector [Gt],
%        Hgs  = guessed H,
%        t    = n*1 vector contains times,
%        s    = relaxation modes,
%
%        SmoothFac = Indirect way of controlling lambda_C, Set between -1 
%                   (lowest lambda explored) and 1 (highest lambda explored);
%                   When set to 0, using lambda_C determined from the
%                   L-curve,
%
% Output: lamC and 3 vectors of size npoints*1 contains a range of lambda, rho
% and eta. "Elbow"  = lamC is estimated using a heuristic.
%
% Plot the L-curve. Can call the routine "corner.m" to find the elbow again
%

function [lamC lam rho eta] = lcurve(Gexp, Hgs, t, s, SmoothFac)

	npoints  = 40;
	lam_min  = 1e-12;
	lam_max  = 1e+3;

	hlam    = (lam_max/lam_min)^(1/(npoints-1));
	lam     = zeros(npoints,1);
	lam     = lam_min * hlam.^[0:npoints-1]';

	eta     = zeros(npoints,1);
	rho     = zeros(npoints,1);

	%
	% This step can be "parfor"ed if you have multicore
	%
	
	for i = 1:length(lam)
		lambda  = lam(i);
		H       = LevenMarq(lambda, Gexp, Hgs, t, s);
		rho(i)  = norm((1 - kernel(H,t,s)./Gexp));
		eta(i)  = norm(diff(H,2));
	end

	% Strategy 1: From original version
	%~ i    = max(find(rho<1.01*min(rho)));
	
	% Another reasonable strategy
	%~ slope = abs(diff(eta)./diff(rho));
	%~ i = max(find(slope > 2.0)) - 1;

	%
	% 11/4/2015: New strategy to find corner!
	%
	er = rho/max(rho)+eta/max(eta);
	[ermin, i] = min(er);    

	lamC = lam(i);

	%
	% Dialling in the Smoothness Factor
	%
	if(SmoothFac > 0)
		lamC = exp(log(lamC) + SmoothFac*(log(lam_max) - log(lamC)));
	elseif(SmoothFac < 0)
		lamC = exp(log(lamC) + SmoothFac*(log(lamC) - log(lam_min)));
	end

end


