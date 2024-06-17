%
% Function: MaxwellModes(input)
%
% Solves the linear least squares problem to obtain the DRS
%
% Input: z = points distributed according to the density,
%        w = n*1 vector contains frequencies,
%        Gp, Gpp
%        Prune = Avoid modes with -ve weights (=1), or don't care (0) 
%
% Output: g, tau = spectrum 
%         error = relative error between the input data and the G*(w) inferred from the DRS
%         condKp = condition number
%

function [g tau error condKp] = MaxwellModes(z, w, Gp, Gpp, prune)

	N      = length(z);
	tau    = exp(z);
	n      = length(w);
	Gexp   = [Gp; Gpp];
	
	%
	% Prune small -ve weights g(i)
	%
	if(nargin < 4)
		prune = 0;
	end

	if(~prune)

		[g error condKp] = LLS(w,tau,Gexp);

	else

		tau1 = tau;
		[g1 error1 condKp1] = LLS(w,tau1,Gexp);

		%
		% Run through negative "g"s and mark for
		% deletion if the weight is small, and
		% repeat above calculation
		%
		ineg       = find(g1<0);
		gneg       = g1(ineg);

		tau2       = tau1;
		tau2(ineg) = [];

		[g2 error2 condKp2] = LLS(w,tau2,Gexp);

		if(condKp2 < condKp1 && (error2/error1 - 1) < 0.05)
			error  = error2;
			condKp = condKp2;
			g      = g2;
			tau    = tau2;
		else
			error  = error1;
			condKp = condKp1;
			g      = g1;
			tau    = tau1;
		end

	end

end

%
% Subfunction which does the actual LLS problem
%

function [g error condKp] = LLS(w,tau,Gexp)

	n       = length(Gexp)/2;
	[X, Y]  = meshgrid(tau,w);
	ws      = X.*Y;
	ws2     = ws.^2;
	clear X Y;

	K       = [ (ws2./(1+ws2)) ;(ws./(1+ws2)) ];
	Kp      = diag(1./Gexp) * K;  	% gets (Gp/GpE - 1)^2, instead of  (Gp -  GpE)^2


	clear K;

	condKp  = cond(Kp);
	g       = Kp\ones(size(Gexp));

	clear Kp;

	GpM   = (ws2./(1+ws2)) * g;
	GppM  = (ws./(1+ws2))  * g;
	error = sum((GpM./Gexp(1:n) - 1).^2 + (GppM./Gexp(n+1:2*n) - 1).^2);

end
