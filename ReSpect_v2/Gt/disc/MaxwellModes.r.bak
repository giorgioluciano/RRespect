%
% Function: MaxwellModes(input)
%
% Solves the linear least squares problem to obtain the DRS
%
% Input: z = points distributed according to the density,
%        t  = n*1 vector contains times,
%        Gt = n*1 vector contains G(t),
%
%        Prune = Avoid modes with -ve weights (=1), or don't care (0) 
%
% Output: g, tau = spectrum  (array)
%         error = relative error between the input data and the G(t) inferred from the DRS
%         condKp = condition number
%

function [g tau error condKp] = MaxwellModes(z, t, Gt, prune)

	N      = length(z);
	tau    = exp(z);
	n      = length(t);
	Gexp   = Gt;

	%
	% Prune small -ve weights g(i)
	%

	if(nargin < 4)
		prune = 0;
	end

	if(~prune)

		[g error condKp] = LLS(t,tau,Gexp);

	else

		tau1 = tau;
		[g1 error1 condKp1] = LLS(t,tau1,Gexp);

		%
		% Run through negative "g"s and mark for
		% deletion if the weight is small, and
		% repeat above calculation
		%
		ineg       = find(g1<0);
		gneg       = g1(ineg);

		tau2       = tau1;
		tau2(ineg) = [];

		[g2 error2 condKp2] = LLS(t,tau2,Gexp);

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

function [g error condKp] = LLS(t,tau,Gexp)

	n       = length(Gexp);
	[S, T]  = meshgrid(tau,t);
	K       = exp(-T./S);
	clear S T;

	% gets (Gt/GtE - 1)^2, instead of  (Gt -  GtE)^2

	Kp      = diag(1./Gexp) * K;  
	condKp  = cond(Kp);
	g       = Kp\ones(size(Gexp));

	clear Kp;

	GtM   = K * g;
	error = sum((GtM./Gexp - 1).^2);

end
