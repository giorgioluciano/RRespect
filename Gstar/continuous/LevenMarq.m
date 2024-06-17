%
% Function LevenMarq(input)
%
% Purpose: Given a lambda, this function finds the H_lambda(s) that minimizes V(lambda)
%
%          V(lambda) := 1/n * ||Gst - kernel(H)||^2 + lambda/nl * ||L H||^2
%
% Input  : lambda = regularization parameter ,
%          Gst = experimental data,
%          H = guessed H,
%          w = n*1 vector contains frequencies,
%          s = relaxation modes,
%
% Output : H_lambda
% 
%          Uses the Levenberg-Marquardt method to solve the nonlinear
%          minimization problem. Following the algorithm in the source:
%
% "A Brief Description of the Levenberg-Marquardt Algorithm Implemened by levmar
%  Manolis I. A. Lourakis"
%

function [Hlambda] = LevenMarq(lambda, Gst, H, w, s)

	n   = length(w);
	ns  = length(s);
	nl  = ns - 2;
	hs  = s(2)/s(1);

	tau = 1e-3;          % scaling parameter 1
	nu  = 2;             % scaling parameter 2

	r   = zeros(2*n + nl,1);
	Jr  = zeros(2*n+nl,ns);
	i   = 0;
	pf  = sqrt(lambda);


	% L is a nl*ns tridiagonal matrix with 1
	% -2 and 1 on its diagonal.

	L  = diag(ones(ns-1,1),1)+diag(ones(ns-1,1),-1)+diag(-2*ones(ns,1));
	L  = L(2:nl+1,:);

	[r, Jr] = GetResidualJacobian(pf, L, Gst, H, w, s);
	mu      = tau * max(diag(Jr' * Jr));               % scaling parameter 3

	%
	% Doing the LM stuff
	%
	ContinueCriteria = 1;

	while(ContinueCriteria && i<5000)

		i=i+1;  % iteration count

		Delta = (Jr' * Jr + mu * eye(ns))\(-Jr' * r);   
		Hnew  = H + Delta;

		rnew = GetResidualJacobian(pf, L, Gst, Hnew, w, s);
		rho  = (norm(r)^2 - norm(rnew)^2)/(Delta' * (mu * Delta - Jr'*r));

		if(rho > 0)

			H       = Hnew;     
			[r, Jr] = GetResidualJacobian(pf, L, Gst, H, w, s);

			mu      = mu * max(1/3,1-(2*rho-1)^3);
			nu      = 2;

			% 
			% Do we need to stop now?
			% 

			if( ((norm(Delta) < 1e-6 ) &&  ...
			(norm(Jr'*r, inf) < 1e-6)     ||   ...
			(norm(Delta) < 1e-6 * norm(H))) )

			%~ if( (norm(Jr'*r, inf) < 1e-6) || ...
			%~ ((norm(Delta) < 1e-6 )     ||   ...
			%~ (norm(Delta) < 1e-6 * norm(H))) )

				ContinueCriteria = 0;

			end

		else

			mu      = mu * nu;
			nu      = 2  * nu;

		end

	end

	Hlambda = H;
  
end

%
% HELPER FUNCTION: Gets Residuals r and Jacobian J
% If called with only one output variable, only returns residual
%
function [r Jr] = GetResidualJacobian(pf, L, Gst, H, w, s)

	n   = length(w);
	ns  = length(s);
	nl  = ns - 2;

	r   = zeros(2*n + nl,1);
	Jr  = zeros(2*n+nl,ns);

	% 
	% Get the residual vector first
	% r = vector of size (2n+nl,1)
	%

	r(1:2*n)        = (1 - kernel(H,w,s)./Gst)/sqrt(n);  % the Gp and Gpp
	r(2*n+1:2*n+nl) = pf * diff(H,2)/sqrt(nl);            % second derivative

	%
	% Furnish the Jacobian Jr
	% (2n+nl)*ns matrix
	%

	if(nargout == 2)

		Kmatrix               =  (1./Gst) * ones(1,ns)/sqrt(n);
		Jr(1:2*n,1:ns)        = -kernelD(H,w,s) .* Kmatrix;
		Jr(2*n+1:2*n+nl,1:ns) = pf * L/sqrt(nl);

	end

end

