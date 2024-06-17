%
% Function: InitializeH(input)
%
% Input:  Gexp = n*1 vector [Gt],
%         t = n*1 vector contains times,
%         s = relaxation modes,
%
% Output: H = guessed H
%

function H = InitializeH(Gexp,t,s)

	%
	% To guess spectrum, pick a negative Hgs and a large value of lambda to get a
	% solution that is most determined by the regularization, then use that as
	% the next guess. 
	%

	H       = -5.0 + sin(pi * s);

	lambda  = 1e0;
	Hlam    = LevenMarq(lambda, Gexp, H, t, s);

	%
	% Successively improve the initial guess until you have are reasonably good
	% guess for low lambda
	%

	lambda  = 1e-8;
	H       = LevenMarq(lambda, Gexp, Hlam, t, s);

end
