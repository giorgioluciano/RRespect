%
% Function: InitializeH(input)
%
% Input:  Gexp = 2n*1 vector [Gp; Gpp],
%         w = n*1 vector contains frequencies,
%         s = relaxation modes,
%
% Output: H = guessed H
%

function H = InitializeH(Gexp,w,s)

	%
	% To guess spectrum, pick a negative Hgs and a large value of lambda to get a
	% solution that is most determined by the regularization, then use that as
	% the next guess. 
	%

	%  H       = -5*ones(size(s));
	
	H       = -5.0 + sin(pi * s);
	lambda  = 1e0;
	Hlam    = LevenMarq(lambda, Gexp, H, w, s);

	%
	% Successively improve the initial guess until you have are reasonably good
	% guess for low lambda
	%
	
	lambda  = 1e-8;
	H       = LevenMarq(lambda, Gexp, Hlam, w, s);
   
end
