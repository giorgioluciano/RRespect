%
% Function: GetExpData(input)
%
% Reads in the experimental data from the input file
%
% Input:  fname = name of file that contains G(t) in 2 columns [t Gt]
%
% Output: A n*1 vector "t", and a n*1 vector Gt
%

function [t Gt] = GetExpData(fname)

	data = load(fname);  % 2 columns, t - Gt
	to   = data(:,1);    % (n*1) vector
	Gto  = data(:,2);

	% any repeated "time" values
	[to, i, j] = unique(to);
	Gto        = Gto(i);

	% Sanitize the input by spacing it out. Using linear interpolation

	t  =  transpose(logspace(log10(min(to)),log10(max(to)), numel(to)));
	Gt =  interp1(to, Gto, t, 'linear','extrap');

end
