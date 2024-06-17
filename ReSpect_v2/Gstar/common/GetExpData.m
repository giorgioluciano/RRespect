%
% Function: GetExpData(input)
%
% Reads in the experimental data from the input file
%
% Input:  fname = name of file that contains G*(w) in 3 columns [w Gp Gpp]
%         Space it evenly on a log scale
%
% Output: A n*1 vector "w", and a 2n*1 vector Gexp = [Gp; Gpp]
%

function [w Gexp] = GetExpData(fname)

	data = load(fname);  % 3 columns, w - Gp - Gpp
	wo   = data(:,1);    % (n*1) vector
	Gpo  = data(:,2);
	Gppo = data(:,3);

	% any repeated frequency values
	
	[wo, i, j] = unique(wo);
	Gpo        = Gpo(i);
	Gppo       = Gppo(i);

	% Sanitize the input by spacing it out. Using linear interpolation

	w  =  transpose(logspace(log10(min(wo)),log10(max(wo)), numel(wo)));
	Gp =  interp1(wo, Gpo, w, 'linear','extrap');
	Gpp = interp1(wo, Gppo, w, 'linear','extrap');

	% Use super-smoother to clean it up further. May be optional
	Gp =  supsmu(w, Gp);
	Gpp = supsmu(w, Gpp);

	Gexp = [Gp;Gpp];     % Gp followed by Gpp (2n*1)

end
