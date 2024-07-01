%
% Function PlotMaxwellModes(input)
%
% Plots and compares experimental G(t) with that obtained from the
% corresponding corresponding DRS
%
% Input: g, tau = spectrum (array)
%        t  = n*1 vector contains times,
%        Gt = n*1 vector contains G(t),
%
function PlotMaxwellModes(g, tau, t, Gt);

	N  = length(g);

	[S, T]  = meshgrid(tau,t);
	K       = exp(-T./S);
	GtM     = K * g;

	loglog(t,Gt,'bo',t,GtM,'k-','LineWidth',2)

end
