%
% Function: discSpec(par)
%
% Uses the continuous relaxation spectrum extracted using contSpec()
% to determine an approximate discrete approximation.
%
% Input: Communicated by the datastructure "par"
%
%        fGstFile = name of file that contains G(t) in 2 columns [t Gt]
%                   default: 'Gt.dat' is assumed.
%        verbose  = 1, then prints onscreen messages, and prints datafiles
%        plotting = 1, then plots to stdio.
%
%        prune    = 1, then tries to kill modes with -ve g(i)
%  
%        Nopt     = optional argument, if you want to get a spectrum for a
%                   specified number of modes. If absent it will use some
%                   heuristic algorithm to figure out an optimum.
%        
% In addition par.BaseDistWt and par.condWt control the blending of the
% flat profile, and the weight of the "condition number" in determining
% the optimal Nopt. Both are set to 0.5 by default.
%
%
% Output: Nopt    = optimum number of discrete modes
%         [g tau] = spectrum
%         error   = error norm of the discrete fit
%        
%         dmodes.dat : Prints the [g tau] for the particular Nopt
%         Nopt.dat   : If Nopt not supplied, then optimum [N error(N) cond(N)]
%         Gfitd.dat  : The discrete G(t) for Nopt [t Gt] 
%

function [Nopt g tau error] =  discSpec(par)

	%
    % Add appropriate subdirectories to search path
    %
	addpath('./disc:./common')

	if nargin == 0
		par = SetParameters();  % Load in global settings
	end

	if(par.verbose)
		fprintf('\n(*) Start\n(*) Loading Data Files: ...');
	end

	[t, Gt, s, H] = ReadData(par.GstFile,'output/H.dat');

	n   = length(t);
	ns  = length(s);

	%
	% Find the distribution of nodes you need
	%
	wt  = GetWeights(H,t,s);
	wt  = wt/trapz(log(s),wt);
	wt  = (1 - par.BaseDistWt) * wt + (par.BaseDistWt*mean(wt))*ones(size(wt));

	% 
	% Try different N: number of Maxwell modes
	%

	Nmax  = min(floor(3*log(max(t)/min(t))),n/4); % maximum N
	Nmin  = max(floor(0.5*log10(max(t)/min(t))),3);   % minimum N

	npts  = Nmax - Nmin + 1;
	Nv    = Nmin:1:Nmax;
	ev    = zeros(size(Nv));
	condN = ev;

	for i = 1:npts
		N = Nv(i);
		[z hz]                 = GridDensity(log(s),wt,N);            % Select "tau" Points
		[g tau ev(i) condN(i)] = MaxwellModes(z, t, Gt, par.prune);   % Get g_i
	end

	%
	% Use supplied number of modes or
	%
	if(par.Nopt > 0)

		Nopt = par.Nopt;

		if(par.verbose)
			fprintf('\n(*) Using %d number of discrete modes\n',Nopt);
		end

	else      

		emin     = min(ev);
		condNmin = min(condN);

		cost     = (1 - par.condWt) * (ev - emin).^2 + ...
					par.condWt  * (log(condN/condNmin)).^2;
		Nopt     = Nv(find(cost==min(cost)))

		if(par.verbose)
			fprintf('\n(*) Number of optimum nodes = %d\n',Nopt);
		end

	end

	%
	% Send the best data-set stats
	%

	[z hz]            = GridDensity(log(s),wt,Nopt);             % Select "tau" Points
	[g tau error cKp] = MaxwellModes(z, t, Gt, par.prune);   % Get g_i

	%
	% Some Plotting
	%

	if(par.plotting)

		subplot(2,1,1)

		loglog(tau,g,'o-')
		xlabel('tau')
		ylabel('g')

		subplot(2,1,2)

		PlotMaxwellModes(g, tau, t, Gt)

		hold on 
		K = kernel(H,t,s);
		loglog(t, K, 'r-','LineWidth',2);
		xlabel('t')
		ylabel('G(t)')
		hold off

	end
  

  
	%
	% Some Printing
	%

	if(par.verbose)

		f1 = fopen('output/dmodes.dat','w');
		fprintf('(*) Condition number of matrix equation: %e\n',cKp);

		fprintf('\n\t\tModes\n\t\t-----\n\n');
		fprintf('i \t    g(i) \t    tau(i)\n');
		fprintf('---------------------------------------\n');
		for i = 1:length(g)
			fprintf('%d \t %9.5e \t %9.5e\n',i,g(i),tau(i));
			fprintf(f1,'%d \t %9.5e \t %9.5e\n',i,g(i),tau(i));
		end
		fprintf('\n');
		fclose(f1);

		if(par.Nopt == 0)
			f2 = fopen('output/Nopt.dat','w');
			for i = 1:npts
				fprintf(f2,'%e\t%e\t%e\n',Nv(i),ev(i),condN(i));
			end
			fclose(f2);
		end

		f3 = fopen('output/Gfitd.dat','w');
		K   = kernel(H,t,s);
		for i = 1:n
			fprintf(f3,'%e\t%e\t%e\n',t(i),K(i));
		end
		fclose(f3);
	end

end

% 
% Read Data Files
%
function   [t, Gt, s, H] = ReadData(fNameGst,fNameH)

	% Read input data
	[t Gt] = GetExpData(fNameGst);

	% Read continuous spectrum
	data    = load(fNameH);
	s       = data(:,1);
	H       = data(:,2);

end


