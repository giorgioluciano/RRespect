%
% Function: contSpec(par)
%
% Using a simplified L-curve method to compute the continuous relaxation 
% spectra H(s) given G(t) from an input file
%
% Uses parameters from SetParameters.m file:
% These parameters are prefixed by par.XXXXX, where XXXXX is the paramter
% Definitions of these parameters are provided in the SetParameters.m file
%
% Output:   Yields the fitted H(s), and optionally critical lambda "lamC"
%           If verbose is on: rho-eta.dat (used to determine lamC)
%                             H.dat (the spectrum)
%                             Gfit.dat (G(t) implied by the extracted H(s))
%

function  [H lamC] = contSpec(par)

    %
    % Add appropriate subdirectories to search path
    %
	addpath('./continuous:./common')

	if nargin == 0
		par = SetParameters();  % Load in global settings
	end
  
	%
	% Load the data in
	% 
	if(par.verbose)
		fprintf('\n(*) Start\n(*) Loading Data File: %s...',par.GstFile);
	end

	[t Gexp] = GetExpData(par.GstFile);

	if(par.verbose)
		fprintf('done\n(*) Initial Set up...');
	end

	tic
  
	%
	% Set up some internal variables
	%

	n    = length(t);
	ns   = par.ns;    % discretization of 'tau'

	tmin = t(1);
	tmax = t(n);

	switch par.FreqEnd
		case 1
			smin = exp(-pi/2) * tmin; smax = exp(pi/2) * tmax;
		case 2
			smin = tmin; smax = tmax;
		case 3
			smin = exp(+pi/2) * tmin; smax = exp(-pi/2) * tmax;
	end 


	hs   = (smax/smin)^(1/(ns-1));
	s    = smin * hs.^[0:ns-1]';

	Hgs  = InitializeH(Gexp,t,s);

	te = toc;

	tic

	%
	% Find Optimum Lambda with 'lcurve'
	%

	if(par.verbose)
		fprintf('(%5.1f seconds)\n(*) Building the L-curve, ...',te);
	end

	if(par.lamC == 0)
		[lamC lam rho eta] = lcurve(Gexp, Hgs, t, s, par.SmFacLam);
	else
		lamC = par.lamC;
	end

	te = toc;

	if(par.verbose)
		fprintf('%e (%5.1f seconds)\n(*) Extracting the continuous spectrum, ...',lamC,te);
	end

	tic;

	%
	% Get the spectrum
	%
	
	H = LevenMarq(lamC, Gexp, Hgs, t, s);    
	te = toc;

	%
	% Print some datafiles
	%
	
	if(par.verbose)
		fprintf('done (%5.1f seconds)\n(*) Writing and Printing, ...',te);
	end

	if(par.verbose)

		if(par.lamC == 0)
			f1  = fopen('output/rho-eta.dat','w');
			for i = 1:length(lam)
				fprintf(f1,'%e\t%e\t%e\n',lam(i),rho(i),eta(i));
			end
			fclose(f1);
		end

		f2  = fopen('output/H.dat','w');
		for i = 1:ns
			fprintf(f2,'%e\t%e\n',s(i),H(i));
		end
		fclose(f2);

		f3  = fopen('output/Gfit.dat','w');
		K   = kernel(H,t,s);
		for i = 1:n
			fprintf(f3,'%e\t%e\n',t(i),K(i));
		end
		fclose(f3);

	end

	%
	% Graphing
	%    

	if(par.plotting)

		subplot(2,1,1)
		semilogx(s,H,'o-')
		xlabel('s')
		ylabel('H(s)')
		title('H')

		subplot(2,1,2);
		K = kernel(H,t,s);
		loglog(t,Gexp,'o',t, K, 'k-');
		xlabel('w')
		ylabel('Gt(exp), Gt(fit)')
		title('G(t)')

	end

	if(par.verbose)
		fprintf('done\n(*) End\n');
	end

end
