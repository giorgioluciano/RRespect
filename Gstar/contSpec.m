%
% Function: contSpec(par)
%
% Using a simplified L-curve method to compute the continuous relaxation 
% spectra H(s) given G*(w) from an input file
%
% Uses parameters from SetParameters.m file:
% These parameters are prefixed by par.XXXXX, where XXXXX is the paramter
% Definitions of these parameters are provided in the SetParameters.m file
%
% Output:   Yields the fitted H(s), and optionally critical lambda "lamC"
%           If verbose is on: rho-eta.dat (used to determine lamC)
%                             H.dat (the spectrum)
%                             Gfit.dat (G* implied by the extracted H(s))
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

	[w Gexp] = GetExpData(par.GstFile);

	if(par.verbose)
		fprintf('done\n(*) Initial Set up...');
	end

	tic

	%
	% Set up some internal variables
	%
	n    = length(w);
	ns   = par.ns;    % discretization of 'tau'

	wmin = w(1);
	wmax = w(n);

	switch par.FreqEnd
	case 1
		smin = exp(-pi/2)/wmax; smax = exp(pi/2)/wmin;
	case 2
		smin = 1/wmax; smax = 1/wmin;
	case 3
		smin = exp(+pi/2)/wmax; smax = exp(-pi/2)/wmin;
	end 

	hs   = (smax/smin)^(1/(ns-1));
	s    = smin * hs.^[0:ns-1]';

	Hgs  = InitializeH(Gexp,w,s);

	t = toc;

	tic
	%
	% Find Optimum Lambda with 'lcurve'
	%
	if(par.lamC == 0)
		[lamC lam rho eta] = lcurve(Gexp, Hgs, w, s, par.SmFacLam);
	else
		lamC = par.lamC;
	end

	t = toc;

	if(par.verbose)
		fprintf('%e (%5.1f seconds)\n(*) Extracting the continuous spectrum, ...',lamC,t);
	end

	tic;
	
	% Get the spectrum

	H = LevenMarq(lamC, Gexp, Hgs, w, s);    
	t = toc;

	%
	% Print some datafiles
	%
	if(par.verbose)
		fprintf('done (%5.1f seconds)\n(*) Writing and Printing, ...',t);
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
		K   = kernel(H,w,s);
		for i = 1:n
			fprintf(f3,'%e\t%e\t%e\n',w(i),K(i),K(n+i));
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
		K = kernel(H,w,s);
		loglog(w,Gexp(1:n),'o',w,K(1:n),'k-',w,Gexp(n+1:2*n),'s',w,K(n+1:2*n),'k-');
		xlabel('w')
		ylabel('G*(exp), G*(fit)')
		title('Gp')

	end

	if(par.verbose)
		fprintf('done\n(*) End\n');
	end

end
