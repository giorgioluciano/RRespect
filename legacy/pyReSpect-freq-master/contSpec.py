#
# Help to find continuous spectrum: 
#
# March 2019 major update:
# (i)   added plateau modulus G0 (also in pyReSpect-time) calculation
# (ii)  following Hansen Bayesian interpretation of Tikhonov to extract p(lambda): 
# (iii) simplifying lcurve (starting from high lambda to low) - super cost savings!
# (iv)  still showing lamdaC from prior method (although abridged lambda scale) comparison
#

from common import *

# HELPER FUNCTIONS

def InitializeH(Gexp, wexp, s, kernMat,  *argv):
    """
     Function: InitializeH(input)

     Input:  Gexp    = 2n*1 vector [G\'\';G"],
              wexp    = 2n*1 vector of weights
             s       = relaxation modes,
             kernMat = matrix for faster kernel evaluation
             G0      = optional; if plateau is nonzero	

     Output: H = guessed H
              G0 = optional guess if *argv is nonempty	
    """

    H    = -5.0 * np.ones(len(s)) + np.sin(np.pi * s)
    lam  = 1e0

    if len(argv) > 0:
        G0       = argv[0]
        Hlam, G0 = getH(lam, Gexp, wexp, H, kernMat, G0)		
        return Hlam, G0
    else:
        Hlam     = getH(lam, Gexp, wexp, H, kernMat)
        return Hlam


def getAmatrix(ns):
    """Generate symmetric matrix A = L\'\' * L required for error analysis:
       helper function for lcurve in error determination"""
    nl = ns - 2
    L  = np.diag(np.ones(ns-1), 1) + np.diag(np.ones(ns-1),-1) + np.diag(-2. * np.ones(ns))
    L  = L[1:nl+1,:]
    return np.dot(L.T, L)

def getBmatrix(H, kernMat, Gexp, wexp, *argv):
    """get the Bmatrix required for error analysis; helper for lcurve()"""
    n   = int(len(Gexp)/2);
    ns  = len(H);
    nl  = ns - 2;
    r   = np.zeros(n);

    Kmatrix = np.dot((wexp/Gexp).reshape(2*n,1), np.ones((1,ns)));
    Jr      = -kernelD(H, kernMat) * Kmatrix;    

    if len(argv) > 0:
        G0 = argv[0]
        r  = wexp * (1. - kernel_prestore(H, kernMat, G0)/Gexp)
    else:
        r = wexp * (1. - kernel_prestore(H, kernMat)/Gexp)

    B = np.dot(Jr.T, Jr) + np.diag(np.dot(r.T, Jr))
    return B

def lcurve(Gexp, wexp, Hgs, kernMat, par, *argv):
    """
     Function: lcurve(input)
     March 2019: starting from large lambda to small cuts calculation time by a lot
    """

    if par['plateau']:
        G0 = argv[0]

    npoints = int(par['lamDensity'] * (np.log10(par['lam_max']) - np.log10(par['lam_min'])))
    hlam    = (par['lam_max']/par['lam_min'])**(1./(npoints-1.))	
    lam     = par['lam_min'] * hlam**np.arange(npoints)
    eta     = np.zeros(npoints)
    rho     = np.zeros(npoints)
    logP    = np.zeros(npoints)
    H       = Hgs.copy()
    n       = len(Gexp)
    ns      = len(H)
    nl      = ns - 2
    logPmax = -np.inf
    Hlambda = np.zeros((ns, npoints))

    Amat       = getAmatrix(len(H))
    _, LogDetN = np.linalg.slogdet(Amat)

    for i in reversed(range(len(lam))):

        lamb    = lam[i]

        if par['plateau']:
            H, G0   = getH(lamb, Gexp, wexp, H, kernMat, G0)			
            rho[i]  = np.linalg.norm(wexp * (1. - kernel_prestore(H, kernMat, G0)/Gexp))
            Bmat    = getBmatrix(H, kernMat, Gexp, wexp, G0)			
        else:
            H       = getH(lamb, Gexp, wexp, H, kernMat)
            rho[i]  = np.linalg.norm(wexp * (1. - kernel_prestore(H,kernMat)/Gexp))
            Bmat    = getBmatrix(H, kernMat, Gexp, wexp)

        eta[i]       = np.linalg.norm(np.diff(H, n=2))
        Hlambda[:,i] = H

        _, LogDetC = np.linalg.slogdet(lamb*Amat + Bmat)
        V          =  rho[i]**2 + lamb * eta[i]**2		
        logP[i]    = -V + 0.5 * (LogDetN + ns*np.log(lamb) - LogDetC) - lamb

        if(logP[i] > logPmax):
            logPmax = logP[i]
        elif(logP[i] < logPmax - 18):
            break		

    lam  = lam[i:]
    logP = logP[i:]
    eta  = eta[i:]
    rho  = rho[i:]
    logP = logP - max(logP)
    Hlambda = Hlambda[:,i:]

    plam = np.exp(logP); plam = plam/np.sum(plam)
    lamM = np.exp(np.sum(plam*np.log(lam)))

    if par['SmFacLam'] > 0:
        lamM = np.exp(np.log(lamM) + par['SmFacLam']*(max(np.log(lam)) - np.log(lamM)));
    elif par['SmFacLam'] < 0:
        lamM = np.exp(np.log(lamM) + par['SmFacLam']*(np.log(lamM) - min(np.log(lam))));

    if par['plotting']:
        plt.clf()
        plt.axvline(x=lamM, c='gray', label=r'$\lambda_m$ = {0:0.2e}'.format(lamM))
        plt.ylim(-20,1)
        plt.plot(lam, logP, 'o-')
        plt.xscale('log')
        plt.xlabel(r'$\lambda$')
        plt.ylabel(r'$\log\,p(\lambda)$')
        plt.legend(loc='upper left', fontsize=14)
        plt.tight_layout()
        plt.savefig('output/logP.pdf')

    return lamM, lam, rho, eta, logP, Hlambda

def getH(lam, Gexp, wexp, H, kernMat, *argv):
    """minimize_H  V(lambda)"""
    if len(argv) > 0:
        Hplus= np.append(H, argv[0])
        res_lsq = least_squares(residualLM, Hplus, jac=jacobianLM, args=(lam, Gexp, wexp, kernMat))
        return res_lsq.x[:-1], res_lsq.x[-1]
    else:
        res_lsq = least_squares(residualLM, H, jac=jacobianLM, args=(lam, Gexp, wexp, kernMat))			
        return res_lsq.x

def residualLM(H, lam, Gexp, wexp, kernMat):
    """HELPER FUNCTION: Gets Residuals r"""

    n   = int(kernMat.shape[0]/2);
    ns  = kernMat.shape[1];
    nl  = ns - 2;
    r   = np.zeros(2*n + nl);

    if len(H) > ns:
        G0       = H[-1]
        H        = H[:-1]
        r[0:2*n] = wexp * (1. - kernel_prestore(H, kernMat, G0)/Gexp)
    else:
        r[0:2*n] = wexp * (1. - kernel_prestore(H,kernMat)/Gexp)

    r[2*n:2*n+nl] = np.sqrt(lam) * np.diff(H, n=2)

    return r

def jacobianLM(H, lam, Gexp, wexp, kernMat):
    """HELPER FUNCTION: Gets Jacobian J"""

    n   = int(kernMat.shape[0]/2);
    ns  = kernMat.shape[1];
    nl  = ns - 2;

    L  = np.diag(np.ones(ns-1), 1) + np.diag(np.ones(ns-1),-1) + np.diag(-2. * np.ones(ns))	     
    L  = L[1:nl+1,:]

    Kmatrix = np.dot((wexp/Gexp).reshape(2*n,1), np.ones((1,ns)))

    if len(H) > ns:
        G0     = H[-1]
        H      = H[:-1]
        Jr  = np.zeros((2*n + nl,ns+1))	
        Jr[0:2*n, 0:ns]     = -kernelD(H, kernMat) * Kmatrix
        Jr[0:n, ns]         = -wexp[:n]/Gexp[:n]
        Jr[2*n:2*n+nl,0:ns] = np.sqrt(lam) * L
        Jr[2*n:2*n+nl, ns]  = np.zeros(nl)
    else:
        Jr  = np.zeros((2*n + nl,ns))	
        Jr[0:2*n, 0:ns]     = -kernelD(H, kernMat) * Kmatrix
        Jr[2*n:2*n+nl,0:ns] = np.sqrt(lam) * L

    return Jr

def kernelD(H, kernMat):
    """outputs the 2n*ns dimensional vector DK(H)(w)"""
    n   = int(kernMat.shape[0]/2);
    ns  = kernMat.shape[1];
    Hsuper  = np.dot(np.ones((2*n,1)), np.exp(H).reshape(1, ns))       
    DK      = kernMat * Hsuper
    return DK

def getContSpec(par):

    if par['verbose']:
        print('\n(*) Start\n(*) Loading Data File: {}...'.format(par['GexpFile']))

    w, Gexp, wexp = GetExpData(par['GexpFile'])

    if par['verbose']:
        print('(*) Initial Set up...', end="")

    n    = len(w)
    ns   = par['ns']

    wmin = w[0];
    wmax = w[n-1];

    if par['FreqEnd'] == 1:
        smin = np.exp(-np.pi/2)/wmax; smax = np.exp(np.pi/2)/wmin		
    elif par['FreqEnd'] == 2:
        smin = 1./wmax; smax = 1./wmin				
    elif par['FreqEnd'] == 3:
        smin = np.exp(+np.pi/2)/wmax; smax = np.exp(-np.pi/2)/wmin

    hs   = (smax/smin)**(1./(ns-1))
    s    = smin * hs**np.arange(ns)

    kernMat = getKernMat(s, w)
    tic  = time.time()

    if par['plateau']:
        Hgs, G0  = InitializeH(Gexp, wexp, s, kernMat, np.min(Gexp))		
    else:
        Hgs      = InitializeH(Gexp, wexp, s, kernMat)	

    if par['verbose']:
        te   = time.time() - tic
        print('\t({0:.1f} seconds)\n(*) Building the L-curve ...'.format(te), end="")	
        tic  = time.time()

    if par['lamC'] == 0:
        if par['plateau']:
            lamC, lam, rho, eta, logP, Hlam = lcurve(Gexp, wexp, Hgs, kernMat, par, G0)
        else:
            lamC, lam, rho, eta, logP, Hlam = lcurve(Gexp, wexp, Hgs, kernMat, par)
    else:
        lamC = par['lamC']

    if par['verbose']:
        te = time.time() - tic
        print('({0:.1f} seconds)\n(*) Extracting CRS, ...\n\t... lamM = {1:0.3e}; '.
              format(te, lamC), end="")
        tic = time.time()

    if par['plateau']:
        H, G0  = getH(lamC, Gexp, wexp, Hgs, kernMat, G0);
        print('G0 = {0:0.3e} ...'.format(G0), end="")
    else:
        H  = getH(lamC, Gexp, wexp, Hgs, kernMat);

    #----------------------
    # Print some datafiles
    #----------------------
    if par['verbose']:
        te = time.time() - tic
        print('done ({0:.1f} seconds)\n(*) Writing and Printing, ...'.format(te), end="")

        # FIX: salva H.dat sempre (con 3 colonne: s, H, dH=0)
        # dH e' disponibile solo con plotting=True; usiamo zeros come placeholder
        if par['plateau']:
            K = kernel_prestore(H, kernMat, G0)
            np.savetxt('output/H.dat', np.c_[s, H, np.zeros(len(s))], fmt='%e',
                       header='G0 = {0:0.3e}'.format(G0))
        else:
            K = kernel_prestore(H, kernMat)
            np.savetxt('output/H.dat', np.c_[s, H, np.zeros(len(s))], fmt='%e')

        np.savetxt('output/Gfit.dat', np.c_[w, K[:n], K[n:]], fmt='%e')

        # print Hlam, rho-eta, and logP if lcurve has been visited
        if par['lamC'] == 0:

            if os.path.exists("output/Hlam.dat"):
                os.remove("output/Hlam.dat")

            fHlam = open('output/Hlam.dat','ab')
            for i, lamb in enumerate(lam):
                np.savetxt(fHlam, Hlam[:,i])	
            fHlam.close()	

            np.savetxt('output/logPlam.dat', np.c_[lam, logP])
            np.savetxt('output/rho-eta.dat', np.c_[lam, rho, eta], fmt='%e')

    #------------
    # Graphing
    #------------

    if par['plotting']:

        plt.clf()
        plt.semilogx(s, H, 'o-')
        plt.xlabel(r'$s$')
        plt.ylabel(r'$H(s)$')

        if par['lamC'] == 0:
            plam = np.exp(logP); plam = plam/np.sum(plam)			
            Hm   = np.zeros(len(s))
            Hm2  = np.zeros(len(s))
            cnt  = 0
            for i in range(len(lam)):	
                if plam[i] > 0.1:
                    Hm   += Hlam[:,i]
                    Hm2  += Hlam[:,i]**2
                    cnt  += 1

            Hm = Hm/cnt
            dH = np.sqrt(Hm2/cnt - Hm**2)

            plt.semilogx(s, H+2.5*dH, c='gray', alpha=0.5)
            plt.semilogx(s, H-2.5*dH, c='gray', alpha=0.5)

            if par['verbose']:			
                if par['plateau']:
                    K = kernel_prestore(Hm, kernMat, G0)	
                    np.savetxt('output/H.dat', np.c_[s, H, dH], fmt='%e',
                               header='G0 = {0:0.3e}'.format(G0))
                else:
                    K = kernel_prestore(Hm, kernMat)	
                    np.savetxt('output/H.dat', np.c_[s, H, dH], fmt='%e')

                np.savetxt('output/Gfit.dat', np.c_[w, K[:n], K[n:]], fmt='%e')

        plt.tight_layout()
        plt.savefig('output/H.pdf')

        plt.clf()

        if par['plateau']:
            K = kernel_prestore(H, kernMat, G0)	
        else:
            K = kernel_prestore(H, kernMat)

        plt.loglog(w, Gexp[:n], 'o')
        plt.loglog(w, K[:n], 'k-')
        plt.loglog(w, Gexp[n:], 'o')
        plt.loglog(w, K[n:], 'k-')
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'$G^{*}$')
        plt.savefig('output/Gfit.pdf')

        try:
            lam
        except NameError:
            print("lamC prespecified, so not printing rho-eta.pdf/dat")
        else:
            plt.clf()
            plt.plot(rho, eta)
            rhost = np.interp(lamC, lam, rho)
            etast = np.interp(lamC, lam, eta)
            plt.plot(rhost, etast, 'o', c='k')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$\rho$')
            plt.ylabel(r'$\eta$')
            plt.savefig('output/rho-eta.pdf')

    if par['verbose']:
        print('done\n(*) End\n')

    if not par['plateau']:
        G0 = 0.

    return s, H, G0, lamC

def guiFurnishGlobals(par):

    from matplotlib import rcParams

    w, Gexp, wexp = GetExpData(par['GexpFile'])

    if par['verbose']:
        print('(*) Initial Set up...', end="")

    n    = len(w)
    ns   = par['ns']

    wmin = w[0];
    wmax = w[n-1];

    if par['FreqEnd'] == 1:
        smin = np.exp(-np.pi/2)/wmax; smax = np.exp(np.pi/2)/wmin		
    elif par['FreqEnd'] == 2:
        smin = 1./wmax; smax = 1./wmin				
    elif par['FreqEnd'] == 3:
        smin = np.exp(+np.pi/2)/wmax; smax = np.exp(-np.pi/2)/wmin

    hs   = (smax/smin)**(1./(ns-1))
    s    = smin * hs**np.arange(ns)

    kernMat = getKernMat(s, w)

    par['verbose']  = False
    par['plotting'] = False

    lam, rho, eta = np.loadtxt('output/rho-eta.dat', unpack=True)

    rcParams['axes.labelsize']  = 14 
    rcParams['xtick.labelsize'] = 12
    rcParams['ytick.labelsize'] = 12 
    rcParams['legend.fontsize'] = 12
    rcParams['lines.linewidth'] = 2

    if par['plotting']: plt.clf()

    return s, w, kernMat, Gexp, wexp, par, lam, rho, eta


if __name__ == '__main__':
    par = readInput('inp.dat')
    s, H, G0, lamC = getContSpec(par)
