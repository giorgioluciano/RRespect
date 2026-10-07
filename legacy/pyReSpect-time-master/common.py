# common.py — fixed for scipy >= 1.12
# cumtrapz → cumulative_trapezoid

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import nnls, minimize, least_squares

# compatibilità scipy vecchio/nuovo
try:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except ImportError:
    from scipy.integrate import cumtrapz

import time
import os

plt.style.use('ggplot')

def readInput(fname='inp.dat'):
    par = {}
    for line in open(fname):
        li = line.strip()
        if len(li) > 0 and not li.startswith("#"):
            li = line.rstrip('\n').split(':')
            key = li[0].strip()
            tmp = li[1].strip()
            val = eval(tmp)
            par[key] = val
    if not os.path.exists("output"):
        os.makedirs("output")
    return par

def GetExpData(fname):
    try:
        data = np.loadtxt(fname)
        cols = data.shape[1]
        if cols == 2:
            to = data[:,0]; Gto = data[:,1]
        else:
            to = data[:,0]; Gto = data[:,1]; wGo = data[:,2]
    except OSError:
        print('*Error*: G(t) data file not found or incorrectly formatted')
        quit()

    to, indices = np.unique(to, return_index=True)
    Gto = Gto[indices]
    if cols == 3:
        wGo = wGo[indices]

    if cols == 2:
        f  = interp1d(to, Gto, fill_value="extrapolate")
        t  = np.geomspace(np.min(to), np.max(to), 100)
        Gt = f(t)
        return t, Gt, np.ones(len(Gt))
    else:
        return to, Gto, wGo

def getKernMat(s, t):
    ns  = len(s)
    hsv = np.zeros(ns)
    hsv[0]      = 0.5 * np.log(s[1]/s[0])
    hsv[ns-1]   = 0.5 * np.log(s[ns-1]/s[ns-2])
    hsv[1:ns-1] = 0.5 * (np.log(s[2:ns]) - np.log(s[0:ns-2]))
    S, T = np.meshgrid(s, t)
    return np.exp(-T/S) * hsv

def kernel_prestore(H, kernMat, *argv):
    G0 = argv[0] if len(argv) > 0 else 0.
    return np.dot(kernMat, np.exp(H)) + G0
