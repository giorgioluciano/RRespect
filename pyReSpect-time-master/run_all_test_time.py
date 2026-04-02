#!/usr/bin/env python3
# run_all_test_time.py  v2
# Batch runner pyReSpect-time — legge risultati direttamente dai file output

import os, sys, shutil
import numpy as np
import pandas as pd

base_dir  = os.path.dirname(os.path.abspath(__file__))
tests_dir = os.path.join(base_dir, "tests")
out_root  = os.path.join(base_dir, "output_PY_time")
os.makedirs(out_root, exist_ok=True)

sys.path.insert(0, base_dir)
from contSpec import getContSpec
from discSpec import getDiscSpecMagic

test_files = sorted([f for f in os.listdir(tests_dir) if f.endswith(".dat")])
print(f"Trovati {len(test_files)} file: {test_files}")

base_par = {
    "ns": 100, "lamC": 0, "SmFacLam": 0, "FreqEnd": 1,
    "verbose": True, "plotting": False,
    "lam_min": 1e-10, "lam_max": 1e3, "lamDensity": 3,
    "deltaBaseWeightDist": 0.2, "minTauSpacing": 1.25,
    "MaxNumModes": 0, "plateau": False,
}

results = []

for fname in test_files:
    fpath  = os.path.join(tests_dir, fname)
    outdir = os.path.join(out_root, fname.replace(".dat",""))
    os.makedirs(outdir, exist_ok=True)
    shutil.copy(fpath, os.path.join(outdir, fname))

    par = base_par.copy()
    par["GexpFile"] = fname

    orig_dir = os.getcwd()
    os.chdir(outdir)
    if not os.path.exists("output"):
        os.makedirs("output")

    lamC = None; Nopt = None; error = None; log10cond = None

    try:
        print(f"\n{'='*40}\n{fname}\n{'='*40}")
        H, lamC = getContSpec(par)
        Nopt, g, tau, error = getDiscSpecMagic(par)

        # condition number dal kernel discreto
        from common import GetExpData, getKernMat
        t, Gt, wt = GetExpData(fname)
        if len(tau) > 0:
            S, T      = np.meshgrid(tau, t)
            K         = np.exp(-T/S)
            Kp        = (wt/Gt).reshape(-1,1) * K
            log10cond = float(np.log10(np.linalg.cond(Kp)))

        print(f"[PY] lamC={lamC:.3e}  Nopt={Nopt}  err={error:.4f}  log10cond={log10cond:.2f}")

    except Exception as e:
        print(f"[PY] ERRORE: {e}")

    os.chdir(orig_dir)

    results.append({
        "test": fname, "lamC": lamC, "Nopt": Nopt,
        "error": error, "log10cond": log10cond
    })

df = pd.DataFrame(results)
print("\n===== RISULTATI BATCH PYTHON =====")
print(df.to_string(index=False))
df.to_csv(os.path.join(out_root, "results_summary.csv"), index=False)
print(f"\n(*) Salvato: {out_root}/results_summary.csv")
