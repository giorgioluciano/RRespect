
run_all_python = '''#!/usr/bin/env python3
# run_all_test_time.py
# Batch runner per pyReSpect-time su tutti i test set
# Uso: metti questo script nella stessa cartella di contSpec.py e discSpec.py
#      e lancia: python run_all_test_time.py

import os, sys, shutil
import numpy as np
import pandas as pd

# aggiungi la cartella corrente al path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from contSpec import getContSpec
from discSpec import getDiscSpecMagic

# ── Cartella base (dove sta questo script) ────────────────────────────────────
base_dir = os.path.dirname(os.path.abspath(__file__))
tests_dir = os.path.join(base_dir, "tests")
out_root  = os.path.join(base_dir, "output_PY_time")
os.makedirs(out_root, exist_ok=True)

# ── Lista test ────────────────────────────────────────────────────────────────
test_files = sorted([f for f in os.listdir(tests_dir) if f.endswith(".dat")])

if not test_files:
    print("Nessun file .dat trovato in", tests_dir)
    sys.exit(1)

print(f"Trovati {len(test_files)} file di test: {test_files}")

# ── Parametri base (equivalenti al SetParameters R) ──────────────────────────
base_par = {
    "ns"          : 100,
    "lamC"        : 0,
    "SmFacLam"    : 0,
    "FreqEnd"     : 1,
    "verbose"     : True,
    "plotting"    : False,   # metti True se vuoi i PDF
    "lam_min"     : 1e-10,
    "lam_max"     : 1e3,
    "lamDensity"  : 3,
    "deltaBaseWeightDist": 0.2,
    "minTauSpacing"      : 1.25,
    "MaxNumModes"        : 0,
    "plateau"     : False,
}

results = []

for fname in test_files:
    fpath = os.path.join(tests_dir, fname)
    print(f"\\n{'='*40}")
    print(f"Running Python time-domain on: {fname}")
    print(f"{'='*40}")

    # cartella output per questo test
    outdir = os.path.join(out_root, fname.replace(".dat", ""))
    os.makedirs(outdir, exist_ok=True)

    # copia il .dat in outdir (come fa lo script R)
    shutil.copy(fpath, os.path.join(outdir, fname))

    # costruisci par specifico per questo test
    par = base_par.copy()
    par["GexpFile"] = fname   # solo il nome, lo trova in outdir

    # vai in outdir per leggere/scrivere
    orig_dir = os.getcwd()
    os.chdir(outdir)

    if not os.path.exists("output"):
        os.makedirs("output")

    try:
        H, lamC = getContSpec(par)

        Nopt, g, tau, error = getDiscSpecMagic(par)

        # leggi condKp da dmodes (non restituito direttamente)
        # calcoliamo il condition number dal kernel
        from common import GetExpData, getKernMat
        import numpy as np
        t, Gexp, wexp = GetExpData(fname)
        s = np.loadtxt("output/H.dat")[:, 0]
        kernMat = getKernMat(s, t)
        from scipy.optimize import nnls
        S, T = np.meshgrid(tau, t)
        K = np.exp(-T/S)
        Kp = (wexp/Gexp).reshape(-1,1) * K
        condKp = np.linalg.cond(Kp)

        print(f"[PY] {fname:15s}: lamC={lamC:.3e}  Nopt={Nopt}  err={error:.4f}  cond={condKp:.3e}")

        results.append({
            "test"      : fname,
            "lamC"      : lamC,
            "Nopt"      : Nopt,
            "error"     : error,
            "log10cond" : np.log10(condKp)
        })

        # copia H.dat e dmodes.dat nella cartella outdir root (gia' ci siamo)
        # sono gia' in output/ quindi ok

    except Exception as e:
        print(f"[PY] ERRORE su {fname}: {e}")
        results.append({
            "test": fname, "lamC": None,
            "Nopt": None, "error": None, "log10cond": None
        })

    os.chdir(orig_dir)

# ── Salva risultati ────────────────────────────────────────────────────────────
df = pd.DataFrame(results)
print("\\n===== RISULTATI BATCH PYTHON =====")
print(df.to_string(index=False))

out_csv = os.path.join(out_root, "results_summary.csv")
df.to_csv(out_csv, index=False)
print(f"\\n(*) Salvato: {out_csv}")
'''

with open("output/run_all_test_time.py", "w") as f:
    f.write(run_all_python)
print("run_all_test_time.py written")
