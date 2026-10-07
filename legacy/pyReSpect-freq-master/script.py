
run_freq_py_v3 = '''#!/usr/bin/env python3
# run_all_test_freq_v3.py

import os, sys, shutil
import numpy as np
import pandas as pd

base_dir  = os.path.dirname(os.path.abspath(__file__))
tests_dir = os.path.join(base_dir, "tests")
out_root  = os.path.join(base_dir, "output_PY_freq")
os.makedirs(out_root, exist_ok=True)

sys.path.insert(0, base_dir)

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
    print(f"\\n{'='*40}\\n{fname}\\n{'='*40}")

    fpath  = os.path.join(tests_dir, fname)
    outdir = os.path.join(out_root, fname.replace(".dat",""))
    os.makedirs(outdir, exist_ok=True)

    # copia il .dat in outdir
    shutil.copy(fpath, os.path.join(outdir, fname))

    # vai in outdir PRIMA di importare/chiamare qualsiasi cosa
    orig_dir = os.getcwd()
    os.chdir(outdir)
    if not os.path.exists("output"):
        os.makedirs("output")

    lamC = None; Nopt = None; error = None; log10cond = None

    try:
        # import dentro il loop cosi' usa la working dir corrente
        from contSpec import getContSpec
        from discSpec import getDiscSpecMagic
        from common   import GetExpData

        par = base_par.copy()
        par["GexpFile"] = fname

        # contSpec scrive output/H.dat qui dentro outdir
        s, H, G0, lamC = getContSpec(par)

        # discSpec legge output/H.dat da outdir — ora lo trova!
        Nopt, g, tau, error = getDiscSpecMagic(par)

        # condition number
        w, Gexp, wexp = GetExpData(fname)
        if len(tau) > 0:
            S, W  = np.meshgrid(tau, w)
            ws    = S * W; ws2 = ws**2
            K     = np.vstack([ws2/(1+ws2), ws/(1+ws2)])
            Kp    = (wexp/Gexp).reshape(-1,1) * K
            log10cond = float(np.log10(np.linalg.cond(Kp)))

        print(f"[PY-FREQ] lamC={lamC:.3e}  Nopt={Nopt}  err={error:.4f}  log10cond={log10cond:.2f}")

    except Exception as e:
        import traceback
        print(f"[PY-FREQ] ERRORE su {fname}: {e}")
        traceback.print_exc()

    os.chdir(orig_dir)

    results.append({
        "test": fname, "lamC": lamC, "Nopt": Nopt,
        "error": error, "log10cond": log10cond
    })

df = pd.DataFrame(results)
print("\\n===== RISULTATI PYTHON FREQUENCY =====")
print(df.to_string(index=False))
df.to_csv(os.path.join(out_root, "results_summary.csv"), index=False)
print(f"\\n(*) Salvato: {out_root}/results_summary.csv")
'''

with open("output/run_all_test_freq_v3.py", "w") as f:
    f.write(run_freq_py_v3)
print("Done")
