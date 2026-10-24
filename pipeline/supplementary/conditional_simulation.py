"""Conditional-design Gaussian simulation, distinct from the count pipeline."""
import numpy as np
import pandas as pd
import math
from scipy import linalg, stats
from .common import need
def write(frame, out, name):
    frame.to_csv(out / name, index=False)
def simulate_designs(designs, out, reps, seed):
    """Gaussian design benchmark, not a calibrated RNA-seq or limma simulator."""
    rows = []
    rng = np.random.default_rng(seed)
    for model, (x, xf, cages) in designs.items():
        for a in [x, xf]:
            need(np.linalg.matrix_rank(a) == a.shape[1] and len(a)>a.shape[1]+2, 'Simulation design not estimable')
        codes, levels = pd.factorize(cages)
        zc = np.eye(len(levels))[codes]
        cage_g = np.array([x[codes==i,1].mean() for i in range(len(levels))])
        fits = []
        for label, a in [('OLS_no_cage',x),('OLS_cage_fixed',xf)]:
            pinv = linalg.pinv(a)
            fits.append((label,a,pinv,float((pinv@pinv.T)[1,1]),len(a)-a.shape[1]))
        for rho in [0., .05, .20]:
            for corr in [0., .5]:
                for truth in [0., .3]:
                    # Total random variance is one. Correlated deterministic cage
                    # means are a separate, deliberately misspecified scenario.
                    errors = math.sqrt(1-rho)*rng.normal(size=(len(x),reps))
                    errors += math.sqrt(rho)*zc@rng.normal(size=(len(levels),reps))
                    y = truth*x[:,1,None]+errors+corr*zc@(cage_g-cage_g.mean())[:,None]
                    for label,a,pinv,v,df in fits:
                        co = pinv@y
                        residual = y-a@co
                        se = np.sqrt((residual**2).sum(axis=0)/df*v)
                        est = co[1]
                        good = np.isfinite(est)&np.isfinite(se)&(se>0)
                        crit = stats.t.ppf(.975,df)
                        covered = np.abs(est[good]-truth)<=crit*se[good]
                        sig = 2*stats.t.sf(np.abs(est[good]/se[good]),df)<.05
                        rate = float(covered.mean()) if good.any() else np.nan
                        rows.append(dict(model=model, estimator=label, random_cage_ICC=rho,
                            correlated_cage_shift=corr, true_beta=truth, replicates=reps,
                            valid=int(good.sum()), failure_rate=1-good.mean(),
                            bias=np.mean(est[good]-truth), empirical_SD=np.std(est[good],ddof=1),
                            RMS_estimated_SE=np.sqrt(np.mean(se[good]**2)), coverage=rate,
                            coverage_MCSE=math.sqrt(rate*(1-rate)/good.sum()),
                            mean_CI_width=np.mean(2*crit*se[good]), rejection_rate=sig.mean(),
                            rejection_MCSE=math.sqrt(sig.mean()*(1-sig.mean())/good.sum()),
                            error_scope='pointwise_not_FDR', simulation='Gaussian_design_only_not_limma'))
    write(pd.DataFrame(rows), out, '322_design_simulation.csv')
