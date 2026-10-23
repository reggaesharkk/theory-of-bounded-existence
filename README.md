# The Theory of Bounded Existence — Appendix B v1.0 (Lean 4 Boundary-Verified, Multi-Solver Poincaré & VRK v0.3.1 Edition)

[![Verify Appendix B v1.0 (Lean 4 Lake & Python Suite)](https://github.com/reggaesharkk/theory-of-bounded-existence/actions/workflows/verify_ci.yml/badge.svg)](https://github.com/reggaesharkk/theory-of-bounded-existence/actions/workflows/verify_ci.yml)

* **Author:** Prince Upadhyay (`starxshark@gmail.com`, ORCID: [0009-0002-1170-4487](https://orcid.org/0009-0002-1170-4487))
* **Release Date:** October 23, 2026
* **GitHub Repository & Passing CI:** [reggaesharkk/theory-of-bounded-existence](https://github.com/reggaesharkk/theory-of-bounded-existence) | [GitHub Actions CI (`verify_ci.yml`)](https://github.com/reggaesharkk/theory-of-bounded-existence/actions/workflows/verify_ci.yml)
* **Upstream Lineage DOIs:**
  * IntentSeal VRK v0.3.1: [10.5281/zenodo.23018459](https://doi.org/10.5281/zenodo.23018459) (Lineage: [10.5281/zenodo.23016445](https://doi.org/10.5281/zenodo.23016445); Ecosystem: [10.5281/zenodo.23474505](https://doi.org/10.5281/zenodo.23474505))
  * Physics v1.0 / Appendix A v1.5 Soft-Wall Completion: [10.5281/zenodo.23534449](https://doi.org/10.5281/zenodo.23534449) | [10.5281/zenodo.23501526](https://doi.org/10.5281/zenodo.23501526) | [10.17605/OSF.IO/N3DQU](https://doi.org/10.17605/OSF.IO/N3DQU)
  * Operationalizing Self-Modeling Organization ($\Psi$): [10.17605/OSF.IO/8CFQW](https://doi.org/10.17605/OSF.IO/8CFQW) ([Internet Archive](https://archive.org/details/psi-self-modeling-organization) | [GitHub](https://github.com/reggaesharkk/psi-self-modeling-benchmark))
  * Low-Rank Spectral Compression (LRSC v1.2, $r_{\text{LRSC}} = 11/99 = 1/9$): [10.17605/OSF.IO/NM5BW](https://doi.org/10.17605/OSF.IO/NM5BW)
* **Tolerance-Free VRK v0.3.1 Archival Consistency Digest (SHA-256):**
  `597b9ae3679dcd2e2bf17cf7c1616c463f48424ff1fd8453395936ecba066bbc`
  *(Seals the exact rational Hopf parameters, polynomial coefficients, boolean test gates, and SHA-256 hashes of `ExistenceVerification.lean` and `lakefile.lean`; mathematical correctness rests on the symbolic, numerical, and Lean proofs below.)*

## Overview, Scope & Six-Layer Verification Architecture

This archive contains **14 files at the root directory** (13 release files verified by `MANIFEST_SHA256.txt`, plus `MANIFEST_SHA256.txt` itself) for **Appendix B v1.0: The Theory of Bounded Existence**, an illustrative, dimensionless 4D dynamical system with unit-normalized structural constants ($a=b=c=f=g=\mathcal{E}_{\text{avail}}=1$) designed to study the mathematical mechanics of bounded self-referential feedback:

1. **Globally Smooth ($C^\infty$) Quadratic Self-Feedback Ratio (`R_self_sq`):**
   Every component (`ExistenceVerification.lean`, `verification_existence_v1_0.py`, `dashboard.html`, `launch_existence_sandbox.py`, and `theory_of_existence_v1_0.tex`) uses the smooth quadratic ratio
   $$R_{\text{self}}(s) = \frac{\text{fb}(s)^2}{\text{fb}(s)^2 + (f\mathcal{P}_0)^2 + \sigma_0^2}, \qquad \text{fb}(s) = e\tanh(\dot{I}/\epsilon_I)\,\dot{H}_{\text{obs}}, \quad \sigma_0 = 0.05.$$
   At any equilibrium ($\dot{H}_{\text{obs}}^* = 0$), $R_{\text{self}}^* = 0$, $\Psi^* = 0$, and $\nabla R_{\text{self}}|_ * = (0,0,0,0)$ identically (so $\Psi(t)$ measures non-equilibrium feedback activity, crossing $0.50$ when $\kappa_{\text{meta}}R_{\text{self}}(1-\Psi)$ exceeds decay $\mu_\Psi\Psi$). Consequently, row 4 of the $4\times 4$ equilibrium Jacobian $J_4$ is $[0, 0, 0, -\mu_\Psi]$ and $\lambda_{\text{lucid}}$ decouples from linear equilibrium stability (`ANALYTIC_THEOREM`). Only the product $\Lambda_{\text{eff}} = \lambda_{\text{sim}}r_{\text{LRSC}}$ enters the ODEs (fixing $r_{\text{LRSC}} = 1/9$ sets the scale where $\lambda_{\text{sim}} = 7.20$ corresponds to $\Lambda_{\text{eff}} = 0.80$).
2. **Closed-Form $e_c(\mathcal{P}_0)$ Domain ($\mathcal{P}_0 > 10/99$), Two-Gate Error Bound ($\mathcal{P}_0 \ge 0.105$), & Global Monotonicity ($de_c/d\mathcal{P}_0 > 0$):**
   In the $\epsilon_H, \epsilon_I \to 0^+$ limit on the unregularized interior positivity domain $\mathcal{P}_0 > \mathcal{P}_{\text{crit}} = 10/99 = 0.\overline{10} \approx 0.10101$ (where $H_{\text{obs},\text{v1.4}}^* = 100 - (10+\mathcal{P}_0)/\mathcal{P}_0 > 0$), the open-loop Hopf threshold reduces to the rational function $e_c(\mathcal{P}_0) = \mathcal{P}_0(\mathcal{P}_0+20)(11\mathcal{P}_0+10)(\mathcal{P}_0^2+120\mathcal{P}_0+100)/[100(\mathcal{P}_0+10)(\mathcal{P}_0^2+130\mathcal{P}_0+200)]$ with strictly positive derivative $de_c/d\mathcal{P}_0 > 0$ (all 8 numerator and 7 denominator polynomial coefficients strictly positive). At $\epsilon_H = \epsilon_I = 0.05$, $e_c^{\text{CL}}(\mathcal{P}_0, 0) = [e_c(\mathcal{P}_0) + O(\epsilon_H^{-1}\operatorname{sech}^2(H_{\text{obs},\text{sw}}^*/\epsilon_H))]/\tanh(\mathcal{P}_0/\epsilon_I)$: for $\mathcal{P}_0 \ge 0.105$, $H_{\text{obs},\text{sw}}^* \ge 3.7619\text{ bits} > 75\epsilon_H$ suppresses the entropy-wall term $\operatorname{sech}^2(H_{\text{obs},\text{sw}}^*/\epsilon_H) \le 1.78 \times 10^{-65}$ to machine precision, leaving only the explicit factor $1/\tanh(\mathcal{P}_0/\epsilon_I)$ ($|e_c^{\text{CL}} - e_c| \le 3.74 \times 10^{-4}$ at $\mathcal{P}_0 = 0.105$ and $< 10^{-16}$ at $\mathcal{P}_0 \ge 1.0$). For $0 < \mathcal{P}_0 \le 10/99$ (including the Minimum Slider corner $\mathcal{P}_0 = 0.10$), the soft-wall threshold $e_c^{\text{CL}}(0.10, 0) = 0.099132$ governs at the boundary-layer root $H_{\text{obs},\text{sw}}^* = 0.1295\text{ bits}$.
3. **Zero-Drive ($\mathcal{P}_0 = 0$) Origin Stability, Transcritical Threshold $\Lambda_{\text{eff}} = 1$, and Bistability:**
   At $\mathcal{P}_0 = 0$, the origin $x_{\text{origin}}^* = (0,0,0,0)$ is an exact equilibrium with Jacobian eigenvalues $\lambda_{1,3} = [-(b+g) \pm \sqrt{(b-g)^2 + 4a\mathcal{E}_{\text{avail}}f\Lambda_{\text{eff}}}]/2$, $\lambda_2 = -\dot{I}_{\text{erase}}/\epsilon_H = -200$, and $\lambda_4 = -\mu_\Psi = -0.60$. For $0 \le \Lambda_{\text{eff}} < bg/(a\mathcal{E}_{\text{avail}}f) = 1$ (including Regime IVb where $\Lambda_{\text{eff}} = 0.80 < 1$, giving eigenvalues $\{-0.105573, -1.894427, -200.0, -0.60\}$), the origin is locally asymptotically stable while undergoing a transcritical bifurcation at $\Lambda_{\text{eff}} = 1$. Thus Regime IVb is bistable: near-origin initial conditions (e.g., $x(0) = (10^{-3})^4$) decay to $(0,0,0,0)$ ($\|x(320)\|_2 < 10^{-12}$), whereas trajectories in the macroscopic interior basin (verified across 3 interior initial conditions) converge to the periodic limit cycle ($T = 19.387103\text{ s}$, $\overline{\mathcal{P}}_{\text{eff},T} = 2.190321$, 3-IC spread $< 10^{-8}$, tightened Floquet $\rho_{\text{FD,tight}}(DP_3) < 10^{-7}$).
4. **Lean 4 Pointwise Boundary/Slope Formalization (`ExistenceVerification.lean` + `BoundaryVerification.lean`) & Classical ODE Step:**
   `lakefile.lean` registers both `BoundaryVerification` and `ExistenceVerification` as default library targets (`0 sorry`, `0 nlinarith`, `#print axioms` audit lines included). `ExistenceVerification.lean` machine-checks 15 pointwise vector-field theorems over Mathlib's $\mathbb{R}$: `R_self_sq_nonneg`, the pointwise linear comparison inequality `dgamma_dt_4D_ge_neg_decay` ($d\gamma/dt \ge -b\gamma$), the derived ceiling bounds `dHobs_dt_4D_le_M_H` and `dIdot_dt_4D_nonpos_at_ceiling` (at the conservative worst-case ceiling $\dot{I}_{\max,4\mathrm{D}} = (e M_H + f\mathcal{P}_{\text{eff},\max})/g$), the unified 8-face certificate `existence_box_8face_nagumo_certificate` ($\mathcal{P}_0 > 0$), and the zero-drive Face 5 theorem `dIdot_dt_4D_nonneg_at_boundary_zero_drive` ($\mathcal{P}_0 \ge 0$). Global forward-invariance of $\mathcal{B}_4$ for continuous-time $C^\infty$ ODE trajectories and strict interior positivity $\gamma(t) \ge \gamma(0)e^{-bt} > 0$ follow by classical application of the Nagumo–Bony–Brezis theorem and Grönwall's lemma. Layer 6 of `verification_existence_v1_0.py` performs an offline static source scan so Python users can run the script without a local Lean toolchain, while full Lean 4 compilation (`lake build`) is executed in GitHub Actions CI (`verify_ci.yml`).
5. **Multi-Tolerance Poincaré, Whole-Period Integral, Floquet Inequalities & Min/Max Extremes Audit:**
   Across $\text{rtol} \in \{10^{-8}, 10^{-10}, 10^{-12}\}$, `LSODA` and `DOP853` agree on the Poincaré section $\Sigma = \{\gamma = 2.0, \dot{\gamma} > 0\}$, with 6-digit converged whole-period integral $\overline{\mathcal{P}}_{\text{eff},T} = 3.266554$ in Regime IVa ($2.190321$ in Regime IVb). Because Liouville trace contraction $\exp(\int_0^T \operatorname{div} F\,dt) \sim 10^{-3010} \ll 1$ places the $(H_{\text{obs}}, \dot{I})$ return-map eigenvalues at the double-precision finite-difference noise floor ($\sim 10^{-11}$ to $10^{-9}$), finite-difference Floquet spectral radii and 3-IC spreads are certified strictly as upper-bound inequalities ($\Delta_{\text{IC}} < 10^{-8}$, $\rho_{\text{FD,tight}}(DP_{2,H\dot{I}}) < 10^{-7}$ across Regimes II, IIIb, IVa, and IVb), alongside the exact analytical integral $\mu_{\Psi,\text{exact}} = \exp(-\int_0^T (\kappa_{\text{meta}}R_{\text{self}} + \mu_\Psi)\,dt) = 2.537698 \times 10^{-6}$ in Regime II, and verified Minimum/Maximum slider bounds (`diff_min < 1e-7` at `rtol=1e-10` and `< 1e-8` at `rtol=1e-12`).

## Quick Verification Commands

```bash
# 1. Verify SHA-256 checksums of all 13 payload files listed in MANIFEST_SHA256.txt:
sha256sum -c MANIFEST_SHA256.txt

# 2. Build both Lean 4 libraries (BoundaryVerification and ExistenceVerification):
lake exe cache get
lake build

# 3. Run the 6-layer symbolic, Poincaré, Floquet, Min/Max & VRK v0.3.1 verification suite:
pip install -r requirements.txt
python3 verification_existence_v1_0.py

# 4. Launch the interactive 4D ODE simulator locally:
python3 launch_existence_sandbox.py
```
