# Bounded Self-Referential Dynamics in a 4D Soft-Wall ODE System — Appendix B v1.0

[![Verify Appendix B v1.0 (Lean 4 Lake & Python Suite)](https://github.com/reggaesharkk/theory-of-bounded-existence/actions/workflows/verify_ci.yml/badge.svg)](https://github.com/reggaesharkk/theory-of-bounded-existence/actions/workflows/verify_ci.yml)

* **Author:** Prince Upadhyay (`starxshark@gmail.com`, ORCID: [0009-0002-1170-4487](https://orcid.org/0009-0002-1170-4487))
* **Release Date:** October 23, 2026
* **GitHub Repository & Passing CI:** [reggaesharkk/theory-of-bounded-existence](https://github.com/reggaesharkk/theory-of-bounded-existence) | [GitHub Actions CI (`verify_ci.yml`)](https://github.com/reggaesharkk/theory-of-bounded-existence/actions/workflows/verify_ci.yml)
* **Upstream Lineage DOIs:**
  * Physics v1.0 / Appendix A v1.5 Soft-Wall Completion: [10.5281/zenodo.23534449](https://doi.org/10.5281/zenodo.23534449) | [10.5281/zenodo.23501526](https://doi.org/10.5281/zenodo.23501526) | [10.17605/OSF.IO/N3DQU](https://doi.org/10.17605/OSF.IO/N3DQU)
  * IntentSeal VRK v0.3.1: [10.5281/zenodo.23018459](https://doi.org/10.5281/zenodo.23018459) (Lineage: [10.5281/zenodo.23016445](https://doi.org/10.5281/zenodo.23016445); Ecosystem: [10.5281/zenodo.23474505](https://doi.org/10.5281/zenodo.23474505))
  * Operationalizing Self-Modeling Organization ($\Psi$): [10.17605/OSF.IO/8CFQW](https://doi.org/10.17605/OSF.IO/8CFQW) ([Internet Archive](https://archive.org/details/psi-self-modeling-organization) | [GitHub](https://github.com/reggaesharkk/psi-self-modeling-benchmark))
  * Low-Rank Spectral Compression (LRSC v1.2, $r_{\text{LRSC}} = 11/99 = 1/9$): [10.17605/OSF.IO/NM5BW](https://doi.org/10.17605/OSF.IO/NM5BW)
* **Quantitative Archival Reproducibility Digest (SHA-256):**
  `de998b2fa755830ef45864614776390eba573aa1323c61a9d8c65f1c0688bd18`

## Overview & Six-Layer Verification Architecture

This archive contains **14 files at the root directory** (13 release files verified by `MANIFEST_SHA256.txt`, plus `MANIFEST_SHA256.txt` itself) for **Appendix B v1.0: Bounded Self-Referential Dynamics in a 4D Soft-Wall ODE System**, a dimensionless, unit-normalized 4D dynamical model $x(t) = (\gamma(t), H_{\text{obs}}(t), \dot{I}(t), \Psi(t))$ with structural constants $a=b=c=f=g=\mathcal{E}_{\text{avail}}=1$:

1. **Globally Smooth ($C^\infty$) Quadratic Self-Feedback Ratio (`R_self_sq`) & Equilibrium Decoupling:**
   Every component (`ExistenceVerification.lean`, `verification_existence_v1_0.py`, `dashboard.html`, `launch_existence_sandbox.py`, and `theory_of_existence_v1_0.tex`) uses the smooth quadratic ratio
   $$R_{\text{self}}(s) = \frac{\text{fb}(s)^2}{\text{fb}(s)^2 + (f\mathcal{P}_0)^2 + \sigma_0^2}, \qquad \text{fb}(s) = e\tanh(\dot{I}/\epsilon_I)\,\dot{H}_{\text{obs}}, \quad \sigma_0 = 0.05.$$
   At any steady state ($\dot{H}_{\text{obs}}^* = 0$), both $R_{\text{self}}^* = 0$ and $\nabla R_{\text{self}}|_ * = (0,0,0,0)$ vanish by quadratic construction, so row 4 of the $4\times 4$ equilibrium Jacobian $J_4$ is $[0, 0, 0, -\mu_\Psi]$ and $\lambda_{\text{lucid}}$ decouples from linear fixed-point stability. Only the product $\Lambda_{\text{eff}} = \lambda_{\text{sim}}r_{\text{LRSC}}$ enters the ODEs (fixing $r_{\text{LRSC}} = 1/9$ sets the scale where $\lambda_{\text{sim}} = 7.20$ corresponds to $\Lambda_{\text{eff}} = 0.80$).
2. **Closed-Form $e_c(\mathcal{P}_0)$ Domain ($\mathcal{P}_0 > 10/99$) & Global Monotonicity ($de_c/d\mathcal{P}_0 > 0$):**
   In the $\epsilon_H, \epsilon_I \to 0^+$ limit on the unregularized interior positivity domain $\mathcal{P}_0 > \mathcal{P}_{\text{crit}} = 10/99 = 0.\overline{10} \approx 0.10101$ (where $H_{\text{obs},\text{v1.4}}^* = 100 - (10+\mathcal{P}_0)/\mathcal{P}_0 > 0$), the open-loop Hopf threshold reduces to the rational function $e_c(\mathcal{P}_0) = \mathcal{P}_0(\mathcal{P}_0+20)(11\mathcal{P}_0+10)(\mathcal{P}_0^2+120\mathcal{P}_0+100)/[100(\mathcal{P}_0+10)(\mathcal{P}_0^2+130\mathcal{P}_0+200)]$ with strictly positive derivative $de_c/d\mathcal{P}_0 > 0$ (all 8 primitive numerator coefficients `[11, 3085, 250000, 5542000, 49530000, 120550000, 94000000, 20000000]` and 7 denominator coefficients strictly positive). For $\mathcal{P}_0 > 0$, the closed-loop threshold is given by the Routh–Hurwitz locus $e_c^{\text{CL}}(\mathcal{P}_0, \Lambda_{\text{eff}})$.
3. **Zero-Drive ($\mathcal{P}_0 = 0$) $\Psi$-Memory Persistence, $2\times 2$ Ablation & Saddle-Node of Limit Cycles (SNLC):**
   At $\mathcal{P}_0 = 0$ and $0 \le \Lambda_{\text{eff}} < 1$, the equilibrium quadratic $\Lambda_{\text{eff}}(\gamma^*)^2 + 10(1-\Lambda_{\text{eff}})\gamma^* = 0$ has roots $\gamma^* = 0$ and $\gamma^* = -10(1-\Lambda_{\text{eff}})/\Lambda_{\text{eff}} < 0$ ($-2.5 < 0$ at $\Lambda_{\text{eff}} = 0.80$), so **no positive interior equilibrium or Hopf bifurcation exists** (`N/A (SNLC)`), and the origin $(0,0,0,0)$ is locally asymptotically stable (eigenvalues $\{-0.105573, -1.894427, -200.0, -0.60\}$ at $\Lambda_{\text{eff}} = 0.80$). A $2\times 2$ ablation of $(\lambda_{\text{sim}}, \lambda_{\text{lucid}})$ at $\mathcal{P}_0 = 0$ shows:
   * **Memoryless Proportional Loop Alone Collapses ($\lambda_{\text{sim}} = 7.20, \lambda_{\text{lucid}} = 0$):** Because $\gamma(t)$ drops to $\gamma_{\min} \approx 3.6 \times 10^{-4}$ during the deep erasure trough, $\Lambda_{\text{eff}}\gamma(t)$ vanishes at the trough and the trajectory decays to $(0,0,0,0)$ ($\|x(400)\|_2 = 1.84 \times 10^{-16} < 10^{-12}$).
   * **Slow $\Psi$-Memory Reservoir Bridges the Trough ($\lambda_{\text{sim}} = 0, \lambda_{\text{lucid}} = 0.30$, Regime IVc):** Because $\Psi(t)$ decays slowly at rate $\mu_\Psi = 0.60$, it stays above $\Psi_{\min} = 0.7305$ across the trough so that $\lambda_{\text{lucid}}\Psi(t) \ge 0.219$ acts as a stored capacitive drive that sustains the periodic cycle ($T = 17.936755$, $\overline{\mathcal{P}}_{\text{eff},T} = 0.239590$, $\max \dot{I} = 265.85$).
   * **Saddle-Node of Limit Cycles (SNLC) Continuation & Synergistic Flux Amplification (Regime IVb):** Bisection in $\lambda_{\text{lucid}}$ at $\mathcal{P}_0 = 0$ locates the SNLC birth of the periodic cycle at $\lambda_{\text{lucid,SNLC}}(0) = 0.080341$ when $\Lambda_{\text{eff}} = 0$. Coupling $\Lambda_{\text{eff}} = 0.80$ ($\lambda_{\text{sim}} = 7.20$, Regime IVb) amplifies the recovery stroke, boosting the period-averaged flux by $9.14\times$ to $\overline{\mathcal{P}}_{\text{eff},T} = 2.190321$ ($T = 19.387103$) and lowering the critical SNLC threshold by $27.8\%$ to $\lambda_{\text{lucid,SNLC}}(0.80) = 0.057974$.
4. **20D Variational Monodromy Matrix $\Phi(T)$ & Saltation-Projected Transverse Floquet Multipliers:**
   Layer 4 of `verification_existence_v1_0.py` integrates the 20D first-variation monodromy equations $\dot{x} = F_{\text{exist}}(x), \dot{\Phi} = J_4(x)\Phi(t)$ ($\Phi(0) = I_4$, with `jacobian_4d` verified against `SymPy` to $1.65 \times 10^{-17}$) over one period $T$ on $\Sigma = \{\gamma = 2.0, \dot{\gamma} > 0\}$. Across Regimes II, IIIb, IVa, IVb, and IVc, the monodromy matrix $M = \Phi(T)$ satisfies the tangent eigenvector identity $\|M F(x_\Sigma) - F(x_\Sigma)\|_2 / \|F(x_\Sigma)\|_2 \le 1.72 \times 10^{-9}$ (confirming $\mu_{\text{tangent}} = 1$), and the saltation-projected $3\times 3$ return-map Jacobian $DP_{\text{var}} = (P_\Sigma \Phi(T))_{1:4,1:4}$ yields $\rho_{\text{var}}(DP_{2,H\dot{I}}) \in [5.03 \times 10^{-12}, 6.64 \times 10^{-12}] < 10^{-9}$, while in Regime II $(DP_{\text{var}})_{3,3}$ matches the analytical orbital integral $\mu_{\Psi,\text{exact}} = \exp(-\int_0^T (\kappa_{\text{meta}}R_{\text{self}} + \mu_\Psi)\,dt) = 2.537698 \times 10^{-6}$ to $< 10^{-11}$.
5. **Lean 4 Pointwise Boundary & Linear-Floor Algebra (`ExistenceVerification.lean`, 17 Theorems):**
   `lakefile.lean` registers both `BoundaryVerification` and `ExistenceVerification` as default library targets (`0 sorry`, `0 nlinarith`, `#print axioms` audit lines included). `ExistenceVerification.lean` machine-checks 17 pointwise algebraic theorems over Mathlib's $\mathbb{R}$: `R_self_sq_nonneg`, the upper and lower entropy-rate bounds `dHobs_dt_4D_le_M_H` and `dHobs_dt_4D_ge_neg_M_H` ($-M_H \le dH_{\text{obs}}/dt \le M_H$), the pointwise linear comparison lower bounds `dgamma_dt_4D_ge_neg_decay` ($d\gamma/dt \ge -b\gamma$) and `dIdot_dt_4D_ge_linear_floor` ($d\dot{I}/dt \ge f\mathcal{P}_0 - (g + e M_H/\epsilon_I)\dot{I}$), the eight boundary-face theorems (including `dIdot_dt_4D_nonneg_at_boundary_zero_drive` for $\mathcal{P}_0 \ge 0$), and the unified 8-face certificate `existence_box_8face_nagumo_certificate`. Continuous-time trajectory invariance $x(t) \in \mathcal{B}_4$, Grönwall positivity $\gamma(t) \ge \gamma(0)e^{-bt} > 0$, and uniform positive floors $\liminf_{t\to\infty}\dot{I}(t) \ge \dot{I}_{\text{floor}}(\mathcal{P}_0) > 0$ and $\liminf_{t\to\infty}\gamma(t) \ge \gamma_{\text{floor}}(\mathcal{P}_0) > 0$ (for $\mathcal{P}_0 > 0$) follow by classical Nagumo and ODE comparison theorems.

## Quick Verification Commands

```bash
# 1. Verify SHA-256 checksums of all 13 payload files listed in MANIFEST_SHA256.txt:
sha256sum -c MANIFEST_SHA256.txt

# 2. Build both Lean 4 libraries (BoundaryVerification and ExistenceVerification):
lake exe cache get
lake build

# 3. Run the 6-layer symbolic, Poincaré, 20D variational Floquet, SNLC & Min/Max verification suite:
pip install -r requirements.txt
python3 verification_existence_v1_0.py

# 4. Launch the interactive 4D ODE simulator locally:
python3 launch_existence_sandbox.py
```
