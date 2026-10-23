#!/usr/bin/env python3
"""
Appendix B v1.0 — Interactive 4D Soft-Wall & Psi-Memory ODE Sandbox Server
Author: Prince Upadhyay

Serves `dashboard.html` on http://localhost:8085/ (and provides a `/api/regimes` JSON endpoint
computed via SciPy LSODA using the smooth C^infty quadratic ratio R_self = fb^2 / (fb^2 + (f*P_0)^2 + sigma_0^2)).
"""
import json
from pathlib import Path
import sys
import webbrowser
from http.server import SimpleHTTPRequestHandler
from socketserver import TCPServer
from threading import Timer
import numpy as np
from scipy.integrate import solve_ivp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

PARAMS = {
    "a": 1.0, "b": 1.0, "c": 1.0, "e": 3.60, "f": 1.0, "g": 1.0,
    "E_avail": 1.0, "I_erase": 10.0, "gamma_max": 10.0, "H_max": 100.0,
    "eps_H": 0.05, "eps_I": 0.05, "sigma_0": 0.05, "r_lrsc": 1.0 / 9.0, "mu_Psi": 0.60,
}


def existence_4d_field(t, y, P_0, lam_sim, lam_lucid, kappa_meta):
    gamma, H_obs, I_dot, Psi = y
    P_eff = P_0 + lam_sim * PARAMS["r_lrsc"] * gamma + lam_lucid * Psi
    dgamma = PARAMS["a"] * I_dot * PARAMS["E_avail"] * (1.0 - gamma / PARAMS["gamma_max"]) - PARAMS["b"] * gamma
    dHobs = PARAMS["c"] * gamma * (PARAMS["H_max"] - H_obs) - PARAMS["I_erase"] * np.tanh(H_obs / PARAMS["eps_H"])
    fb = PARAMS["e"] * np.tanh(I_dot / PARAMS["eps_I"]) * dHobs
    dIdot = fb + PARAMS["f"] * P_eff - PARAMS["g"] * I_dot
    R_self = (fb * fb) / (fb * fb + (PARAMS["f"] * P_0) ** 2 + PARAMS["sigma_0"] ** 2)
    dPsi = kappa_meta * R_self * (1.0 - Psi) - PARAMS["mu_Psi"] * Psi
    return [dgamma, dHobs, dIdot, dPsi]


def simulate_benchmark_regimes():
    t_eval = np.linspace(0, 160, 1200)
    y0_common = [2.05, 95.0, 4.0, 0.05]
    configs = {
        "I_waking":         (5.00, 0.00, 0.00, 0.20),
        "II_low_psi":       (4.00, 0.00, 0.00, 0.20),
        "IIIa_readout":     (4.00, 0.00, 0.00, 2.50),
        "IIIb_coupled":     (4.00, 0.00, 0.30, 2.50),
        "IVa_low_drive":    (1.00, 7.20, 0.30, 2.50),
        "IVb_zero_coupled": (0.00, 7.20, 0.30, 2.50),
        "IVc_pure_psi_mem": (0.00, 0.00, 0.30, 2.50),
    }
    out = {"t": np.round(t_eval, 3).tolist()}
    for key, (P0, lsim, lluc, kap) in configs.items():
        sol = solve_ivp(
            existence_4d_field, (0, 160), y0_common,
            args=(P0, lsim, lluc, kap), t_eval=t_eval,
            method="LSODA", rtol=1e-9, atol=1e-12
        )
        out[key] = {
            "gamma": np.round(sol.y[0], 4).tolist(),
            "H_obs": np.round(sol.y[1], 4).tolist(),
            "I_dot": np.round(sol.y[2], 4).tolist(),
            "Psi":   np.round(sol.y[3], 4).tolist(),
        }
    return out


class ExistenceDashboardHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        base_dir = Path(__file__).resolve().parent
        if self.path in ("/", "/index.html", "/dashboard.html"):
            dash_path = base_dir / "dashboard.html"
            self.send_response(200)
            self.send_header("Content-type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(dash_path.read_bytes())
        elif self.path == "/api/regimes":
            self.send_response(200)
            self.send_header("Content-type", "application/json; charset=utf-8")
            self.end_headers()
            self.wfile.write(json.dumps(simulate_benchmark_regimes()).encode("utf-8"))
        else:
            self.send_error(404, "File Not Found")


if __name__ == "__main__":
    print("Launching 4D Soft-Wall & Psi-Memory Sandbox on http://localhost:8085/")
    Timer(1.0, lambda: webbrowser.open_new("http://localhost:8085/")).start()
    TCPServer.allow_reuse_address = True
    with TCPServer(("", 8085), ExistenceDashboardHandler) as server:
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nShutting down sandbox safely.")
