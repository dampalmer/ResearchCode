"""
Numerical propagation of chief/deputy relative dynamics with J2 secular drift.

Propagates delta-alpha = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy].
RHS uses instantaneous deputy (a, e, i, omega) from the current state so mean-longitude
rates track evolving eccentricity; chief geometry uses the same first-order secular model.
"""

from __future__ import annotations

from typing import Iterable, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from j2_roe_core import (
    ClassicalElements,
    build_j2_params,
    chief_secular_at_time,
    reconstruct_deputy_from_relative,
)

# Re-export for callers that imported J2OrbitParams from this module
from j2_roe_core import J2OrbitParams  # noqa: F401


def _rhs_delta_alpha(
    t: float,
    state: np.ndarray,
    chief_epoch: ClassicalElements,
    *,
    chief_mean_anomaly_rate: str,
) -> np.ndarray:
    chief_t = chief_secular_at_time(chief_epoch, t, chief_mean_anomaly_rate=chief_mean_anomaly_rate)
    deputy_t = reconstruct_deputy_from_relative(chief_t, state)
    chief_p = build_j2_params(chief_t)
    deputy_p = build_j2_params(deputy_t)
    omega_c = chief_t.omega
    omega_d = deputy_t.omega
    i_c = chief_t.i

    d_delta_a = 0.0
    # J2 secular delta_lam_dot only. Keplerian (n_d - n_c) along-track coupling enters the
    # position map via -3/2 delta_a u; including it here double-counts that drift.
    d_delta_lambda = (deputy_p.eta * deputy_p.p * deputy_p.k - chief_p.eta * chief_p.p * chief_p.k) + (
        deputy_p.k * deputy_p.q - chief_p.k * chief_p.q
    ) * np.cos(i_c)
    d_delta_ex = -deputy_p.e * deputy_p.k * deputy_p.q * np.sin(omega_d) + chief_p.e * chief_p.k * chief_p.q * np.sin(
        omega_c
    )
    d_delta_ey = deputy_p.e * deputy_p.k * deputy_p.q * np.cos(omega_d) - chief_p.e * chief_p.k * chief_p.q * np.cos(
        omega_c
    )
    d_delta_ix = 0.0
    d_delta_iy = -2.0 * (deputy_p.k * np.cos(deputy_p.i) - chief_p.k * np.cos(chief_p.i)) * np.sin(chief_p.i)

    return np.array(
        [d_delta_a, d_delta_lambda, d_delta_ex, d_delta_ey, d_delta_ix, d_delta_iy],
        dtype=float,
    )


def propagate_relative_dynamics_j2(
    initial_delta_alpha: Iterable[float],
    chief_epoch: ClassicalElements,
    t_final: float,
    dt: float,
    *,
    chief_mean_anomaly_rate: str = "j2_averaged",
    method: str = "RK45",
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Propagate relative dynamics under J2 secular model.

    Args:
        initial_delta_alpha: [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy].
        chief_epoch: Mean-line chief at t=0.
        t_final: Final time [s].
        dt: Output step [s].
        chief_mean_anomaly_rate: ``j2_averaged`` (default) uses n + eta*P*k so chief
            mean anomaly carries the same first-order J2 secular rate as omega/Omega.
            ``kepler`` uses sqrt(mu/a^3) only. Legacy alias ``kepler_42`` is accepted.
            Relative EOMs and chief omega/Omega always include J2 regardless of mode.
    """
    if t_final < 0.0:
        raise ValueError("t_final must be non-negative.")
    if dt <= 0.0:
        raise ValueError("dt must be positive.")

    y0 = np.asarray(initial_delta_alpha, dtype=float)
    if y0.shape != (6,):
        raise ValueError("initial_delta_alpha must be length 6.")

    num_steps = int(np.ceil(t_final / dt))
    times = np.linspace(0.0, t_final, num_steps + 1, dtype=float)

    solution = solve_ivp(
        fun=lambda t, y: _rhs_delta_alpha(
            t, y, chief_epoch, chief_mean_anomaly_rate=chief_mean_anomaly_rate
        ),
        t_span=(0.0, t_final),
        y0=y0,
        method=method,
        t_eval=times,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"J2 relative dynamics integration failed: {solution.message}")

    return solution.t, solution.y.T
