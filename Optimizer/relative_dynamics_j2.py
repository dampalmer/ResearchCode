"""
Numerical propagation of chief/deputy relative dynamics with J2 secular drift.

This module propagates relative elements directly:
    delta-alpha = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class J2OrbitParams:
    """Per-orbit constants for the secular J2 relative dynamics model."""

    e: float  # eccentricity
    i: float  # inclination [rad]
    eta: float  # sqrt(1 - e^2)
    p: float  # 3 cos^2(i)
    q: float  # 5 cos^2(i) - 1
    k: float  # 3 n J2 Re^2 / (4 a^2 eta^4)
    omega0: float  # argument of periapsis at t=0 [rad]

    @property
    def omega_dot(self) -> float:
        """Secular argument-of-periapsis rate under this model."""
        return self.k * self.q


def _rhs_delta_alpha(
    t: float,
    state: np.ndarray,
    chief: J2OrbitParams,
    deputy: J2OrbitParams,
) -> np.ndarray:
    """
    Time derivative of delta-alpha using the provided J2 secular model.
    """
    _ = state  # Model is explicit in time; state kept for ODE interface.
    omega_c_t = chief.omega0 + chief.omega_dot * t
    omega_d_t = deputy.omega0 + deputy.omega_dot * t

    d_delta_a = 0.0
    d_delta_lambda = (deputy.eta * deputy.p * deputy.k - chief.eta * chief.p * chief.k) + (
        deputy.k * deputy.q - chief.k * chief.q
    ) * np.cos(chief.i)
    d_delta_ex = -deputy.e * deputy.k * deputy.q * np.sin(omega_d_t) + chief.e * chief.k * chief.q * np.sin(
        omega_c_t
    )
    d_delta_ey = deputy.e * deputy.k * deputy.q * np.cos(omega_d_t) - chief.e * chief.k * chief.q * np.cos(
        omega_c_t
    )
    d_delta_ix = 0.0
    d_delta_iy = -2.0 * (deputy.k * np.cos(deputy.i) - chief.k * np.cos(chief.i)) * np.sin(chief.i)

    return np.array(
        [d_delta_a, d_delta_lambda, d_delta_ex, d_delta_ey, d_delta_ix, d_delta_iy],
        dtype=float,
    )


def propagate_relative_dynamics_j2(
    initial_delta_alpha: Iterable[float],
    chief: J2OrbitParams,
    deputy: J2OrbitParams,
    t_final: float,
    dt: float,
    *,
    method: str = "RK45",
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Numerically propagate relative dynamics for chief/deputy under J2 secular drift.

    Args:
        initial_delta_alpha: Initial relative state (length 6):
            [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy].
        chief: Chief J2 secular parameters.
        deputy: Deputy J2 secular parameters.
        t_final: Final propagation time [s], must be >= 0.
        dt: Output sampling interval [s], must be > 0.
        method: solve_ivp method (default: "RK45").
        rtol: Relative tolerance passed to solve_ivp.
        atol: Absolute tolerance passed to solve_ivp.

    Returns:
        times: np.ndarray with shape (N,)
        states: np.ndarray with shape (N, 6) for
            [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]
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
        fun=lambda t, y: _rhs_delta_alpha(t, y, chief, deputy),
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

