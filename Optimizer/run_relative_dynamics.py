"""
Base script to propagate and plot relative orbital-element dynamics.

Inputs:
- Chief classical elements
- Initial relative state:
    [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]
- Time span
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from propagator_42.runner import ClassicalElements as ClassicalElements42
from propagator_42.runner import propagate_relative_hill_42
from relative_dynamics_j2 import J2OrbitParams, propagate_relative_dynamics_j2


# Earth constants (SI)
MU_EARTH = 3.986004418e14  # m^3/s^2
J2_EARTH = 1.08262668e-3
R_EARTH = 6378137.0  # m


@dataclass(frozen=True)
class ClassicalElements:
    """Classical orbital elements."""

    a: float  # semi-major axis [m]
    e: float  # eccentricity [-]
    i: float  # inclination [rad]
    omega: float  # argument of periapsis [rad]
    Omega: float  # RAAN [rad]
    M: float  # mean anomaly [rad]


def build_j2_params(elements: ClassicalElements) -> J2OrbitParams:
    """Compute secular J2 model terms from classical elements."""
    eta = np.sqrt(1.0 - elements.e**2)
    n = np.sqrt(MU_EARTH / elements.a**3)
    p = 3.0 * np.cos(elements.i) ** 2
    q = 5.0 * np.cos(elements.i) ** 2 - 1.0
    k = 3.0 * n * J2_EARTH * R_EARTH**2 / (4.0 * elements.a**2 * eta**4)
    return J2OrbitParams(e=elements.e, i=elements.i, eta=eta, p=p, q=q, k=k, omega0=elements.omega)


def reconstruct_deputy_from_relative(
    chief: ClassicalElements,
    delta_alpha0: np.ndarray,
) -> ClassicalElements:
    """
    Reconstruct deputy COEs from chief COEs and initial relative elements.

    This uses the relative-element definitions:
    delta_ex = ed*cos(omegad) - ec*cos(omegac)
    delta_ey = ed*sin(omegad) - ec*sin(omegac)
    delta_a  = ad/ac - 1
    delta_ix = id - ic
    delta_iy = (Omegad - Omegac) * sin(ic)
    """
    delta_a, _delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy = delta_alpha0

    ex_d = delta_ex + chief.e * np.cos(chief.omega)
    ey_d = delta_ey + chief.e * np.sin(chief.omega)
    e_d = float(np.hypot(ex_d, ey_d))
    omega_d = float(np.arctan2(ey_d, ex_d))

    a_d = chief.a * (1.0 + delta_a)
    i_d = chief.i + delta_ix
    Omega_d = chief.Omega + delta_iy / np.sin(chief.i)

    # M is not needed by the current J2 RHS, so keep chief M as placeholder.
    return ClassicalElements(a=a_d, e=e_d, i=i_d, omega=omega_d, Omega=Omega_d, M=chief.M)


def propagate_from_chief_and_relative(
    chief: ClassicalElements,
    delta_alpha0: np.ndarray,
    t_final: float,
    output_dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build chief/deputy J2 parameters and propagate relative dynamics.
    """
    deputy = reconstruct_deputy_from_relative(chief, delta_alpha0)
    chief_params = build_j2_params(chief)
    deputy_params = build_j2_params(deputy)
    return propagate_relative_dynamics_j2(
        initial_delta_alpha=delta_alpha0,
        chief=chief_params,
        deputy=deputy_params,
        t_final=t_final,
        dt=output_dt,
        method="RK45",
        rtol=1e-12,
        atol=1e-12,
    )


def chief_argument_of_latitude_history(chief: ClassicalElements, times: np.ndarray) -> np.ndarray:
    """
    Compute chief argument of latitude history u_c(t) = omega_c(t) + M_c(t).

    Uses secular J2 drift with linear-in-time omega and M rates.
    """
    chief_params = build_j2_params(chief)
    n = np.sqrt(MU_EARTH / chief.a**3)

    # Secular rates under J2 model
    omega_dot = chief_params.k * chief_params.q
    m_dot = n + chief_params.eta * chief_params.p * chief_params.k

    omega_t = chief.omega + omega_dot * times
    m_t = chief.M + m_dot * times
    return omega_t + m_t


def relative_position_hill_approx(
    chief: ClassicalElements,
    states: np.ndarray,
    u_c: np.ndarray,
) -> np.ndarray:
    """
    Approximate deputy relative position in chief Hill frame from ROE states.

    r_approx = a_c * [
        roe1 - roe3*cos(u_c) - roe4*sin(u_c),
        roe2 + 2*roe3*sin(u_c) - 2*roe4*cos(u_c),
        roe5*sin(u_c) - roe6*cos(u_c)
    ]
    """
    cos_u = np.cos(u_c)
    sin_u = np.sin(u_c)

    x = chief.a * (states[:, 0] - states[:, 2] * cos_u - states[:, 3] * sin_u)
    y = chief.a * (states[:, 1] + 2.0 * states[:, 2] * sin_u - 2.0 * states[:, 3] * cos_u)
    z = chief.a * (states[:, 4] * sin_u - states[:, 5] * cos_u)
    return np.column_stack((x, y, z))


def plot_relative_history(times: np.ndarray, states: np.ndarray) -> plt.Figure:
    """Plot each relative state component versus time."""
    labels = [
        "delta_a [-]",
        "delta_lambda [rad]",
        "delta_ex [-]",
        "delta_ey [-]",
        "delta_ix [rad]",
        "delta_iy [rad]",
    ]

    t_hours = times / 3600.0
    fig, axes = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    axes = axes.ravel()
    for idx, label in enumerate(labels):
        axes[idx].plot(t_hours, states[:, idx], linewidth=1.5)
        axes[idx].set_ylabel(label)
        axes[idx].grid(True, alpha=0.3)
    axes[-2].set_xlabel("time [hours]")
    axes[-1].set_xlabel("time [hours]")
    fig.suptitle("Relative Dynamics Under J2 Secular Drift")
    fig.tight_layout()
    return fig


def plot_hill_trajectory_3d(r_hill: np.ndarray) -> plt.Figure:
    """Plot approximate relative trajectory in chief Hill frame."""
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(r_hill[:, 0], r_hill[:, 1], r_hill[:, 2], color="k", linewidth=0.8)
    ax.scatter(r_hill[0, 0], r_hill[0, 1], r_hill[0, 2], s=30, label="start")
    ax.scatter(r_hill[-1, 0], r_hill[-1, 1], r_hill[-1, 2], s=30, label="end")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("T [m]")
    ax.set_zlabel("N [m]")
    ax.set_title("Approximate Relative Trajectory in Chief Hill Frame")
    ax.legend(loc="best")
    return fig


def plot_hill_comparison(times: np.ndarray, r_approx: np.ndarray, r_42: np.ndarray) -> plt.Figure:
    """Overlay approximate and 42 trajectories and show component errors."""
    t_hours = times / 3600.0
    err = r_approx - r_42

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax3d = fig.add_subplot(2, 2, 1, projection="3d")
    ax3d.plot(r_approx[:, 0], r_approx[:, 1], r_approx[:, 2], color="k", linewidth=0.8, label="ROE approx")
    ax3d.plot(r_42[:, 0], r_42[:, 1], r_42[:, 2], color="tab:blue", linewidth=0.8, label="42 detailed")
    ax3d.set_xlabel("R [m]")
    ax3d.set_ylabel("T [m]")
    ax3d.set_zlabel("N [m]")
    ax3d.set_title("Hill Trajectory Overlay")
    ax3d.legend(loc="best")

    labels = ["R error [m]", "T error [m]", "N error [m]"]
    for idx, ax in enumerate([axes[0, 1], axes[1, 0], axes[1, 1]]):
        ax.plot(t_hours, err[:, idx], linewidth=1.0)
        ax.set_ylabel(labels[idx])
        ax.set_xlabel("time [hours]")
        ax.grid(True, alpha=0.3)

    fig.suptitle("ROE Approximation vs 42 Detailed Propagation")
    fig.tight_layout()
    return fig


def save_figure_outputs(fig: plt.Figure, output_stem: str = "relative_dynamics") -> Path:
    """
    Save quick-look PNG image.

    Returns:
        png_path
    """
    figures_dir = Path(__file__).resolve().parent / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    png_path = figures_dir / f"{output_stem}.png"

    fig.savefig(png_path, dpi=160)
    return png_path


def main() -> None:
    # Example input block: replace these values as needed.
    chief = ClassicalElements(
        a=7000e3, # value in meters
        e=0.0010,
        i=np.deg2rad(55.0),
        omega=np.deg2rad(40.0),
        Omega=np.deg2rad(20.0),
        M=np.deg2rad(10.0),
    )
    # INPUT RELATIVE ELEMENTS, VALUES ARE IN RADIANS OR NONDIMENSIONAL UNITS
    # delta_a = a_d/a_c - 1 -> Should be 0 if no relative motion is desired
    # delta_lambda = (M_d - M_c) + (omega_d - omega_c) + (Omega_d - Omega_c) * cos(i_c)
    # delta_ex = e_d * cos(omega_d) - e_c * cos(omega_c)
    # delta_ey = e_d * sin(omega_d) - e_c * sin(omega_c)
    # delta_ix = i_d - i_c
    # delta_iy = (Omega_d - Omega_c) * sin(i_c)
    delta_alpha0 = np.array([0, 2e-6, -3e-7, 4e-7, -2e-6, 1e-6], dtype=float)

    t_final = 30.0 * 86400.0  # 30 days [s]
    output_dt = 60.0  # output every minute [s]

    times, states = propagate_from_chief_and_relative(
        chief=chief,
        delta_alpha0=delta_alpha0,
        t_final=t_final,
        output_dt=output_dt,
    )
    u_c = chief_argument_of_latitude_history(chief, times)
    r_hill = relative_position_hill_approx(chief, states, u_c)
    chief_42 = ClassicalElements42(
        a=chief.a,
        e=chief.e,
        i=chief.i,
        omega=chief.omega,
        Omega=chief.Omega,
        M=chief.M,
    )
    r_hill_42 = propagate_relative_hill_42(chief_42, delta_alpha0, times)

    fig = plot_relative_history(times, states)
    fig_hill = plot_hill_trajectory_3d(r_hill)
    fig_compare = plot_hill_comparison(times, r_hill, r_hill_42)

    png_path_states = save_figure_outputs(fig, output_stem="relative_dynamics")
    png_path_hill = save_figure_outputs(fig_hill, output_stem="relative_hill_trajectory")
    png_path_compare = save_figure_outputs(fig_compare, output_stem="relative_hill_comparison_42")
    print(f"Saved quick-look image: {png_path_states}")
    print(f"Saved 3D Hill-frame trajectory: {png_path_hill}")
    print(f"Saved comparison vs 42 detailed: {png_path_compare}")
    plt.show()


if __name__ == "__main__":
    main()
