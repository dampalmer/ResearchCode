"""
Shared J2 / ROE helpers: classical elements, mean-line initialization, reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MU_EARTH = 398600.4418e9
J2_EARTH = 1.08262668e-3
R_EARTH = 6378137.0


@dataclass(frozen=True)
class ClassicalElements:
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M: float


@dataclass(frozen=True)
class J2OrbitParams:
    e: float
    i: float
    eta: float
    p: float
    q: float
    k: float
    omega0: float


def build_j2_params(elements: ClassicalElements) -> J2OrbitParams:
    eta = np.sqrt(1.0 - elements.e**2)
    n = np.sqrt(MU_EARTH / elements.a**3)
    p_coef = 3.0 * np.cos(elements.i) ** 2
    q_coef = 5.0 * np.cos(elements.i) ** 2 - 1.0
    k = 3.0 * n * J2_EARTH * R_EARTH**2 / (4.0 * elements.a**2 * eta**4)
    return J2OrbitParams(elements.e, elements.i, eta, p_coef, q_coef, k, elements.omega)


def mean_motion_kepler_from_mean_sma(mean_sma_m: float) -> float:
    return float(np.sqrt(MU_EARTH / mean_sma_m**3))


def chief_mean_anomaly_dot(chief0: ClassicalElements, *, rate_mode: str = "j2_averaged") -> float:
    """Chief secular mean-anomaly rate.

    ``j2_averaged`` (default): n + eta*P*k, matching first-order J2 secular Ṁ.
    ``kepler`` / ``kepler_42``: Keplerian n = sqrt(mu/a^3) only.
    Note: chief ω̇ and Ω̇ always include J2 via ``chief_secular_at_time``.
    """
    n = mean_motion_kepler_from_mean_sma(chief0.a)
    if rate_mode in ("kepler", "kepler_42"):
        return n
    if rate_mode == "j2_averaged":
        pp = build_j2_params(chief0)
        return n + pp.eta * pp.p * pp.k
    raise ValueError(f"Unknown chief mean-anomaly rate mode: {rate_mode!r}")


def true_anomaly_from_mean(M: np.ndarray, e: float) -> np.ndarray:
    Mw = np.mod(np.asarray(M, dtype=float) + np.pi, 2.0 * np.pi) - np.pi
    E = np.where(np.abs(e) < 0.8, Mw, np.pi)
    for _ in range(80):
        f = E - e * np.sin(E) - Mw
        fp = 1.0 - e * np.cos(E)
        d = f / fp
        E = E - d
        if float(np.max(np.abs(d))) < 1e-14:
            break
    beta = np.sqrt((1.0 + e) / (1.0 - e))
    return 2.0 * np.arctan(beta * np.tan(E / 2.0))


def _true_from_mean_scalar(M: float, e: float) -> float:
    return float(true_anomaly_from_mean(np.array([M], dtype=float), e)[0])


def mean_anomaly_from_true_scalar(nu: float, e: float) -> float:
    nu_f = float(nu)
    e_f = float(e)
    E = np.arctan2(np.sqrt(1.0 - e_f * e_f) * np.sin(nu_f), e_f + np.cos(nu_f))
    return float(E - e_f * np.sin(E))


def classical_from_osculating(
    sma_m: float,
    ecc: float,
    inc_rad: float,
    raan_rad: float,
    argp_rad: float,
    nu_rad: float,
) -> ClassicalElements:
    return ClassicalElements(
        float(sma_m),
        float(ecc),
        float(inc_rad),
        float(argp_rad),
        float(raan_rad),
        mean_anomaly_from_true_scalar(nu_rad, ecc),
    )


def osculating_kepler_to_mean_line(
    a_osc: float,
    e: float,
    i: float,
    raan: float,
    argp: float,
    nu_true: float,
    *,
    mu: float = MU_EARTH,
    j2: float = J2_EARTH,
    rw: float = R_EARTH,
) -> ClassicalElements:
    e2 = e * e
    sin2i = np.sin(i) ** 2
    sinw = np.sin(argp + nu_true)
    sin2w = sinw * sinw
    cosnu = np.cos(nu_true)
    g = ((1.0 + e * cosnu) / (1.0 - e2)) ** 3 * (1.0 - 3.0 * sin2i * sin2w)
    disc = a_osc * a_osc - 4.0 * j2 * rw * rw * g
    if disc <= 0.0:
        raise ValueError("Osculating-to-mean SMA: negative discriminant (check inputs).")
    mean_sma = 0.5 * (a_osc + np.sqrt(disc))
    E = np.arctan2(np.sqrt(1.0 - e2) * np.sin(nu_true), e + np.cos(nu_true))
    mean_anom = float(E - e * np.sin(E))
    return ClassicalElements(float(mean_sma), float(e), float(i), float(argp), float(raan), mean_anom)


def orbit_averaged_osculating_coe_hist(
    coe_hist: dict[str, np.ndarray],
    times: np.ndarray,
    *,
    window_hours: float = 2.0,
    mu: float = MU_EARTH,
) -> dict[str, np.ndarray]:
    """Rolling mean of osculating COE histories (suppresses short-period ripple)."""
    t = np.asarray(times, dtype=float).reshape(-1)
    if t.size < 2:
        raise ValueError("times must have at least two samples.")
    dt = float(np.median(np.diff(t)))
    if dt <= 0.0:
        raise ValueError("times must be strictly increasing.")

    a_ref = float(np.median(np.asarray(coe_hist["sma_m"], dtype=float)))
    period_s = 2.0 * np.pi * np.sqrt(a_ref**3 / mu)
    win_from_hours = max(int(round(window_hours * 3600.0 / dt)), 1)
    win_from_orbits = max(int(round(1.0 * period_s / dt)), 1)
    win = max(win_from_hours, win_from_orbits)

    def _roll(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float).reshape(-1)
        if win <= 1 or x.size <= win:
            return x.copy()
        kernel = np.ones(win, dtype=float) / float(win)
        full = np.convolve(x, kernel, mode="valid")
        left = np.full(win // 2, full[0])
        right = np.full(x.size - left.size - full.size, full[-1])
        return np.concatenate([left, full, right])

    raan_u = np.unwrap(np.asarray(coe_hist["raan_rad"], dtype=float))
    argp_u = np.unwrap(np.asarray(coe_hist["argp_rad"], dtype=float))
    return {
        "sma_m": _roll(np.asarray(coe_hist["sma_m"], dtype=float)),
        "ecc": _roll(np.asarray(coe_hist["ecc"], dtype=float)),
        "inc_rad": _roll(np.asarray(coe_hist["inc_rad"], dtype=float)),
        "raan_rad": _roll(raan_u),
        "argp_rad": _roll(argp_u),
        "nu_rad": np.asarray(coe_hist["nu_rad"], dtype=float).copy(),
    }


def mean_chief_for_roe_propagation(
    cyg_osc: dict[str, np.ndarray],
    times_s: np.ndarray,
    *,
    epoch_index: int = 0,
) -> ClassicalElements:
    """
    ROE chief IC: orbit-mean osculating elements → mean-line Kepler elements.

    Averages one complete forward orbit of onboard osculating samples to remove
    short-period content, then applies ``osculating_kepler_to_mean_line``
    (ROE has no short-period model).
    """
    t = np.asarray(times_s, dtype=float).reshape(-1)
    n = t.size
    if n < 2:
        raise ValueError("Need at least two onboard element samples for orbit averaging.")
    k = int(epoch_index)
    if not (0 <= k < n):
        raise ValueError(f"epoch_index {k} out of range for {n} samples.")

    a_ref = float(cyg_osc["sma_m"][k])
    period_s = 2.0 * np.pi * np.sqrt(a_ref**3 / MU_EARTH)
    dt = float(np.median(np.diff(t)))
    n_orbit = max(int(round(1.0 * period_s / dt)), 2)
    sl = slice(k, min(n, k + n_orbit))

    a_osc = float(np.mean(np.asarray(cyg_osc["sma_m"][sl], dtype=float)))
    e_osc = float(np.mean(np.asarray(cyg_osc["ecc"][sl], dtype=float)))
    i_osc = float(np.mean(np.asarray(cyg_osc["inc_rad"][sl], dtype=float)))
    raan_osc = float(np.mean(np.unwrap(np.asarray(cyg_osc["raan_rad"][sl], dtype=float))))
    argp_osc = float(np.mean(np.unwrap(np.asarray(cyg_osc["argp_rad"][sl], dtype=float))))
    nu_epoch = float(cyg_osc["nu_rad"][k])

    return osculating_kepler_to_mean_line(
        a_osc,
        e_osc,
        i_osc,
        raan_osc,
        argp_osc,
        nu_epoch,
    )


def classical_elements_to_eci(
    elements: ClassicalElements,
    *,
    mu: float = MU_EARTH,
    nu_rad: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    e = float(elements.e)
    a = float(elements.a)
    if nu_rad is None:
        nu = float(true_anomaly_from_mean(np.array([elements.M], dtype=float), e)[0])
    else:
        nu = float(nu_rad)
    i = float(elements.i)
    omega = float(elements.omega)
    Omega = float(elements.Omega)
    p = a * (1.0 - e * e)
    r_pf = p / (1.0 + e * np.cos(nu))
    r_pqw = np.array([r_pf * np.cos(nu), r_pf * np.sin(nu), 0.0])
    v_pqw = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(i), np.sin(i)
    cw, sw = np.cos(omega), np.sin(omega)
    rot = (
        np.array([[cO, sO, 0.0], [-sO, cO, 0.0], [0.0, 0.0, 1.0]])
        @ np.array([[1.0, 0.0, 0.0], [0.0, ci, si], [0.0, -si, ci]])
        @ np.array([[cw, sw, 0.0], [-sw, cw, 0.0], [0.0, 0.0, 1.0]])
    )
    return rot @ r_pqw, rot @ v_pqw


def reconstruct_deputy_from_relative(
    chief: ClassicalElements,
    delta_alpha: np.ndarray,
) -> ClassicalElements:
    delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy = (
        float(x) for x in delta_alpha[:6]
    )
    ex_d = delta_ex + chief.e * np.cos(chief.omega)
    ey_d = delta_ey + chief.e * np.sin(chief.omega)
    e_d = float(np.hypot(ex_d, ey_d))
    omega_d = float(np.arctan2(ey_d, ex_d))
    a_d = chief.a * (1.0 + delta_a)
    i_d = chief.i + delta_ix
    sc = np.sin(chief.i)
    if abs(sc) < 1e-10 and abs(delta_iy) > 1e-10:
        raise ValueError("sin(i_c)≈0 with non-zero δiy is ill-conditioned.")
    omega_d_raan = chief.Omega if abs(sc) < 1e-10 else chief.Omega + delta_iy / sc
    M_d = (
        chief.M
        + delta_lambda
        - (omega_d - chief.omega)
        - (omega_d_raan - chief.Omega) * np.cos(chief.i)
    )
    return ClassicalElements(a_d, e_d, float(i_d), omega_d, float(omega_d_raan), float(M_d))


def argument_of_latitude(
    elements: ClassicalElements,
    *,
    use_true_anomaly: bool = False,
) -> float:
    """Mean argument of latitude u = omega + M (or omega + nu for osculating samples)."""
    if use_true_anomaly:
        nu = float(true_anomaly_from_mean(np.array([elements.M], dtype=float), elements.e)[0])
        return float(elements.omega + nu)
    return float(elements.omega + elements.M)


def delta_lambda_from_deputy(
    chief: ClassicalElements,
    deputy: ClassicalElements,
    *,
    use_true_anomaly: bool = False,
) -> float:
    """delta_lambda = (u_d - u_c) + (Omega_d - Omega_c) cos(i_c)."""
    u_c = argument_of_latitude(chief, use_true_anomaly=use_true_anomaly)
    u_d = argument_of_latitude(deputy, use_true_anomaly=use_true_anomaly)
    return float((u_d - u_c) + (deputy.Omega - chief.Omega) * np.cos(chief.i))


def inverse_delta_alpha_from_deputy(
    chief: ClassicalElements,
    deputy: ClassicalElements,
    *,
    use_true_anomaly: bool = False,
) -> np.ndarray:
    delta_a = deputy.a / chief.a - 1.0
    delta_ex = deputy.e * np.cos(deputy.omega) - chief.e * np.cos(chief.omega)
    delta_ey = deputy.e * np.sin(deputy.omega) - chief.e * np.sin(chief.omega)
    delta_ix = deputy.i - chief.i
    sc = np.sin(chief.i)
    delta_iy = 0.0 if abs(sc) < 1e-10 else (deputy.Omega - chief.Omega) * sc
    delta_lambda = delta_lambda_from_deputy(
        chief, deputy, use_true_anomaly=use_true_anomaly
    )
    return np.array(
        [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy],
        dtype=float,
    )


def chief_secular_at_time(
    chief0: ClassicalElements,
    t: float,
    *,
    chief_mean_anomaly_rate: str = "j2_averaged",
) -> ClassicalElements:
    chief_params = build_j2_params(chief0)
    omega_dot = chief_params.k * chief_params.q
    m_dot = chief_mean_anomaly_dot(chief0, rate_mode=chief_mean_anomaly_rate)
    raan_dot = -2.0 * chief_params.k * np.cos(chief0.i)
    return ClassicalElements(
        chief0.a,
        chief0.e,
        chief0.i,
        float(chief0.omega + omega_dot * t),
        float(chief0.Omega + raan_dot * t),
        float(chief0.M + m_dot * t),
    )


def chief_secular_classical_history(
    chief0: ClassicalElements,
    times: np.ndarray,
    *,
    chief_mean_anomaly_rate: str = "j2_averaged",
) -> dict[str, np.ndarray]:
    chief_params = build_j2_params(chief0)
    omega_dot = chief_params.k * chief_params.q
    m_dot = chief_mean_anomaly_dot(chief0, rate_mode=chief_mean_anomaly_rate)
    raan_dot = -2.0 * chief_params.k * np.cos(chief0.i)
    times = np.asarray(times, dtype=float)
    omega_t = chief0.omega + omega_dot * times
    M_t = chief0.M + m_dot * times
    Omega_t = chief0.Omega + raan_dot * times
    nu_t = true_anomaly_from_mean(M_t, chief0.e)
    return {
        "sma_m": np.full_like(times, chief0.a, dtype=float),
        "ecc": np.full_like(times, chief0.e, dtype=float),
        "inc_rad": np.full_like(times, chief0.i, dtype=float),
        "raan_rad": Omega_t,
        "argp_rad": omega_t,
        "nu_rad": nu_t,
        "M_rad": M_t,
    }
