"""ROE propagation helpers used by orbit validation scripts."""

from __future__ import annotations

import numpy as np

from j2_roe_core import (
    ClassicalElements,
    chief_secular_classical_history,
    reconstruct_deputy_from_relative,
    true_anomaly_from_mean,
)
from relative_dynamics_j2 import propagate_relative_dynamics_j2


def propagate_from_chief_and_relative(
    chief: ClassicalElements,
    delta_alpha0: np.ndarray,
    t_final: float,
    output_dt: float,
    *,
    chief_mean_anomaly_rate: str = "j2_averaged",
) -> tuple[np.ndarray, np.ndarray, ClassicalElements, np.ndarray]:
    delta0 = np.asarray(delta_alpha0, dtype=float).reshape(6)
    times, states = propagate_relative_dynamics_j2(
        delta0,
        chief,
        t_final,
        output_dt,
        chief_mean_anomaly_rate=chief_mean_anomaly_rate,
    )
    return times, states, chief, delta0


def roe_deputy_classical_history(
    chief: ClassicalElements,
    times: np.ndarray,
    states: np.ndarray,
    *,
    chief_mean_anomaly_rate: str = "j2_averaged",
) -> dict[str, np.ndarray]:
    chief_hist = chief_secular_classical_history(
        chief, times, chief_mean_anomaly_rate=chief_mean_anomaly_rate
    )
    n = times.size
    sma = np.zeros(n)
    ecc = np.zeros(n)
    inc = np.zeros(n)
    raan = np.zeros(n)
    argp = np.zeros(n)
    M_arr = np.zeros(n)
    nu = np.zeros(n)
    for k in range(n):
        c = ClassicalElements(
            float(chief_hist["sma_m"][k]),
            float(chief_hist["ecc"][k]),
            float(chief_hist["inc_rad"][k]),
            float(chief_hist["argp_rad"][k]),
            float(chief_hist["raan_rad"][k]),
            float(chief_hist["M_rad"][k]),
        )
        d = reconstruct_deputy_from_relative(c, states[k])
        sma[k] = d.a
        ecc[k] = d.e
        inc[k] = d.i
        raan[k] = d.Omega
        argp[k] = d.omega
        M_arr[k] = d.M
        nu[k] = float(true_anomaly_from_mean(np.array([d.M]), d.e)[0])
    return {
        "sma_m": sma,
        "ecc": ecc,
        "inc_rad": inc,
        "raan_rad": raan,
        "argp_rad": argp,
        "nu_rad": nu,
        "M_rad": M_arr,
    }


def roe_classical_element_histories(
    chief: ClassicalElements,
    times: np.ndarray,
    states: np.ndarray,
    *,
    chief_mean_anomaly_rate: str = "j2_averaged",
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    chief_hist = chief_secular_classical_history(
        chief, times, chief_mean_anomaly_rate=chief_mean_anomaly_rate
    )
    deputy_hist = roe_deputy_classical_history(
        chief, times, states, chief_mean_anomaly_rate=chief_mean_anomaly_rate
    )
    return chief_hist, deputy_hist
