from __future__ import annotations

import ctypes
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

MU_EARTH = 3.986004418e14  # m^3/s^2
J2_EARTH = 1.08262668e-3
R_EARTH = 6378137.0  # m


@dataclass(frozen=True)
class ClassicalElements:
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M: float


class OrbitElementsJ2(ctypes.Structure):
    _fields_ = [
        ("mu", ctypes.c_double),
        ("eccentricity", ctypes.c_double),
        ("inclination", ctypes.c_double),
        ("right_ascension_epoch", ctypes.c_double),
        ("argument_of_perigee_epoch", ctypes.c_double),
        ("mean_anomaly_epoch", ctypes.c_double),
        ("mean_motion", ctypes.c_double),
        ("mean_semi_major_axis", ctypes.c_double),
        ("epoch", ctypes.c_double),
        ("raan_rate", ctypes.c_double),
        ("arg_perigee_rate", ctypes.c_double),
        ("j2_rw2_over_a", ctypes.c_double),
    ]


class OrbitState(ctypes.Structure):
    _fields_ = [
        ("position", ctypes.c_double * 3),
        ("velocity", ctypes.c_double * 3),
        ("true_anomaly", ctypes.c_double),
        ("semi_major_axis", ctypes.c_double),
        ("semi_latus_rectum", ctypes.c_double),
        ("right_ascension", ctypes.c_double),
        ("argument_of_perigee", ctypes.c_double),
    ]


def reconstruct_deputy_from_relative(chief: ClassicalElements, delta_alpha0: np.ndarray) -> ClassicalElements:
    delta_a, _delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy = delta_alpha0

    ex_d = delta_ex + chief.e * np.cos(chief.omega)
    ey_d = delta_ey + chief.e * np.sin(chief.omega)
    e_d = float(np.hypot(ex_d, ey_d))
    omega_d = float(np.arctan2(ey_d, ex_d))
    a_d = chief.a * (1.0 + delta_a)
    i_d = chief.i + delta_ix
    Omega_d = chief.Omega + delta_iy / np.sin(chief.i)

    M_d = (
        chief.M
        + delta_alpha0[1]
        - (omega_d - chief.omega)
        - (Omega_d - chief.Omega) * np.cos(chief.i)
    )
    return ClassicalElements(a=a_d, e=e_d, i=i_d, omega=omega_d, Omega=Omega_d, M=float(M_d))


def _build_propagator_lib(lib_path: Path, repo_root: Path) -> None:
    src_dir = repo_root / "SOCRATES" / "propagators" / "src"
    include_dir = repo_root / "SOCRATES" / "propagators" / "include"
    cmd = [
        "gcc",
        "-shared",
        "-fPIC",
        "-O2",
        "-I",
        str(include_dir),
        str(src_dir / "vector_math.c"),
        str(src_dir / "orbit_propagator.c"),
        "-lm",
        "-o",
        str(lib_path),
    ]
    subprocess.run(cmd, check=True, cwd=str(repo_root))


def _load_lib() -> ctypes.CDLL:
    here = Path(__file__).resolve().parent
    repo_root = here.parents[1]
    lib_path = here / "libpropagators.so"
    if not lib_path.exists():
        _build_propagator_lib(lib_path, repo_root)

    lib = ctypes.CDLL(str(lib_path))
    lib.propagate_orbit_j2.argtypes = [
        ctypes.POINTER(OrbitElementsJ2),
        ctypes.c_double,
        ctypes.POINTER(OrbitState),
    ]
    lib.propagate_orbit_j2.restype = None
    return lib


def _j2_elements_from_classical(elements: ClassicalElements) -> OrbitElementsJ2:
    eta = np.sqrt(1.0 - elements.e**2)
    n = np.sqrt(MU_EARTH / elements.a**3)
    p = 3.0 * np.cos(elements.i) ** 2
    q = 5.0 * np.cos(elements.i) ** 2 - 1.0
    k = 3.0 * n * J2_EARTH * R_EARTH**2 / (4.0 * elements.a**2 * eta**4)

    raan_rate = -2.0 * k * np.cos(elements.i)
    argp_rate = k * q
    mean_motion = n + eta * p * k

    return OrbitElementsJ2(
        mu=MU_EARTH,
        eccentricity=elements.e,
        inclination=elements.i,
        right_ascension_epoch=elements.Omega,
        argument_of_perigee_epoch=elements.omega,
        mean_anomaly_epoch=elements.M,
        mean_motion=mean_motion,
        mean_semi_major_axis=elements.a,
        epoch=0.0,
        raan_rate=raan_rate,
        arg_perigee_rate=argp_rate,
        j2_rw2_over_a=J2_EARTH * R_EARTH**2 / elements.a,
    )


def _state_to_arrays(state: OrbitState) -> tuple[np.ndarray, np.ndarray]:
    r = np.array([state.position[0], state.position[1], state.position[2]], dtype=float)
    v = np.array([state.velocity[0], state.velocity[1], state.velocity[2]], dtype=float)
    return r, v


def _relative_in_chief_hill(r_c: np.ndarray, v_c: np.ndarray, r_d: np.ndarray) -> np.ndarray:
    dr = r_d - r_c
    r_hat = r_c / np.linalg.norm(r_c)
    h_vec = np.cross(r_c, v_c)
    n_hat = h_vec / np.linalg.norm(h_vec)
    t_hat = np.cross(n_hat, r_hat)
    return np.array([np.dot(dr, r_hat), np.dot(dr, t_hat), np.dot(dr, n_hat)], dtype=float)


def propagate_relative_hill_42(chief: ClassicalElements, delta_alpha0: np.ndarray, times: np.ndarray) -> np.ndarray:
    """
    Propagate chief/deputy with SOCRATES 42-aligned J2 propagator and return
    deputy relative position in chief Hill (RTN) frame for each time sample.
    """
    deputy = reconstruct_deputy_from_relative(chief, delta_alpha0)
    chief_j2 = _j2_elements_from_classical(chief)
    deputy_j2 = _j2_elements_from_classical(deputy)
    lib = _load_lib()

    r_hill = np.zeros((times.size, 3), dtype=float)
    state_c = OrbitState()
    state_d = OrbitState()
    for k, t in enumerate(times):
        lib.propagate_orbit_j2(ctypes.byref(chief_j2), float(t), ctypes.byref(state_c))
        lib.propagate_orbit_j2(ctypes.byref(deputy_j2), float(t), ctypes.byref(state_d))
        r_c, v_c = _state_to_arrays(state_c)
        r_d, _v_d = _state_to_arrays(state_d)
        r_hill[k] = _relative_in_chief_hill(r_c, v_c, r_d)
    return r_hill

