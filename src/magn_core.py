import numpy as np
mu_0 = 1.25663706127*1e-6 # NA^-2

def calc_r(q: float) -> float:
    return np.sqrt(q/np.pi) #m


def calc_l(N: int, r: float) -> float:
    return 2*np.pi*r*N #m


def get_U_max(arr: np.ndarray) -> float:
    return np.max(arr) #V


def calc_H(U_Hs: np.ndarray, I: float, N: int, l: float) -> np.ndarray:
    U_max: float = get_U_max(U_Hs)
    return U_Hs * I / U_max  * N / l  # A/m


def calc_M(U_Ms: np.ndarray, q:float, N: int, nu: float = 50.0) -> np.ndarray:
    return U_Ms / (47*N*q*nu*mu_0)  # A/m


def calc_H_Ms(Us: np.ndarray, I: float, N: int, q: float, nu: float = 50.0) -> np.ndarray:
    r: float = calc_r(q)
    l: float = calc_l(N, r)
    Hs = calc_H(Us[:, 0], I, N, l)
    Ms = calc_M(Us[:, 1], q, N, nu)
    H_Ms = np.column_stack((Hs, Ms))
    return H_Ms