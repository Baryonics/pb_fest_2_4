import hysteresis
from data_parser import DataParser as dp
import numpy as np
import magn_core as mc
from matplotlib import pyplot as plt
import matplotlib

path_Us_3A = "data/A_3-3-1/KH_3A"
path_Us_1A = "data/A_3-3-1/KH_1A"
path_Us_300mA = "data/A_3-3-1/KH_300mA"
path_Us_100mA = "data/A_3-3-1/KH_100mA"

plot_path = "res/plots/beheizbar/"
res_datapath = "res/beheizbar/"

# Us_xA[zeile-, spalte| ], Us_xA[row-, column| ]
Us_1A: np.ndarray = dp.parse_to_df(path_Us_1A)
Us_3A: np.ndarray = dp.parse_to_df(path_Us_3A)
Us_300mA: np.ndarray = dp.parse_to_df(path_Us_300mA)
Us_100mA: np.ndarray = dp.parse_to_df(path_Us_100mA)

q: float = 0.9*1e-4 #m^2
N: int = 17
r: float = mc.calc_r(q) #m
l: float = mc.calc_l(N, r) #m


H_Ms_1A: np.ndarray = mc.calc_H_Ms(Us_1A, 1.0, N, q)
H_Ms_3A: np.ndarray = mc.calc_H_Ms(Us_3A, 3.0, N, q)
H_Ms_300mA: np.ndarray = mc.calc_H_Ms(Us_300mA, 0.3, N, q)
H_Ms_100mA: np.ndarray = mc.calc_H_Ms(Us_100mA, 0.1, N, q)

M = np.array([1,2,3])

hys_1A = hysteresis.Hysteresis(H_Ms_1A)
hys_3A = hysteresis.Hysteresis(H_Ms_3A)
hys_300mA = hysteresis.Hysteresis(H_Ms_300mA)
hys_100mA = hysteresis.Hysteresis(H_Ms_100mA)

fig_1A = hys_1A.plot_hysteresis("Hysterese 1A", "H in A/m", "M in A/m")
fig_3A = hys_3A.plot_hysteresis("Hysterese 3A", "H in A/m", "M in A/m")
fig_300mA = hys_300mA.plot_hysteresis("Hysterese 300mA", "H in A/m", "M in A/m")
fig_100mA = hys_100mA.plot_hysteresis("Hysterese 100mA", "H in A/m", "M in A/m")

print(hys_100mA.fit_function_lower(2))
plt.show()