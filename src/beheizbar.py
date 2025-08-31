from hysteresis import Hysteresis
from hysteresis import export_latex
from data_parser import DataParser as dp
import numpy as np
import magn_core as mc
from matplotlib import pyplot as plt
import matplotlib
from uncertainties import ufloat
import pandas as pd

### Data Paths ###
path_Us_3A = "data/A_3-3-1/KH_3A"
path_Us_1A = "data/A_3-3-1/KH_1A"
path_Us_300mA = "data/A_3-3-1/KH_300mA"
path_Us_100mA = "data/A_3-3-1/KH_100mA"


### Output Paths ###
plot_path = "res/plots/beheizbar/"
res_datapath = "res/beheizbar/"


### Parse Data ###
# Us_xA[zeile-, spalte| ], Us_xA[row-, column| ]
Us_1A: np.ndarray = dp.parse_to_np(path_Us_1A)
Us_3A: np.ndarray = dp.parse_to_np(path_Us_3A)
Us_300mA: np.ndarray = dp.parse_to_np(path_Us_300mA)
Us_100mA: np.ndarray = dp.parse_to_np(path_Us_100mA)


### Set Constants ###
q: float = 0.9*1e-4 #m^2
N: int = 17
r: float = mc.calc_r(q) #m
l: float = mc.calc_l(N, r) #m


### Calculate H and M ###
H_Ms_1A: np.ndarray = mc.calc_H_Ms(Us_1A, 1.0, N, q)
H_Ms_3A: np.ndarray = mc.calc_H_Ms(Us_3A, 3.0, N, q)
H_Ms_300mA: np.ndarray = mc.calc_H_Ms(Us_300mA, 0.3, N, q)
H_Ms_100mA: np.ndarray = mc.calc_H_Ms(Us_100mA, 0.1, N, q)


### Generate Hysteresis Objects ###
hys_1A = Hysteresis(H_Ms_1A)
hys_3A = Hysteresis(H_Ms_3A)
hys_300mA = Hysteresis(H_Ms_300mA)
hys_100mA = Hysteresis(H_Ms_100mA)

df_tex_3A = hys_3A.get_data_df_tex(3.0, y_scale=-3)
df_tex_1A = hys_1A.get_data_df_tex(1.0, y_scale=-3)
df_tex_300mA = hys_300mA.get_data_df_tex(0.3, y_scale=-4)
df_tex_100mA = hys_100mA.get_data_df_tex(0.1, y_scale=-4)

print(df_tex_3A)

dfs = [df_tex_3A, df_tex_1A, df_tex_300mA, df_tex_100mA]


export_latex(dfs, res_datapath, "beheizbar_results.tex")


### Plot Hysteresis Curves ###
fig_1A = hys_1A.plot_hysteresis("Hysterese $I=1A$", "$H$ in $A/m$", "M in A/m", "Remanenz $M_r$", "Koerzitivfeldstärke $H_c$")
fig_3A = hys_3A.plot_hysteresis("Hysterese $I=3A$", "$H$ in $A/m$", "M in A/m", "Remanenz $M_r$", "Koerzitivfeldstärke $H_c$")
fig_300mA = hys_300mA.plot_hysteresis("Hysterese $I=300mA$", "$H$ in $A/m$", "M in A/m", "Remanenz $M_r$", "Koerzitivfeldstärke $H_c$")
fig_100mA = hys_100mA.plot_hysteresis("Hysterese $I=100mA$", "$H in A/m$", "M in A/m", "Remanenz $M_r$", "Koerzitivfeldstärke $H_c$")


plt.show()