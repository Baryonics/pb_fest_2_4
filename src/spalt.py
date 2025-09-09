from hysteresis import Hysteresis
from commutation_curve import CommutationCurve
from hysteresis import export_latex
from data_parser import DataParser as dp
import numpy as np
import magn_core as mc
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter
from uncertainties import ufloat
import pandas as pd


### Data Paths ###
path_Us_0mm_940mA:  str         =   "data/A_3-4-1_Kengr/KS_940mA_0mm"
path_Us_0mm_3A:     str         =   "data/A_3-4-2_Entmag/KS_3A_0mm"
path_Us_0_075mm_3A: str         =   "data/A_3-4-2_Entmag/KS_3A_0-075mm"
path_Us_0_125mm_3A: str         =   "data/A_3-4-2_Entmag/KS_3A_0-125mm"
path_Us_0_2mm_3A:   str         =   "data/A_3-4-2_Entmag/KS_3A_0-2mm"
path_Us_0_5mm_3A:   str         =   "data/A_3-4-2_Entmag/KS_3A_0-5mm"
path_Us_1mm_3A:     str         =   "data/A_3-4-2_Entmag/KS_3A_1mm"


### Output Paths ###
plot_path:          str         =   "res/spalt/plots/"
res_datapath:       str         =   "res/spalt/"


### Parse Data ###
Us_940mA:           np.ndarray  =   dp.parse_to_np(path_Us_0mm_940mA)
Us_0mm:             np.ndarray  =   dp.parse_to_np(path_Us_0mm_3A)
Us_0_075mm:         np.ndarray  =   dp.parse_to_np(path_Us_0_075mm_3A)
Us_0_125mm:         np.ndarray  =   dp.parse_to_np(path_Us_0_125mm_3A)
Us_0_2mm:           np.ndarray  =   dp.parse_to_np(path_Us_0_2mm_3A)
Us_0_5mm:           np.ndarray  =   dp.parse_to_np(path_Us_0_5mm_3A)
Us_1mm:             np.ndarray  =   dp.parse_to_np(path_Us_1mm_3A)


### Set Constants ###
q:          float   = 9.0*1e-5 #m^2
N_p:        int     = 54
N_s:        int     = 17
r:          float   = 0.015 #m
l:          float   = mc.calc_l(N_p, r) #m
I:          float   = 3.0 #A
y_scale:    float   = -4.0

### Excercice 3.4.1 ###
#---------------------------------------------------------------------#
### Calculate H and M ###
H_Ms_940mA:         np.ndarray  =   mc.calc_H_Ms(Us_940mA, 0.940, N_p, q, N_s=N_s)
H_Ms_0mm:           np.ndarray  =   mc.calc_H_Ms(Us_0mm, I, N_p, q, N_s=N_s)
H_Ms_0_075mm:       np.ndarray  =   mc.calc_H_Ms(Us_0_075mm, I, N_p, q, N_s=N_s)
H_Ms_0_125mm:       np.ndarray  =   mc.calc_H_Ms(Us_0_125mm, I, N_p, q, N_s=N_s)
H_Ms_0_2mm:         np.ndarray  =   mc.calc_H_Ms(Us_0_2mm, I, N_p, q, N_s=N_s)
H_Ms_0_5mm:         np.ndarray  =   mc.calc_H_Ms(Us_0_5mm, I, N_p, q, N_s=N_s)
H_Ms_1mm:           np.ndarray  =   mc.calc_H_Ms(Us_1mm, I, N_p, q, N_s=N_s)   


### Generate Hysteresis Objects ###
hys_940mA_0mm:       Hysteresis  =   Hysteresis(H_Ms_940mA)
hys_0mm:            Hysteresis  =   Hysteresis(H_Ms_0mm)
hys_0_075mm:        Hysteresis  =   Hysteresis(H_Ms_0_075mm, is_fit=False)
hys_0_125mm:        Hysteresis  =   Hysteresis(H_Ms_0_125mm, is_fit=False)
hys_0_2mm:          Hysteresis  =   Hysteresis(H_Ms_0_2mm, is_fit=False)
hys_0_5mm:          Hysteresis  =   Hysteresis(H_Ms_0_5mm, is_fit=False)
hys_1mm:            Hysteresis  =   Hysteresis(H_Ms_1mm, is_fit=False)   


### Save data to tex files ##   
df_tex_940mA_0mm:    pd.DataFrame=   hys_940mA_0mm.get_data_df_tex(0.940, y_scale=y_scale)
df_tex_0mm:          pd.DataFrame=   hys_0mm.get_data_df_tex(I, y_scale=y_scale)
#df_tex_0_075mm:     pd.DataFrame=   hys_0_075mm.get_data_df_tex(I, y_scale=y_scale)
#df_tex_0_125mm:     pd.DataFrame=   hys_0_125mm.get_data_df_tex(I, y_scale=y_scale)
#df_tex_0_2mm:       pd.DataFrame=   hys_0_2mm.get_data_df_tex(I,y_scale=y_scale)
#df_tex_0_5mm:       pd.DataFrame=   hys_0_5mm.get_data_df_tex(I, y_scale=y_scale)
#df_tex_1mm:         pd.DataFrame=   hys_1mm.get_data_df_tex(I, y_scale=y_scale)
export_latex(
    dfs     =   [df_tex_940mA_0mm],
    cap     =   "Hysteresekurve bei $I=940mA$ und $0mm$ Spaltbreite",
    path    =   res_datapath,
    filename=   "hyst_940mA_0mm.tex"
)
export_latex(
    dfs     =   [df_tex_0mm],
    cap     =   "Hysteresekurve bei $I=3.0A$ und $0mm$ Spaltbreite",
    path    =   res_datapath,
    filename=   "hyst_3A_0mm.tex"
)


fig_940mA_0mm,_ = hys_940mA_0mm.plot_hysteresis(
    title="Hysterese $I=940mA$ und Spaltbreite $d=0mm$",
    xlabel="$H$ in $A/m$",
    ylabel="$M$ in $A/m$",
    rem_label="Remanenz $M_r$",
    coerc_label="Koerzitivfeldstärke $H_c$"
)

fig_hysts, ax_hysts     =   plt.subplots()

ax_hysts.set_title("Hysterese bei $I=3A$ mit verschiedenen Spaltbreiten")
### Plot Hysteresis Curve ###
fig_0mm,ax_0mm          =   hys_0mm.plot_hysteresis(
    label               =   "$0mm$",
    ax_                 =   ax_hysts,
    is_plot_fit         =   False
)

fig_0_075mm,ax_0_075mm  =   hys_0_075mm.plot_hysteresis(
    label               =   "$0.075mm$",
    ax_                 =   ax_hysts
)

fig_125mm,ax_0_125mm    =   hys_0_125mm.plot_hysteresis(
    label               =   "$0.125mm$",      
    ax_                 =   ax_hysts
)

fig_0_2mm,ax_0_2mm      =   hys_0_2mm.plot_hysteresis(
    label               =   "$0.2mm$",
    ax_                 =   ax_hysts
)

fig_0_5mm,ax_0_5mm      =   hys_0_5mm.plot_hysteresis(
    label               =   "$0.5mm$",
    ax_                 =   ax_hysts
)

fig_1mm,ax_1mm          =   hys_1mm.plot_hysteresis(
    label               =   "1.0mm",
    ax_                 =   ax_hysts
)



fig_hysts.savefig(plot_path + "spalt_hysteresen.png")
fig_940mA_0mm.savefig(plot_path + "940mA_0mm.png")

hys_940mA_0mm.save_fit_params_to_tex(
    path                =   res_datapath,
    fname               =   "param_spalt_hysteresis_940mA_0_mm.tex",
    current             =   0.940
)


hys_0mm.save_fit_params_to_tex(
    path                =   res_datapath,
    fname               =   "param_spalt_hysteresis_0_mm.tex",
    current             =   3.0,
)

def H_ent(H_0, H_sp):
    return H_sp - H_0


def N(H_0, H_sp, M_max):
    return H_ent(H_0, H_sp) / M_max

H_0 = hys_0mm.H_at_M_max

def hys_to_N(hys: Hysteresis):
    return N(H_0, hys.H_at_M_max.n, hys.M_max.n)

M_0 = hys_1mm.M_max.n
print("M_0 = ", M_0)

mask = hys_0mm.X_Ys[:,0] > 0
X_pos = hys_0mm.X_Ys[mask]


idx = np.argmin(np.abs(X_pos[:,1] - M_0))
x_val = X_pos[idx, 0]
y_val = X_pos[idx, 1]

print(f"x = {x_val}, y = {y_val}")

Ns = [hys_to_N(hys) for hys in [hys_0mm, hys_0_075mm, hys_0_125mm, hys_0_2mm, hys_0_5mm, hys_1mm]]
Hs = [H for H in [hys_0_075mm.H_at_M_max.n, hys_0_125mm.H_at_M_max.n, hys_0_2mm.H_at_M_max.n, hys_0_5mm.H_at_M_max.n, hys_1mm.H_at_M_max.n]]
#ds = [0.0, 0.075,0.125, 0.2, 0.5, 1.0]
ds = np.linspace(0.0, 1000.0, 1000)
N_theos = [2*d*1e-3/(2*np.pi*r+2*d*1e-3) for d in ds]

fig_theos, ax_theos = plt.subplots()

ax_theos.plot(ds, N_theos)
ax_theos.set_title("$N_{theo}$ in Abhängigkeit von $d$")
ax_theos.set_xlabel("$d$ in $[mm]$")
ax_theos.set_ylabel("$N$")
ax_theos.grid()

fig_theos.savefig(plot_path + "N_theo.png")

plt.show()