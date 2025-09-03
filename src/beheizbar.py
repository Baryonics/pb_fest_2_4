from hysteresis import Hysteresis
from commutation_curve import CommutationCurve
from hysteresis import export_latex
from data_parser import DataParser as dp
import numpy as np
import magn_core as mc
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter
from uncertainties import ufloat


### Data Paths ###
path_Us_3A = "data/A_3-3-1/KH_3A"
path_Us_1A = "data/A_3-3-1/KH_1A"
path_Us_300mA = "data/A_3-3-1/KH_300mA"
path_Us_100mA = "data/A_3-3-1/KH_100mA"

path_Us_komm_3A = "data/A_3-3-2_Susz/KH_3A"
path_Us_komm_100mA = "data/A_3-3-2_Susz/KH_100mA"

path_Us_temp_3A = "data/A_3-3-3_Temp/KH_3A_Temp"


### Output Paths ###
plot_path = "res/beheizbar/plots/"
res_datapath = "res/beheizbar/"


### Parse Data ###
# Us_xA[zeile-, spalte| ], Us_xA[row-, column| ]
Us_1A: np.ndarray = dp.parse_to_np(path_Us_1A)
Us_3A: np.ndarray = dp.parse_to_np(path_Us_3A)
Us_300mA: np.ndarray = dp.parse_to_np(path_Us_300mA)
Us_100mA: np.ndarray = dp.parse_to_np(path_Us_100mA)

Us_komm_3A: np.ndarray = dp.parse_to_np(path_Us_komm_3A)
Us_komm_100mA: np.ndarray = dp.parse_to_np(path_Us_komm_100mA)

Us_temp_3A: np.ndarray = dp.parse_to_np(path_Us_temp_3A)

### Set Constants ###
q: float = 0.9*1e-4 #m^2
N: int = 17
r: float = mc.calc_r(q) #m
l: float = mc.calc_l(N, r) #m

### Excercice 3.3.1 ###
#---------------------------------------------------------------------#
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

### Save data to tex files ###
df_tex_3A = hys_3A.get_data_df_tex(3.0, y_scale=-4)
df_tex_1A = hys_1A.get_data_df_tex(1.0, y_scale=-4)
df_tex_300mA = hys_300mA.get_data_df_tex(0.3, y_scale=-4)
df_tex_100mA = hys_100mA.get_data_df_tex(0.1, y_scale=-4)

# Combine DataFrames
dfs = [df_tex_3A, df_tex_1A, df_tex_300mA, df_tex_100mA]



### Plot Hysteresis Curves ###
fig_1A = hys_1A.plot_hysteresis(
    title="Hysterese $I=1A$",
    xlabel="$H$ in $A/m$",
    ylabel="$M$ in $A/m$",
    rem_label="Remanenz $M_r$",
    coerc_label="Koerzitivfeldstärke $H_c$"
)

fig_3A = hys_3A.plot_hysteresis(
    title="Hysterese $I=3A$",
    xlabel="$H$ in $A/m$",
    ylabel="$M$ in $A/m$",
    rem_label="Remanenz $M_r$",
    coerc_label="Koerzitivfeldstärke $H_c$"
)

fig_300mA = hys_300mA.plot_hysteresis(
    title="Hysterese $I=300mA$",
    xlabel="$H$ in $A/m$",
    ylabel="$M$ in $A/m$",
    rem_label="Remanenz $M_r$",
    coerc_label="Koerzitivfeldstärke $H_c$"
)

fig_100mA = hys_100mA.plot_hysteresis(
    title="Hysterese $I=100mA$",
    xlabel="$H$ in $A/m$",
    ylabel="$M$ in $A/m$",
    rem_label="Remanenz $M_r$",
    coerc_label="Koerzitivfeldstärke $H_c$"
)


### Save plots ###
fig_1A.savefig(plot_path + "1A.png")
fig_3A.savefig(plot_path + "3A.png")
fig_300mA.savefig(plot_path + "300mA.png")
fig_100mA.savefig(plot_path + "100mA.png")


### Save Data to tex ###
# Save results 
export_latex(dfs,"Kenngrößen von Hysteresen bei verschiedenen Stromstärken $I$" ,res_datapath, "beheizbar_kenngr.tex")

# Save fit parameters to .tex
hys_1A.save_fit_params_to_tex(res_datapath, "1A_hyst_fit_param.tex", 1.0)
hys_3A.save_fit_params_to_tex(res_datapath, "3A_hyst_fit_param.tex", 3.0)
hys_300mA.save_fit_params_to_tex(res_datapath, "300mA_hyst_fit_param.tex", 0.3)
hys_100mA.save_fit_params_to_tex(res_datapath, "100mA_hyst_fit_param.tex", 0.1)
#---------------------------------------------------------------------#


### Excercise 3.3.2 ###
#---------------------------------------------------------------------#
### Calculate H adn M values ###
H_Ms_komm_3A =      mc.calc_H_Ms(Us_komm_3A, 3.0, N, q)
H_Ms_komm_100mA =   mc.calc_H_Ms(Us_komm_100mA, 0.1, N, q)


### Create Commutation Curve Objects ###
commcu_3A =     CommutationCurve(H_Ms_komm_3A)
commcu_100mA =  CommutationCurve(H_Ms_komm_100mA)


### Plot Commutation Curves ###
fig_commcu_fit_3A = commcu_3A.plot_commutation_curve(
    title="Kommutierungskurve bei $I_{max}=3A$",
    xlabel="$H$ in $A/m$",
    ylabel="$M$ in $A/m$"
)

fig_commcu_fit_100mA = commcu_100mA.plot_commutation_curve(
    title="Kommutierungskurve bei $I_{max}=100mA$",
    xlabel="$H$ in $A/m$",
    ylabel="$M$ in $A/M$"
)


### Plot suszeptibility Curves ###
fig_dM_dH_100mA = commcu_100mA.plot_deriv(
    title="Suszeptibilität bei $I=100mA$",
    xlabel="$H$ in $A/m$",
    ylabel="Suszeptibilität $dM/dH$"
)

fig_dM_dH_3A = commcu_3A.plot_deriv(
    title="Suszeptibilität bei $I=3A$",
    xlabel="$H$ in $A/m$",
    ylabel="Suszeptibilität $dM/dH$"
)


### save parameters to tex ###
commcu_100mA.save_fit_params_to_tex(
    path=res_datapath,
    fname="komm_fit_params_100mA.tex",
    current=0.1
)

commcu_3A.save_fit_params_to_tex(
    path=res_datapath,
    fname="komm_fit_params_3A.tex",
    current=3.0
)


### Save all Plots ###
fig_commcu_fit_3A.savefig(plot_path + "komm_3A.png")
fig_commcu_fit_100mA.savefig(plot_path + "komm_100mA.png")
fig_dM_dH_3A.savefig(plot_path + "dM_dH_3A.png")
fig_dM_dH_100mA.savefig(plot_path +"dM_dH_100mA.png")



### Excercise 3.3.3 ###
#---------------------------------------------------------------------#
Ms_temp_3A: np.ndarray = mc.calc_M(Us_temp_3A[:,1],q,N)
Ts_temp_3A: np.ndarray = Us_temp_3A[:,0]

fig_temp, ax_temp = plt.subplots()
ax_temp.set_title("Temperaturabhängigkeit der Magnetisierung")
ax_temp.set_xlabel("Temperatur $T$ in $°C$")
ax_temp.set_ylabel("Magnetisierung $M$ in $A/M$")
ax_temp.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
ax_temp.grid()
ax_temp.scatter(Ts_temp_3A, Ms_temp_3A, label="Messdaten", s=0.4)

fig_temp.savefig(res_datapath+"temp.png")

mask = Ms_temp_3A > 0
T = Ts_temp_3A[mask].astype(float)
M = Ms_temp_3A[mask].astype(float)

dT = np.median(np.diff(T))
M_filtered = savgol_filter(M, window_length=21, polyorder=1, deriv=1, delta=dT, mode="interp")
index = np.argmin(M_filtered)
Tc_val = T[index]

w = 20
Tc_err = np.std(T[index-w:index+w])/np.sqrt(2*w)
Tc = ufloat(Tc_val, Tc_err)

Tc_tex = f"$$ (T_c = {Tc_val:.2f} \pm {Tc_err:.2f})°C $$"
with open(res_datapath+"currie_temp.tex", "w", encoding="utf-8") as f:
    f.write(Tc_tex)
plt.show()