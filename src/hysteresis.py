import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from fit_params_tex import FitParametersTex
import fit_params_tex as fptx


def export_latex(dfs, cap: str, path: str, filename: str,):
    combined_df = (
        pd.concat([d.set_index("$I$ in $[A]$") for d in dfs], axis=1) 
        .reset_index()
    )

    tex = combined_df.to_latex(
        index=False,
        escape=False,
        column_format="c|c|c|c|c",
        caption=cap,
        label="tab:hysterese",
        position="H"
    )

    tex = (tex.replace(r"\toprule", r"")
            .replace(r"\midrule", r"\hline")
            .replace(r"\bottomrule", r"")
            .replace(r"+/-", r" \pm "))
    
    with open(path+filename, "w", encoding="utf-8") as f:
        f.write(tex)



class Hysteresis:  
    def __init__(self, X_Ys: np.ndarray):
        self.X_Ys = X_Ys
        self.param_names = ["$L_0$", "$U_0$", "$k_0$", "$x_{c}^0$", "$\nu_0$", "$m_0$"]  
        self.b_pos, self.b_neg = self._split_branch()
        self.b_pos_params, self.b_pos_param_err = self._fit_branch(self.b_pos)
        self.b_neg_params, self.b_neg_param_err = self._fit_branch(self.b_neg)
        self.y_r, self.y_r_lower, self.y_r_upper = self._calc_y_r()
        self.x_c, self.x_c_lower, self.x_c_upper = self._calc_x_c()
        self.M_max, self.M_max_lower, self.M_max_upper = self._calc_M_max()
        
    
    
    def _split_branch_rude(self) -> tuple[np.ndarray, np.ndarray]:
        b_pos_list, b_neg_list = [], []

        for i in range(len(self.X_Ys) - 1):
            row = self.X_Ys[i, :]
            if self.X_Ys[i, 0] - self.X_Ys[i+1, 0] < 0:
                b_neg_list.append(row)
            else:
                b_pos_list.append(row)

        b_pos= np.array(b_pos_list).reshape(-1, 2)
        b_neg = np.array(b_neg_list).reshape(-1, 2)
        return b_pos, b_neg
    
    
    
    def _clean_branches(self, b_pos: np.ndarray, b_neg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        pass
    
    
    
    def _split_branch(self) -> tuple[np.ndarray, np.ndarray]:
        b_pos, b_neg = self._split_branch_rude()
        #b_pos, b_neg = self._cleanbranches(b_pos, b_neg)
        return b_pos, b_neg
    
    
    
    def _richards(self, x, L, U, k, x_c, nu, m):
        return L + (U - L) / (1.0 + np.exp(-k * (x - x_c)))**nu + m*x
    
    
    
    def _fit_branch(self, branch: np.ndarray):
        x = branch[:, 0].astype(float)
        y = branch[:, 1].astype(float)

        x_min, x_max = np.min(x), np.max(x)
        y_min, y_max = np.min(y), np.max(y)
        x_range = max(x_max - x_min, 1e-12)

        L0 = y_min
        U0 = y_max
        k0 = 2.0 / x_range
        x_c0 = np.median(x)
        nu0 = 1.0
        m0 = 0.0

        p0 = [L0, U0, k0, x_c0, nu0, m0]
             
        bounds = ([y_min, y_min, 0.0, x_min, 0.0],
                  [y_max, y_max, np.inf, x_max, np.inf, np.inf])

        popt, pcov = curve_fit(
            self._richards,
            x, y,
            p0=p0,
            #bounds=bounds,
            maxfev=20000
        )
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
  
    
    
    def fit_function_upper(self, x: float) -> float:
        return self._richards(x, *self.b_neg_params)
    
    
    
    def fit_function_lower(self, x: float) -> float:
        y = self._richards(x, *self.b_pos_params)
        return self._richards(x, *self.b_pos_params)
    
    
    
    def _calc_y_r(self):
        y_r_lower = self.fit_function_lower(0.0)
        y_r_upper = self.fit_function_upper(0.0)
        y_r_val = (abs(y_r_lower) + abs(y_r_upper)) / 2.0
        y_r_err = abs((abs(y_r_lower) - abs(y_r_upper)) / 2.0)
        return ufloat(y_r_val, y_r_err), y_r_lower, y_r_upper


    
    def _calc_x_c(self):
        sol_lower = root_scalar(self.fit_function_lower, bracket=[-500, 500], method='bisect')
        sol_upper = root_scalar(self.fit_function_upper, bracket=[-500, 500], method='bisect')
        x_c_lower = sol_lower.root
        x_c_upper = sol_upper.root
        x_c_val = (abs(x_c_lower) + abs(x_c_upper)) / 2.0
        x_c_err = abs((abs(x_c_lower) - abs(x_c_upper)) / 2.0)
        return ufloat(x_c_val, x_c_err), x_c_lower, x_c_upper
    
    
    
    def _calc_M_max(self):
        M_max_upper = max(self.X_Ys[:, 1])
        M_max_lower = min(self.X_Ys[:, 1])
        M_max_val = (abs(M_max_lower) + abs(M_max_upper)) / 2.0
        M_max_err = abs((abs(M_max_lower) - abs(M_max_upper)) / 2.0)
        return ufloat(M_max_val, M_max_err), M_max_lower, M_max_upper
    
    
    
    def plot_hysteresis(self, title: str, xlabel: str, ylabel: str, rem_label: str, coerc_label: str):
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
        ax.grid()
        ax.scatter(self.X_Ys[:, 0], self.X_Ys[:, 1], label="Messdaten", s=3)
        #ax.scatter(self.b_neg[:, 0], self.b_neg[:, 1], label="branch 2", s=3)
        
        x_fit = np.linspace(min(self.X_Ys[:,0]), max(self.X_Ys[:,0]), 1000)
        y_fit_b_pos = self._richards(x_fit, *self.b_pos_params)
        y_fit_b_neg = self._richards(x_fit, *self.b_neg_params)
        
        ax.plot(x_fit, y_fit_b_pos, 'r', label="Fit vorlaufender Ast")
        ax.plot(x_fit, y_fit_b_neg, 'g', label="Fit rücklaufender Ast")
        
        #ax.scatter(self.x_c_lower, 0, color='yellow', s=60, marker='D', edgecolors='black', label="x_c lower")
        #ax.scatter(self.x_c_upper, 0, color='green', s=60, marker='P', edgecolors='black', label="x_c upper")
        #ax.scatter(0, self.y_r_lower, color='magenta', s=60, marker='X', edgecolors='black', label="y_r lower")
        #ax.scatter(0, self.y_r_upper, color='cyan', s=60, marker='*', edgecolors='black', label="y_r upper")
        
        ax.errorbar(0, self.y_r.n, yerr=self.y_r.s, fmt='o', color='cyan', label=rem_label)
        ax.errorbar(self.x_c.n, 0, xerr=self.x_c.s, fmt='o', color='magenta', label=coerc_label)
        ax.legend()
        return fig



    def save_fit_params_to_tex(self, path: str, fname: str, current: float):
        fit_tex_pos = FitParametersTex(
            popt=self.b_pos_params,
            perr=self.b_pos_param_err,
            names=self.param_names,
            lab=f"tab:fit_par_pos_{current}A",
            cap=f"Fit Parameter Hysterese bei $I={current}A$ vorläufige Richtung "
        )
        
        fit_tex_neg = FitParametersTex(
            popt=self.b_neg_params,
            perr=self.b_neg_param_err,
            names=self.param_names,
            lab=f"tab:fit_par_neg_{current}A",
            cap=f"Fit Parameter Hysterese bei $I={current}A$ rückläufige Richtung "
        )
        
        combined = fptx.combine(
            param_a=fit_tex_pos,
            param_b=fit_tex_neg,
            cap=f"Fit Parameter für den Vor- und rücklaufenden Strang bei $I={current}A$",
            lab=f"tab:fit_param_{current}A"
        )

        fptx.params_export_tex(combined, path, fname)
    
    
    def get_data_df_tex(self, I: float, x_scale: float = 1, y_scale: float = 1) -> pd.DataFrame:
        df_tex = pd.DataFrame({
            "$I$ in $[A]$":     [
                                 f"$M_r^+$ in $[A/m]\cdot 10^{{{y_scale*-1}}}$",
                                 f"$M_r^-$ in $[A/m]\cdot 10^{{{y_scale*-1}}}$",
                                 f"$M_r$ in $[A/m]\cdot 10^{{{y_scale*-1}}}$",
                                 f"$H_c^+$ in $[A/m]$",
                                 f"$H_c^-$ in $[A/m]$",
                                 f"$H_c$ in $[A/m]$",
                                 f"$M_{{max}}^+$ in $[A/m]\cdot 10^{{{y_scale*-1}}}$",
                                 f"$M_{{max}}^-$ in $[A/m]\cdot 10^{{{y_scale*-1}}}$",
                                 f"$M_{{max}}$ in $[A/m]\cdot 10^{{{y_scale*-1}}}$"
                                 ],
            f"${I}$":           [
                                 f"${self.y_r_upper*10**y_scale:.2f}$",
                                 f"${self.y_r_lower*10**y_scale:.2f}$",
                                 f"${self.y_r*10**y_scale:.2f}$",
                                 f"${self.x_c_upper*x_scale:.2f}$",
                                 f"${self.x_c_lower*x_scale:.2f}$",
                                 f"${self.x_c*x_scale:.2f}$",
                                 f"${self.M_max_upper*10**y_scale:.2f}$",
                                 f"${self.M_max_lower*10**y_scale:.2f}$",
                                 f"${self.M_max*10**y_scale:.2f}$"
                                ]
        })
        return df_tex