import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt


def export_latex(dfs,  path: str, filename: str,):
    #combined_df = pd.concat(dfs, ignore_index=True)
    combined_df = (
        pd.concat([d.set_index("$I$ in $[A]$") for d in dfs], axis=1) 
        .reset_index()
    )

    tex = combined_df.to_latex(
        index=False,
        escape=False,
        column_format="lrrrrrr",
        caption="Hysterese–Kennwerte",
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
        self.b1, self.b2 = self._split_branch()
        self.b1_params, self.b1_param_err = self._fit_branch(self.b1)
        self.b2_params, self.b2_param_err = self._fit_branch(self.b2)
        self.y_r, self.y_r_lower, self.y_r_upper = self._calc_y_r()
        self.x_c, self.x_c_lower, self.x_c_upper = self._calc_x_c()
        self.M_max, self.M_max_lower, self.M_max_upper = self._calc_M_max()
        
    
    
    def _split_branch_rude(self) -> tuple[np.ndarray, np.ndarray]:
        b1_list, b2_list = [], []

        for i in range(len(self.X_Ys) - 1):
            row = self.X_Ys[i, :]
            if self.X_Ys[i, 0] - self.X_Ys[i+1, 0] < 0:
                b1_list.append(row)
            else:
                b2_list.append(row)

        b1 = np.array(b1_list).reshape(-1, 2)
        b2 = np.array(b2_list).reshape(-1, 2)
        return b1, b2
    
    
    
    def _clean_branches(self, b1: np.ndarray, b2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        pass
    
    
    
    def _split_branch(self) -> tuple[np.ndarray, np.ndarray]:
        b1, b2 = self._split_branch_rude()
        #b1, b2 = self._cleanbranches(b1, b2)
        return b1, b2
    
    
    
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
        return self._richards(x, *self.b2_params)
    
    
    
    def fit_function_lower(self, x: float) -> float:
        y = self._richards(x, *self.b1_params)
        return self._richards(x, *self.b1_params)
    
    
    
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
        ax.grid()
        ax.scatter(self.X_Ys[:, 0], self.X_Ys[:, 1], label="Messdaten", s=3)
        #ax.scatter(self.b2[:, 0], self.b2[:, 1], label="branch 2", s=3)
        
        x_fit = np.linspace(min(self.X_Ys[:,0]), max(self.X_Ys[:,0]), 1000)
        y_fit_b1 = self._richards(x_fit, *self.b1_params)
        y_fit_b2 = self._richards(x_fit, *self.b2_params)
        
        ax.plot(x_fit, y_fit_b1, 'r', label="Fit unterer Ast")
        ax.plot(x_fit, y_fit_b2, 'g', label="Fit oberer Ast")
        
        #ax.scatter(self.x_c_lower, 0, color='yellow', s=60, marker='D', edgecolors='black', label="x_c lower")
        #ax.scatter(self.x_c_upper, 0, color='green', s=60, marker='P', edgecolors='black', label="x_c upper")
        #ax.scatter(0, self.y_r_lower, color='magenta', s=60, marker='X', edgecolors='black', label="y_r lower")
        #ax.scatter(0, self.y_r_upper, color='cyan', s=60, marker='*', edgecolors='black', label="y_r upper")
        
        ax.errorbar(0, self.y_r.n, yerr=self.y_r.s, fmt='o', color='cyan', label=rem_label)
        ax.errorbar(self.x_c.n, 0, xerr=self.x_c.s, fmt='o', color='magenta', label=coerc_label)
        ax.legend()
        return fig



    def get_data_df_tex(self, I: float, x_scale: float = 1, y_scale: float = 1) -> pd.DataFrame:
        #df_tex = pd.DataFrame({
            #"$I$ in $[A]$":             [f"${I:.2f}$"],
            #f"$M_r^+$ in $[A/m]$":      [f"${self.y_r_upper*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$"],
            #f"$M_r^-$ in $[A/m]$":      [f"${self.y_r_lower*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$"],
            #f"$M_r$ in $[A/m]$":        [f"${self.y_r*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$"],
            #f"$H_c^+$ in $[A/m]$":      [f"${self.x_c_upper*x_scale:.2f}$"],
            #f"$H_c^-$ in $[A/m]$":      [f"${self.x_c_lower*x_scale:.2f}$"],
            #f"$H_c$ in $[A/m]$":        [f"${self.x_c*x_scale:.2f}$"],
            #f"$Max^+$ in $[A/m]$":      [f"${self.M_max_upper*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$"],
            #f"$Max^-$ in $[A/m]$":      [f"${self.M_max_lower*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$"],
            #f"$M_{{max}}$ in $[A/m]$":  [f"${self.M_max*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$"]
            
        #})
        
    
        df_tex = pd.DataFrame({
            "$I$ in $[A]$":     [
                                 f"$M_r^+$ in $[A/m]$",
                                 f"$M_r^-$ in $[A/m]$",
                                 f"$M_r$ in $[A/m]$",
                                 f"$H_c^+$ in $[A/m]$",
                                 f"$H_c^-$ in $[A/m]$",
                                 f"$H_c$ in $[A/m]$",
                                 f"$Max^+$ in $[A/m]$",
                                 f"$Max^-$ in $[A/m]$",
                                 f"$M_{{max}}$ in $[A/m]$"
                                 ],
            f"${I}$":           [
                                 f"${self.y_r_upper*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$",
                                 f"${self.y_r_lower*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$",
                                 f"${self.y_r*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$",
                                 f"${self.x_c_upper*x_scale:.2f}$",
                                 f"${self.x_c_lower*x_scale:.2f}$",
                                 f"${self.x_c*x_scale:.2f}$",
                                 f"${self.M_max_upper*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$",
                                 f"${self.M_max_lower*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$",
                                 f"${self.M_max*10**y_scale:.2f}\cdot 10^{{{y_scale}}}$"
                                ]
        })
        return df_tex