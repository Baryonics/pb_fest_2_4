import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from fit_params_tex import FitParametersTex



class CommutationCurve:
    def __init__(self, X_Ys: np.ndarray):
        self.xs = X_Ys[:,0]
        self.ys = X_Ys[:,1]
        self.komm_fit_used, self.komm_fit_param, self.komm_fit_err = self.fit_data()
        print(self.komm_fit_used)
        print(self.komm_param_names)
        print(self.komm_fit_param)
    
    
    
    def fit_func(self, x, y_s, x_c, a, y_off):
        return y_s * np.tanh((x - x_c) * a ) + y_off
    
    
    
    def fit_func_2(self, x, y_s, a):
        return y_s * (1/np.tanh(a*x) - 1/(a*x))
    
    
    
    def fit_data(self):
        try:
            popt, pcov = curve_fit(self.fit_func, self.xs, self.ys, p0=[np.sign(max(self.ys))*max(self.ys), 0, 0.001, 0])
            self.ys_fit = self.fit_func(self.xs, *popt)
            self.komm_param_names = ["$y_s$", "$x_c$", "$a$", "$y_{off}$"]
            f_used = self.fit_func
        except Exception:
            popt, pcov = curve_fit(self.fit_func_2, self.xs, self.ys, p0=[max(self.ys), 0.01])
            self.ys_fit = self.fit_func_2(self.xs, *popt)
            self.komm_param_names = ["$y_s$", "$a$"]
            f_used = self.fit_func_2
        perr = np.sqrt(np.diag(pcov))
        return f_used, popt, perr
    
    

    def plot_commutation_curve(self, title: str, xlabel: str, ylabel: str):
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
        ax.grid()
        ax.scatter(self.xs, self.ys, label="Messdaten", s=3)
        ax.plot(self.xs, self.ys_fit, label="Fit", color="red")
        ax.legend()
        return fig
    
    
    
    def plot_deriv(self, title: str, xlabel: str, ylabel: str ):
        xs_ = np.linspace(0, max(self.xs), 10000)
        ys_ = self.komm_fit_used(xs_, *self.komm_fit_param)
        dy_dx = np.gradient(ys_, xs_)
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
        ax.grid()
        ax.plot(xs_, dy_dx, label="$dM/dH$", color="red")
        ax.legend()
        return fig
    
    
    
    def save_fit_params_to_tex(self, path: str, fname: str, current: float):
        fit_tex_pos = FitParametersTex(
            popt=self.komm_fit_param,
            perr=self.komm_fit_err,
            names=self.komm_param_names,
            lab=f"tab:komm_par_{current}A",
            cap=f"Fit Parameter Kommutierungskurve bei $I={current}A$."
        )
        
        fit_tex_pos.params_export_tex(path, fname)