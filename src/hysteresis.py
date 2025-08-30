import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


class Hysteresis:
    def __init__(self, X_Ys: np.ndarray):
        self.X_Ys = X_Ys
        self.b1, self.b2 = self._split_branch()
        self.b1_params, self.b1_param_err = self._fit_branch(self.b1)
        self.b2_params, self.b2_param_err = self._fit_branch(self.b2)
    
    
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
        #b1, b2 = self._clean_branches(b1, b2)
        return b1, b2
    
    
    def fit_func(self, x: float, x_c: float, y_s: float, a: float, y_off: float) -> float:
        return y_s * np.tanh((x - x_c) * a ) + y_off
    
    
    def _fit_branch(self, branch: np.ndarray):
        x = branch[:, 0].astype(float)
        y = branch[:, 1].astype(float)

        x_min, x_max = np.min(x), np.max(x)
        y_min, y_max = np.min(y), np.max(y)
        x_range = max(x_max - x_min, 1e-12)

        # sinnvolle Startwerte
        x_c0   = np.median(x)                           # Zentrum der S-Kurve
        y_s0   = 0.5 * (y_max - y_min)                  # Sättigung
        y_off0 = 0.5 * (y_max + y_min)                  # Offset
        a0     = 2.0 / x_range                          # Steilheit ~ 1/Breitenmaß

        p0 = [x_c0, y_s0, a0, y_off0]

        # Bounds: x_c im Datenbereich, y_s >= 0, a > 0
        bounds = ([x_min,     0.0,   0.0, -np.inf],
                [x_max,  np.inf, np.inf,  np.inf])

        popt, pcov = curve_fit(
            self.fit_func,
            x, y,
            p0=p0,
            bounds=bounds,
            maxfev=20000
        )
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    
    def plot_hysteresis(self, title: str, xlabel: str, ylabel: str):
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
        ax.scatter(self.b1[:, 0], self.b1[:, 1], label="branch 1", s=3)
        ax.scatter(self.b2[:, 0], self.b2[:, 1], label="branch 2", s=3)
        
        x_fit = np.linspace(min(self.X_Ys[:,0]), max(self.X_Ys[:,0]), 1000)
        y_fit_b1 = self.fit_func(x_fit, *self.b1_params)
        y_fit_b2 = self.fit_func(x_fit, *self.b2_params)
        
        ax.plot(x_fit, y_fit_b1, 'r', label="fit branch 1")
        ax.plot(x_fit, y_fit_b2, 'g', label="fit branch 2")
        
        ax.legend()
        return fig